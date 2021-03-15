#' Add adducts to the output of MSnbase::calculateFragments
#'
#' @param x `data.frame`, output of [MSnbase::calculateFragments()]
#' @param adducts `data.frame`, with 3 columns mass, name, to
#' @return `data.frame`
#' @noRd
.addAdducts <- function(x, adducts) {
    if (!nrow(adducts)) {
        return(x)
    }

    if (!all(c("mass", "name", "to") %in% colnames(adducts))) {
        stop(
            "The 'adducts' data.frame must have the columns: ",
            "'mass', 'name' and 'to'."
        )
    }

    r <- do.call(rbind, lapply(seq_len(nrow(adducts)), function(i) {
        a <- x[x$type == adducts$to[i], , drop = FALSE]
        a$mz <- a$mz + adducts$mass[i]
        a$ion <- paste0(adducts$name[i], a$pos)
        a
    }))
    x <- rbind(x, r)
    rownames(x) <- NULL
    x
}

#' Add adducts to the output of MSnbase::calculateFragments
#'
#' @param x `data.frame`, output of [MSnbase::calculateFragments()].
#' @param n `numeric(1)`, length of sequence.
#' @param modifications `data.frame`, with 4 columns mass, name, location,
#' variable.
#' @return `data.frame`
#' @noRd
.addCustomModifications <- function(x, n, modifications) {
    if (!nrow(modifications)) {
        return(x)
    }

    if (!all(c("mass", "name", "location", "variable") %in%
             colnames(modifications))) {
        stop(
            "The 'customModifications' data.frame must have the columns: ",
            "'mass', 'name', 'location' and 'variable'."
        )
    }

    if (is.character(modifications$location)) {
        modifications$location[tolower(modifications$location) == "n-term"] <-
            0L
        modifications$location[tolower(modifications$location) == "c-term"] <-
            n + 1L
        modifications$location <- as.integer(modifications$location)
    }

    for (i in seq_len(nrow(modifications))) {
        x <- .cusmod(
            x, n=n,
            id=modifications$name[i],
            deltamz=modifications$mass[i],
            location=modifications$location[i],
            variable=modifications$variable[i]
        )
    }
    x
}

#' Calculate Fragments (via MSnbase::calculateFragments)
#'
#' @param sequence `AAString`, peptide sequence
#' @param type `character`, type of fragments
#' @param modification `character`, unimod names
#' @param adducts `data.frame`, with 3 columns mass, name, to
#' @param neutralLoss `list`, neutral loss (see [MSnbase::calculateFragments()])
#' @param sequenceOrder `character`, `c("original", "random", "inverse")`
#' @param verbose `logical`, verbose output?
#' @return `FragmentViews`
#' @noRd
.calculateFragments <- function(sequence, type=c("a", "b", "c", "x", "y", "z"),
                                modifications=c(
                                    "Carbamidomethyl", "Acetyl",
                                    "Met-loss"
                                ),
                                customModifications=data.frame(),
                                adducts=data.frame(),
                                neutralLoss=defaultNeutralLoss(),
                                sequenceOrder=c(
                                    "original",
                                    "random",
                                    "inverse"
                                ),
                                verbose=interactive()) {
    modifications <- match.arg(
        modifications,
        choices=c(
        "", # allow NULL for nothing
        "Carbamidomethyl",
        "Acetyl",
        "Met-loss"
        ),
        several.ok=TRUE
    )
    sequenceOrder <- match.arg(sequenceOrder)

    ## just to be sure if an AAString is given
    csequence <- as.character(sequence)

    ## TODO: replace by unimod package
    if ("Met-loss" %in% modifications) {
        csequence <- .unimod765(csequence)
    }

    ## has to be done after Met-loss but before any other modification
    csequence <- .reorderSequence(csequence, method=sequenceOrder)

    d <- calculateFragments(
        csequence,
        type=type,
        modifications=NULL,
        neutralLoss=neutralLoss,
        verbose=FALSE
    )

    ## add protein sequence to data.frame to calculate modifications
    d <- rbind(
        list(
            mz=.calculateProteinMass(csequence),
            ion="none", type="protein", pos=0, z=0,
            seq=csequence
        ),
        d
    )

    ## TODO: replace by unimod package
    if ("Acetyl" %in% modifications) {
        d <- .unimod1(d, csequence)
    }
    ## TODO: replace by unimod package
    if ("Carbamidomethyl" %in% modifications) {
        d <- .unimod4(d)
    }
    if (all(!nchar(modifications))) {
        modifications <- NULL
    }
    ## TODO: replace by unimod package
    d <- .addCustomModifications(d, nchar(csequence), customModifications)
    d <- .addAdducts(d, adducts)

    ## remove protein sequence from data.frame (just added to calculate
    ## modifications)
    mass <- d$mz[1L]
    d <- d[-1L, ]

    n <- nchar(csequence)
    FragmentViews(
        csequence, mass=d$mz, type=d$type, z=Rle(d$z), names=d$ion,
        start=ifelse(
            startsWith(csequence, d$seq),
            1L,
            ifelse(endsWith(csequence, d$seq), n-d$pos + 1L, d$pos)
        ),
        width=nchar(d$seq),
        metadata=list(modifications=modifications, mass=mass)
    )
}

#' Calculate protein mass
#'
#' TODO: replace by unimod package
#'
#' @param x `character`, sequence
#' @return `numeric`
#' @noRd
.calculateProteinMass <- function(x) {
    aa <- get.amino.acids()[, c("AA", "ResidueMass")]
    aamass <- setNames(aa$ResidueMass, aa$AA)
    x <- strsplit(x, "")[[1L]]
    sum(aamass[x])
}

#' Match fragments and measured mz values.
#'
#' Similar to MALDIquant:::.match.closest but we need to reimplement it here,
#' because @pavel_shliaha asks for handling duplicated matches differently: see
#' https://github.com/sgibb/topdownr/issues/72
#'
#' @param mz `double`, measured mz.
#' @param fmass `double`, fragment mass
#' @param tolerance `double`, allowed tolerance
#' @param redundantIonMatch `character`, a mz could be matched to two or more
#' fragments, it would be removed or matched to the closest fragment.
#' @param redundantFragmentMatch `character`, multiple mz could be matched to
#' the same fragment, these matches would be removed or the closest mz is
#' chosen.
#' @param relative `logical`, relative tolerance?
#' @return `integer`
#' @noRd
.matchFragments <- function(mz, fmass, tolerance=5e-6,
                            redundantIonMatch=c("remove", "closest"),
                            redundantFragmentMatch=c("remove", "closest",
                                                     "ignore"),
                            relative=TRUE) {
    if (!length(mz)) {
        return(integer())
    }
    if (!length(fmass)) {
        return(rep.int(NA_integer_, length(mz)))
    }

    lIdx <- findInterval(mz, fmass, rightmost.closed=FALSE, all.inside=TRUE)
    rIdx <- lIdx + 1L
    lIdx[lIdx == 0L] <- 1L
    lDiff <- abs(fmass[lIdx] - mz)
    rDiff <- abs(fmass[rIdx] - mz)
    tolerance <- rep_len(tolerance, length(fmass))
    if (relative) {
        tolerance <- tolerance * fmass
    }
    lDiff[lDiff > tolerance[lIdx]] <- Inf
    rDiff[rDiff > tolerance[rIdx]] <- Inf

    d <- which(lDiff >= rDiff)
    lIdx[d] <- rIdx[d]

    ## no match at all
    lIdx[is.infinite(lDiff) & is.infinite(rDiff)] <- NA_integer_

    ## multiple mz to one fragment?
    if (anyDuplicated(lIdx)) {
        redundantFragmentMatch <- match.arg(redundantFragmentMatch)
        if (redundantFragmentMatch == "remove") {
            lIdx[duplicated(lIdx) | duplicated(lIdx, fromLast=TRUE)] <- NA_integer_
        } else if (redundantFragmentMatch == "closest") {
            o <- order(abs(mz - fmass[lIdx]))
            m <- lIdx[o]
            m[duplicated(m)] <- NA_integer_
            lIdx[o] <- m
        }
        ## ignore
    }
    ## multiple fragments to one mz (closest is default because of
    ## which(lDiff >= rDiff); "ignore" isn't possible here
    redundantIonMatch <- match.arg(redundantIonMatch)
    if (redundantIonMatch == "remove") {
        lIdx[is.finite(lDiff) & is.finite(rDiff)] <- NA_integer_
    }

    as.integer(lIdx)
}

#' Reorder protein sequence.
#'
#' @param x `character`, sequence
#' @param method `character`, reorder method
#' @return `character`
#' @noRd
.reorderSequence <- function(x, method=c("original", "random", "inverse")) {
    stopifnot(is.character(x))
    method <- match.arg(method)

    if (method != "original") {
        x <- strsplit(x, "", fixed=TRUE)[[1L]]

        if (method == "random") {
            x <- sample(x)
        } else if (method == "inverse") {
            x <- rev(x)
        }
        x <- paste0(x, collapse = "")
    }
    x
}

#' Apply unimod modification 1: Acetyl
#'
#' Acetylation
#'
#' TODO: replace by unimod package
#'
#' @param x `data.frame`, generated by [MSnbase::calculateFragments()]
#' @param s `character`, sequence
#' @return modified `data.frame`
#' @noRd
.unimod1 <- function(x, s) {
    i <- startsWith(s, x$seq)
    x$mz[i] <- x$mz[i] + 42.010565
    x
}

#' Apply unimod modification 4:  Carbamidomethyl
#'
#' Carboxyamidomethylation
#'
#' TODO: replace by unimod package
#'
#' @param x `data.frame`, generated by [MSnbase::calculateFragments()]
#' @return modified `data.frame`
#' @noRd
.unimod4 <- function(x) {
    iCU <- grep("C|U", x$seq)
    x$mz[iCU] <- x$mz[iCU] + 57.021464
    x
}

#' Apply unimod modification 765: Met-loss
#'
#' N-terminal methionine removed if followed by A, C, G, P, S, T, or V
#'
#' TODO: replace by unimod package
#'
#' @param x `character`/`AAString`, sequence
#' @return `character`/`AAString` without ^M
#' @noRd
.unimod765 <- function(x) {
    gsub("^M([ACGPSTV])", "\\1", x)
}

#' Apply custom modification
#'
#' TODO: replace by unimod package
#'
#' @param fragments `data.frame`, generated by [MSnbase::calculateFragments()].
#' @param n `integer(1)`, length of peptide sequence.
#' @param id `character(1)`.
#' @param deltamz `numeric(1)`.
#' @param location `integer(1)`, 0 for "N-term", n+1 for "C-term" or index.
#' @param variable `logical(1)`, if TRUE the unmodified and modified fragments
#' are returned.
#' @return modified `data.frame`
#' @noRd
.cusmod <- function(x, n, id, deltamz, location, variable) {
    nterm <- grepl("^a|^b|^c", x$type) & location <= x$pos
    cterm <- grepl("^x|^y|^z", x$type) & location > n - x$pos

    i <- nterm | cterm

    m <- x
    m$mz[i] <- m$mz[i] + deltamz
    m$ion[i] <- paste0(m$ion[i], "_m", id)
    m$type[i] <- paste0(m$type[i], "_m", id)

    if (variable)
        rbind(x, m[i, ], make.row.names = FALSE)
    else
        m
}
