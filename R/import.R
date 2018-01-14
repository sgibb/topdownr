#' List TopDown files
#'
#' List all TopDown files:
#'  - .fasta (peptide sequence)
#'  - .mzML (spectra)
#'  - .experiments.csv (fragmentation conditions)
#'  - .txt (header information)
#'
#' @param path file path
#' @param pattern filename pattern
#' @return list (splitted by file extension) with file path
#' @noRd
.listTopDownFiles <- function(path, pattern=".*") {
    files <- list.files(
        path,
        pattern=paste0(pattern, "(", .topDownFileExtRx("cfmt"), ")"),
        recursive=TRUE,
        full.names=TRUE
    )
    l <- split(normalizePath(files), .fileExt(files))
    n <- lengths(l)

    ext <- c("csv", "fasta", "mzML", "txt")

    if (!length(n) || any(!ext %in% names(l))) {
        stop(
            "Could not find any ",
            paste0(ext[!ext %in% names(l)], collapse=", "), " files!"
        )
    }

    if (n["fasta"] > 1L) {
        stop("More than one fasta file found. Consider the 'pattern' argument.")
    }

    if (!all(n["csv"] == n[!grepl("fasta", names(n))])) {
        nd <- n[!grepl("fasta", names(n))]
        stop(
            "There have to be the same number of csv, mzML and txt files. ",
            "Found: ", paste0(names(nd), "=", nd, collapse=", ")
        )
    }
    l
}

#' Read ScanHeadMans method (experiments.csv) output.
#'
#' 1 row per condition
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return data.frame
#' @noRd
.readExperimentCsv <- function(file, verbose=interactive()) {
    stopifnot(.fileExt(file) == "csv")
    d <- read.csv(file, na.strings=c("NA", "N/A"), stringsAsFactors=FALSE)
    colnames(d) <- .camelCase(colnames(d))

    .msg(verbose,
        "Reading ", nrow(d), " experiment conditions from file ", basename(file)
    )

    ## drop MS1
    d <- d[d$MsLevel == 2L, ]

    d$Condition <- seq_len(nrow(d))
    d$Mz <- .targetedMassListToMz(d$TargetedMassList)
    d$File <- gsub(.topDownFileExtRx("csv"), "", basename(file))
    d
}

#' Read fasta file.
#'
#' Stupid fasta reader, ignores comments, and keeps just the first sequence.
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return character
#' @noRd
.readFasta <- function(file, verbose=interactive()) {
    aa <- readAAStringSet(file, nrec=1L, use.names=FALSE)[[1L]]
    .msg(verbose, "Reading sequence from fasta file ", basename(file))
    if (!length(aa)) {
        stop("No sequence found.")
    }
    aa
}

#' Read ScanHeadMans header (txt) output.
#'
#' 1 row per scan (could have more rows than experiment.csv, because multiple
#' scans per condition are possible).
#'
#' @param file character, filename
#' @param verbose logical, verbose output?
#' @return data.frame
#' @noRd
.readScanHeadsTable <- function(file, verbose=interactive()) {
    stopifnot(.fileExt(file) == "txt")
    d <- read.csv(file, na.strings=c("NA", "N/A"), stringsAsFactors=FALSE)
    colnames(d) <- .camelCase(colnames(d))

    .msg(verbose,
        "Reading ", nrow(d), " header information from file ", basename(file)
    )

    ## drop MS1
    d <- d[d$MsOrder == 2L, ]

    ## Somehow the FilterString doesn't always contains the right mass
    ## label and we have to correct them.
    ## See also:
    ## - https://github.com/sgibb/topdownr/issues/25
    fixedFilterStrings <- .fixFilterString(d$FilterString)
    if (nFixed <- sum(fixedFilterStrings != d$FilterString)) {
        warning(
            nFixed, " FilterString entries modified because of ",
            "duplicated ID for different conditions.", immediate.=verbose
        )
    }
    d$FilterString <- fixedFilterStrings

    d$Condition <- as.integer(.filterStringToId(d$FilterString))

    ## Sometimes skips happen and the ID is not the same as in the FilterString
    ## (happens often in the missing_scans files)
    ## See also:
    ## - https://github.com/sgibb/topdownr/issues/14
    if (is.unsorted(d$Condition)) {
        warning(
            "ID in FilterString are not sorted in ascending order. ",
            "Introduce own condition ID via 'cumsum'.", immediate.=verbose
        )
        d$Condition <- cumsum(
            c(TRUE, d$FilterString[-1L] != d$FilterString[-nrow(d)])
        )
    }

    activation <- c("ETD", "CID", "HCD", "UVPD")
    activationColumns <- .camelCase(paste0(activation, "Activation"))
    d[, activationColumns] <- NA_real_

    for (i in seq(along=activation)) {
        ## first activation
        sel <- d$Activation1 == activation[i]
        if (sum(sel)) {
            d[sel, activationColumns[i]] <- d$Energy1[sel]
        }
        ## second activation
        if (activation[i] %in% c("CID", "HCD") &&
            "Activation2" %in% names(d)) {
            sel <- d$Activation2 == activation[i]
            sel[is.na(sel)] <- FALSE
            if (sum(sel)) {
                d[sel, activationColumns[i]] <- d$Energy2[sel]
            }
        }
    }

    d$Activation <- .fragmentationMethod(d[, activationColumns])

    d$File <- gsub(.topDownFileExtRx("txt"), "", basename(file))
    d
}

#' Read MS2 Spectra (mzML)
#'
#' @param file character, filename
#' @param fmass double, fragment mass
#' @param \ldots further arguments passed to .matchFragments
#' @param verbose logical, verbose output?
#' @return list (with data.frame for header and sparseMatrix with intensity
#' values)
#' @noRd
.readMzMl <- function(file, scans, fmass, ..., verbose=interactive()) {
    .msg(verbose,
        "Reading spectra information from file ", basename(file),
        appendLF=FALSE
    )

    fh <- openMSfile(file)
    on.exit(close(fh))

    hd <- header(fh)
    i <- which(hd$msLevel == 2L & hd$acquisitionNum %in% scans)
    hd <- hd[i, !grepl("injectionTime", colnames(hd), fixed=TRUE), drop=FALSE]
    colnames(hd)[grepl("acquisitionNum", colnames(hd), fixed=TRUE)] <- "Scan"
    colnames(hd)[grepl("seqNum", colnames(hd), fixed=TRUE)] <- "SpectrumIndex"
    hd$File <- gsub(.topDownFileExtRx("mzml"), "", basename(file))
    colnames(hd) <- .camelCase(colnames(hd))

    nr <- nrow(hd)
    m <- Matrix(0L, nrow=length(fmass), ncol=nr, sparse=TRUE)

    for (j in seq_along(i)) {
        k <- .matchFragments(peaks(fh, i[j])[, 1L], fmass, ...)
        notNA <- !is.na(k)
        if (sum(notNA)) {
            m[k[notNA], j] <- peaks(fh, i[j])[notNA, 2L]
        }
    }

    ## depending on the used processing software the header column totIonCurrent
    ## doesn't take deisotoping and charge state reduction into account; so we
    ## calculate TIC for our own here
    hd$TotIonCurrent <- .vapply1d(peaks(fh)[i], function(ii)sum(ii))

    .msg(verbose,
        sprintf(
            " (%02.1f%%)",
            round(sum(m != 0L) / sum(hd$PeaksCount) * 100, 1L)
        )
    )

    list(hd=hd, m=m)
}

#' Read single (complete) MS2 Spectrum (mzML)
#'
#' @param file `character`, filename
#' @param i `integer`, scan idx (in file)
#' @return `matrix` (`mzR::peaks`)
#' @noRd
.readSpectrum <- function(file, i) {
    fh <- openMSfile(file)
    on.exit(close(fh))

    if (i <= 0L || i > runInfo(fh)$scanCount) {
        stop("Invalid spectrum index. It has to be 1:",
             runInfo(fh)$scanCount, ".")
    }
    peaks(fh, i)
}

#' Merge ScanCondition and HeaderInformation
#'
#' @param sc data.frame, scan conditions
#' @param hi data.frame, header information
#' @return data.frame
#' @noRd
.mergeScanConditionAndHeaderInformation <- function(sc, hi) {
    stopifnot(is(sc, "data.frame"))
    stopifnot(is(hi, "data.frame"))
    d <- merge(
        sc, hi,
        by=c("File", "Condition"),
        suffixes=c(".ScanCondition", ".HeaderInformation")
    )
    naSACe <- is.na(d$SupplementalActivationCe)
    naCid <- is.na(d$CidActivation)
    naHcd <- is.na(d$HcdActivation)
    if (any(d$SupplementalActivationCe[!naSACe] != d$Energy2[!naSACe]) ||
        any(d$SupplementalActivationCe[!(naSACe | naCid)] !=
            d$CidActivation[!(naSACe | naCid)]) ||
        any(d$SupplementalActivationCe[!(naSACe | naHcd)] !=
            d$HcdActivation[!(naSACe | naHcd)])) {
        stop("Merging of header and method information failed. ",
             "Differences in 'SupplementalActivationCe', 'Energy2', ",
             "'CidActivation' and/or 'HcdActivation' found.")
    }
    d
}

#' Merge spectra and ScanConditions/HeaderInformation (into featureData slot)
#'
#' @param mzml data.frame, header from mzML files
#' @param scdm data.frame, header from ScanHeadsman
#' @return merged data.frame
#' @noRd
.mergeSpectraAndHeaderInformation <- function(mzml, scdm) {
    merge(
        mzml, scdm, sort=FALSE,
        by=c("File", "Scan"),
        suffixes=c(".SpectraInformation", ".HeaderInformation")
    )
}
