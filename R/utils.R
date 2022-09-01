#' Test whether all values in a vector are identical
#' (all(x[1] == x) is slower than unique)
#'
#' @param x vector
#' @return `logical`, TRUE/FALSE
#' @noRd
.allIdentical <- function(x) {
    isTRUE(length(unique(x)) == 1L)
}

#' Create camelCase names.
#'
#' @param x `character`
#' @return `character`
#' @noRd
.camelCase <- function(x) {
    ## convert dots from `make.names` to space
    x <- gsub("\\.+", " ", x)
    ## remove all the other punctiation like (, ), /
    x <- gsub("[[:punct:]]+", "", x)
    ## split AGCTarget to AGC Target and/or
    ## SupplementalActivationCE to Supplemental Activation CE
    x <- gsub("([A-Z]+)(?=([a-z]|[0-9]|[A-Z]$))", " \\1", x, perl=TRUE)
    ## Don't allow UPVDTime but
    ## don't split TIC/ETD to single letters (that's why [^ ])
    x <- gsub("(?<=[^ ][A-Z])([A-Z])(?=[a-z])", " \\1", x, perl=TRUE)
    ## just capitalize the first letter in a word
    x <- gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(x), perl=TRUE)
    ## remove whitespace
    x <- gsub(" ", "", x)
    x
}

#' cat0, cat with sep="", similar to paste0
#'
#' @param \ldots arguments passed to cat
#' @noRd
cat0 <- function(...) {
    cat(..., sep="", append=TRUE)
}

#' Convert character with On/Off, True/False to logical
#'
#' @param x `character`
#' @param na.strings `character`
#' @return `logical`/`character`
#' @noRd
.characterToLogical <- function(x, na.strings=c("NA", "N/A")) {
    if ((!is.character(x) && !is(x, "Rle")) ||
        (is(x, "Rle") && !is.character(runValue(x))))
        stop("'x' has to be a 'character'.")
    y <- gsub("^[[:space:]]*on[[:space:]]*$", "TRUE", decode(x),
              ignore.case=TRUE)
    y <- gsub("^[[:space:]]*off[[:space:]]*$", "FALSE", y, ignore.case=TRUE)

    l <- as.logical(y)

    if (all(
        grepl(paste0(na.strings, collapse="|"), y) | is.na(y) == is.na(l))
    ) l else x
}

#' file extension
#'
#' @param x `character`, file path
#' @param compression `logical`, if TRUE gz, bz2, xz is removed
#' @return `character`
#' @noRd
.fileExt <- function(x, compression=TRUE) {
    if (compression)
        x <- sub("[.](gz|bz2|xz|zip)$", "", x)
    file_ext(x)
}

#' The ScanHeadsman output for the header information contains a column
#' FilterString with the format "FTMS + p NSI Full ms2 [0-9]+\.[0-9]+@hcd35.00
#' [xxx-yyy]". This function converts this format to the ID stored in the mass
#' label.
#'
#' @param x `character`
#' @return `double`
#' @noRd
.filterStringToId <- function(x) {
    if (!is.character(x))
        stop("'x' has to be a 'character'.")
    .massLabelToId(gsub("^.*ms2 ([^@]+)\\@.*$", "\\1", x))
}

#' Ensure unique, increasing ids in FilterStrings
#'
#' The ScanHeadsman output contains some FilterStrings that have wrong IDs (the
#' id from the previous or next run), e.g.:
#'
#' FTMS + p NSI Full ms2 1162.0007@cid28.00 [100.0000-2000.0000]
#' FTMS + p NSI Full ms2 1162.0009@hcd28.00 [100.0000-2000.0000]
#' FTMS + p NSI Full ms2 1162.0009@hcd28.00 [100.0000-2000.0000]
#' FTMS + p NSI Full ms2 1162.0009@cid35.00 [100.0000-2000.0000]
#'
#' pavel-shliaha finds out, that the skip is never larger than 1, never more
#' than 1 in a row and both is possible 1, 2, 2, 4 and 1, 3, 3, 4.
#'
#' In general it seems that the fix is as easy as ensuring that the FilterString
#' ids following a sequence (with some drop-outs).
#'
#' See the following issues for details:
#' - https://github.com/sgibb/topdownr/issues/14
#' - https://github.com/sgibb/topdownr/issues/25
#'
#' TODO: replace this with a better approach if we understand why some ID are
#' skipped.
#'
#' @param x `character`, FilterString
#' @return `character`, fixed FilterString with equal length as `x`
#' @noRd
.fixFilterString <- function(x) {
    ufs <- unique(x)
    ids <- .filterStringToId(ufs)
    pos <- regexpr("[0-9]{3}@", ufs)
    i <- match(x, ufs)
    substr(x, pos[i], pos[i] + 3L) <-
        sprintf("%03d@", .fixFilterStringId(ids))[i]
    x
}

#' Workhorse function of .fixFilterString
#'
#' There are two patterns possible: 1, 2, 2, 4 and 1, 3, 3, 4
#'
#' @param x `numeric`, id
#' @return `numeric`, fixed id
#' @noRd
.fixFilterStringId <- function(x) {
    d <- diff(x)
    ## 1, 2, 2, 4: d == d(1, 0, 2)
    l <- c(FALSE, d == 0L & c(d[-1L], Inf) == 2L)
    x[l] <- x[l] + 1L
    ## 1, 3, 3, 4: d == d(2, 0, 1)
    l <- c(FALSE, d == 2L & c(d[-1L], Inf) == 0L)
    x[l] <- x[l] - 1L
    ## edge case: 2, 2, 3, 4
    if (all(d[1L:2L] == c(0L, 1L))) {
        x[1L] <- x[1L] - 1L
    }
    ## edge case: 1, 2, 3, 3
    if (all(d[length(d) - 0L:1L] == c(0L, 1L))) {
        x[length(x)] <- x[length(x)] + 1L
    }
    x
}

#' Flat nested lists to unnested one
#'
#' @param x `list`
#' @return `list`
#' @noRd
.flatten <- function(x) {
    .is.list <- function(x)inherits(x, "list")
    while(any(.vapply1l(x, .is.list))) {
        x <- unlist(
            lapply(x, function(y) { if (.is.list(y)) { y } else { list(y) } }),
            recursive=FALSE
        )
    }
    x
}

#' Convert number to string and prepend zeros
#'
#' @param x `double`
#' @param asInteger `logical`
#' @param nScientific `numeric`
#' @param na2zero `logical`
#' @return `character`
#' @noRd
.formatNumbers <- function(x, asInteger=NA_integer_, nScientific=5L, na2zero=FALSE) {
    ## sprintf("%02.2f", Rle(c(1.2, 2.2), 2:3)) == "01" instead c("1.20", "1.20",
    ## "2.20", "2.20", "2.20"); so we need to decode
    x <- decode(x)
    n <- .ndigits(max(abs(x), 0, na.rm=TRUE))
    na <- is.na(x)

    if (isTRUE(any(x < 0L))) {
        n <- n + 1
    }

    if (na2zero) {
        x[na] <- 0
    }

    asInteger <- isTRUE(asInteger) ||
        (is.na(asInteger) &&
         all(abs(x[!na] - as.integer(x[!na])) < .Machine$double.eps))

    if (asInteger && n <= nScientific) {
        sprintf(paste0("%0", n, "d"), as.integer(x))
    } else if (asInteger && n > nScientific) {
        sprintf(paste0("%.1e"), as.integer(x))
    } else {
        # + 3L = 2 place after the decimal point + the point itself
        sprintf(paste0("%0", n + 3L, ".2f"), x)
    }
}

#' Get fragmentation method from {ETD,CID,HCD}Activation
#'
#' @param x `data.frame`/`matrix`, 3 columns (ETD,CID,HCD)
#' @return `character`, fragmentation method
#' @noRd
.fragmentationMethod <- function(x) {
    methods <- c(
        "None", "ETD", "CID", "ETcid", "HCD", "EThcd", "CID/HCD",
        "ETD/CID/HCD", "UVPD", rep(NA_character_, 6), "All"
    )
    v <- c(EtdActivation=1L, CidActivation=2L, HcdActivation=4L,
           UvpdActivation=8L)
    if (!all(colnames(x) %in% names(v)))
        stop("Unknown activation in 'x'")
    x <- x[, names(v)]
    apply(x, 1L, function(i)methods[sum(v[as.logical(i)], na.rm=TRUE) + 1L])
}

#' Split list/data.frame
#'
#' @param x `data.frame`, e.g. from .ms2Experiments
#' @param cols `character`, colnames used to split
#' @return `list`
#' @noRd
.groupBy <- function(x, cols) {
    split(x, .groupByLabels(x, cols))
}

#' Create group labels from data.frame columns
#'
#' @param x `data.frame`, or `DataFrame` e.g. from .ms2Experiments
#' @param cols `character`, colnames used to split
#' @param sep `character`, separator
#' @return `list`
#' @noRd
.groupByLabels <- function(x, cols=names(x), sep=":") {
    x <- as.data.frame(x)
    if (any(!cols %in% colnames(x))) {
        stop("All 'cols' have to be valid column names of 'x'.")
    }
    if (length(cols) > 1L) {
        ## `interaction` doesn't handle NA values, so use `paste` instead
        do.call(paste, c(x[, cols], sep=sep))
    } else {
        as.character(x[, cols])
    }
}

#' Create group id
#'
#' @param x `data.frame`
#' @param cols `character`, colnames used for grouping
#' @return `integer`
#' @noRd
.groupId <- function(x, cols) {
    if (!is.data.frame(x) && !inherits(x, "DataFrame"))
        stop("'x' has to be a 'data.frame' or 'DataFrame'")
    groups <- .groupByLabels(x, cols)
    ugroups <- unique(groups)
    match(groups, ugroups)
}

#' head/fill/tail of vector elements (e.g. for printing)
#'
#' @param x vector
#' @param fill `character`, or NULL
#' @param n `integer`
#' @return vector, `c(head(x, n), fill, tail(x, n))`
#' @noRd
.hft <- function(x, fill="...", n=3) {
    if (length(x) <= 2 * n) {
        x
    } else {
        c(head(x, n), fill, tail(x, n))
    }
}

#' Test for numeric similarity/equality (similar to all.equal)
#'
#' @param x `numeric`
#' @param y `numeric`
#' @param tolerance `double`
#' @return TRUE/FALSE
#' @noRd
.isEqual <- function(x, y, tolerance=1e-3) {
    if (!is.numeric(x) || !is.numeric(y))
        stop("Both 'x' and 'y' have to be 'numeric'")
    if (length(x) != length(y) && length(x) != 1L && length(y) != 1L)
        stop("Both 'x' and 'y' have to have the same length or length of one.")
    isTRUE(all(abs(y - x) < tolerance))
}

#' Add log message.
#'
#' @param ... arguments passed to `paste0`
#' @return msg with date string added
#' @noRd
.logmsg <- function(...) {
    paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

#' Make names (similar to base::make.names but start appending with 1)
#'
#' @param x `character` (or must be convertable to `character`)
#' @param prefix `character`, prefix
#' @param sep `character`, separator
#' @return `character`
#' @noRd
.makeNames <- function(x, prefix, sep=":") {
    x <- as.character(x)
    if (!missing(prefix)) {
        x <- paste0(prefix, x)
    }
    ave(x, x, FUN=function(xx) {
        if (length(xx) == 1L) {
            xx
        } else {
            sprintf(
                paste0("%s", sep, "%0", .ndigits(length(xx)), "d"),
                xx, seq_along(xx)
            )
        }
    })
}

#' Create mass label
#'
#' Identifying the experiments by the running time/order is complicated.
#' Sometimes the instrument records a new run with the same settings which moves
#' all indicies. Same is true for the start times.
#'
#' We modify the target mass slightly by rounding to one decimal place and use
#' the second to fourth (default) to encode the id. Subsequently it is possible
#' to find the results by their individual mass.
#'
#' @param x `double`, original mass
#' @param id `double`, run id
#' @param divisor `double`, divisor (determines which decimal place)
#' @return `double`, mass label (id encoded in the second to fourth decimal
#' place)
#' @seealso [.massLabelToId()]
#' @noRd
.massLabel <- function(x, id, divisor=10000L) {
    if (any(log10(divisor) <= log10(id) + 1L)) {
        stop("'divisor' has to be at least two digits more than 'id'")
    }
    round(x, 1L) + id / divisor
}

#' Extract ID from mass labels
#'
#' @param x `character`, mass label
#' @param idDigits `integer`, number of digits behind the decimal place that
#' store id information (not mass information) from the last one (e.g. 3 if
#' the id is 123 and the mz is 900.0123)
#' @seealso [.massLabel()]
#' @noRd
.massLabelToId <- function(x, idDigits=3L) {
    # was the following before, but results in round errors ("7" becomes 6L)
    # x <- as.double(x)
    # as.integer((x - round(x, 1L)) * multiplicator)
    n <- nchar(x)
    as.integer(substring(x, n - idDigits + 1L, n))
}

#' verbose output
#'
#' @param \ldots arguments passed to `message`
#' @noRd
.msg <- function(verbose, ...) {
    if (verbose) {
        message(...)
    }
}

#' number of digits
#'
#' @param x `double`
#' @return `integer`
#' @noRd
.ndigits <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x) || all(x == 0L)) {
        1L
    } else {
        trunc(log10(abs(x)) + 1L)
    }
}

#' similar to lengths but for rows
#'
#' @param x `list` of `data.frames`/`matrices`
#' @noRd
.nrows <- function(x) {
    if (!inherits(x, "list"))
        stop("'x' has to b a 'list'")
    .vapply1d(x, nrow)
}

#' Resample rows (Experiments) in a data.frame
#'
#' @param x `data.frame`
#' @param fun `function`, to apply on the number of rows, e.g. use `seq` to do
#' nothing.
#' @return reordered `data.frame`
#' @noRd
.resample <- function(x, fun=sample) {
    fun <- match.fun(fun)
    x[fun(nrow(x)), , drop=FALSE]
}

#' generate valid scan description
#'
#' @param x `data.frame`
#' @param replications `numeric`, number of replications
#' @param prefix `character`, prefix
#' @return unique names
#' @noRd
.scanDescription <- function(n, replications=1L, prefix="C") {
    if (!is.numeric(n) || length(n) != 1L)
        stop("'n' has to be a 'numeric' of length one.")
    if (!is.numeric(replications) || length(replications) != 1L ||
        replications <= 0L)
        stop(
            "'replications' has to be a 'numeric' of length one and ",
            "greater than 0.")
    if (!is.character(prefix) ||
        (length(prefix) != 1L && length(prefix) != n) ||
        any(!nzchar(prefix)))
        stop(
            "'prefix' has to be a 'character' of length one or ",
            "equal the length of 'n'."
        )
    desc <- paste0(
        prefix,
        rep.int(.formatNumbers(seq_len(n)), replications)
    )
    if (replications > 1L)
        desc <- paste0(desc, "R", rep(seq_len(replications), each=n))
    if (any(nchar(desc) > 16L))
        stop("Calcium devices accept just 16 characters for 'ScanDescription'.")
    if (anyDuplicated(desc))
        stop("Duplicated 'ScanDescription' are not allowed.")
    desc
}

#' shortend string to width and place "..." in the middle
#'
#' @param x `character`
#' @param width number of letters allowed/width of terminal
#' @return `character`
#' @noRd
.snippet <- function(x, width=getOption("width")) {
    nc <- nchar(x)
    w <- (width - 2L:3L) %/% 2L

    ifelse(
        nc <= width,
        x,
        paste0(
            substring(x, 1L, w[1L]), "...", substring(x, nc - w[2L] + 1L, nc)
        )
    )
}

#' normalize subset (turn `logical`, `numeric` and `character` to `numeric`)
#'
#' @param i ANY, subset
#' @param n `integer`, length of the original object
#' @param nms `character`, names of the orignal object
#' @return `logical`, vector of length n
#' @noRd
.subset <- function(i, n, nms=NULL) {
    ## use decode to turn Rle into native vectors, does nothing if i is
    ## already a native vector
    i <- decode(i)
    if (anyNA(i))
        stop("Subsetting by 'NA' is not supported.")
    if (!is.null(nms) && length(nms) != n)
        stop("If 'nms' is given it has to have the same length as 'n'")
    switch(
        class(i),
        "character" = .subsetByCharacter(i, nms),
        "logical" = .subsetByLogical(i, n),
        "integer" = ,
        "double" = ,
        "numeric" = .subsetByNumeric(i, n),
        stop("Unknown index class.")
    )
}

#' subset by `character`
#'
#' @param i `character`, subset
#' @param nms `character`, names of the orignal object
#' @return `numeric`
#' @noRd
.subsetByCharacter <- function(i, nms=NULL) {
    if (is.null(nms))
        return(integer())
    if (!is.character(i))
        stop("'i' has to be a 'character'.")
    if (!is.character(nms))
        stop("'nms' has to be a 'character'.")

    ii <- match(i, nms)
    if (anyNA(ii))
        stop("Subscript out of bound: ", paste0("'", i[is.na(ii)], "'",
                                                collapse=", "))
    ii
}

#' subset by `logical`
#'
#' @param i `logical`, subset
#' @param n `integer`, length
#' @return `numeric`
#' @noRd
.subsetByLogical <- function(i, n) {
    if (!is.logical(i))
        stop("'i' has to be a 'logical'.")
    if (!is.numeric(n))
        stop("'n' has to be a 'numeric'.")

    which(rep_len(i, n))
}

#' subset by `numeric`
#'
#' @param i `numeric`, subset
#' @param n `integer`, length
#' @return `numeric`
#' @noRd
.subsetByNumeric <- function(i, n) {
    if (!is.numeric(i))
        stop("'i' has to be a 'numeric'.")
    if (!is.numeric(n))
        stop("'n' has to be a 'numeric'.")

    if (any(n < i))
        stop("Subscript out of bound: ",
            paste0("'", i[n < i], "'", collapse=", "))
    i
}

#' subset filenames
#'
#' @param files `character`, filenames
#' @param selFiles `character`, selectedFiles, basename without extension
#' @return `character`
#' @noRd
.subsetFiles <- function(files, selFiles) {
    gsub(.topDownFileExtRx("cfmt"), "", files) %in% selFiles
}

#' swap file extensions
#'
#' @param x `character`, file name
#' @param ext `character`, new extension
#' @return `character`
#' @noRd
.swapFileExt <- function(x, ext="meth") {
    paste(file_path_sans_ext(x), ext, sep=".")
}

#' The ScanHeadsMan output for the scan conditons contains a column
#' TargetedMassList with the format "(mz=[0-9]+\.[0-9]+ z=[0-9]+ name=)". This
#' function converts this format to truncated (one decimal place) mz values.
#'
#' @param x `character`
#' @return `double`
#' @noRd
.targetedMassListToMz <- function(x) {
    if (!is.character(x))
        stop("'x' has to be a 'character'.")
    trunc(as.double(gsub("^.*mz=([^ ]+) z.*$", "\\1", x)) * 10L) / 10L
}

#' TopDown file extensions.
#'
#' @param type `character`, which file ext
#' @return `character`, regexp for file extensions
#' @noRd
.topDownFileExtRx <- function(type=c("cfmt", "cmt", "csv", "fasta", "mzml",
                                     "txt", "raw", "all")) {
    type <- match.arg(type)
    ext <- c(
        csv="experiments\\.csv",
        fasta="fasta",
        mzml="mz[Mm][Ll]",
        raw="raw",
        txt="txt"
    )
    sel <- switch(
        type,
        "all" = seq_along(ext),
        "cfmt" = c("csv", "fasta", "mzml", "txt"),
        "cmt" = c("csv", "mzml", "txt"),
        type
    )
    paste0("\\.", ext[sel], "(\\.(gz|bz2|xz|zip))?$", collapse="|")
}

#' Find indicies of top N values
#'
#' @param x `double`/`character`, value of interest
#' @param groupBy `character`, grouping variable
#' @return `integer`, indices
#' @noRd
.topIdx <- function(x, groupByLabels, n) {
    if (!is.double(x) && !is.numeric(x) && !is.character(x)) {
        stop("'x' has to be of type 'double' or 'character'")
    }
    if (n < 1L) {
        stop("'n' has to be greater or equal than 1.")
    }
    if (length(x) != length(groupByLabels)) {
        stop("'length(x)' and 'length(groupByLabels)' have to be equal.")
    }
    o <- order(x, decreasing=TRUE, na.last=TRUE)
    i <- unlist(
        lapply(split(o, groupByLabels[o]), "[", seq_len(n)), use.names=FALSE
    )
    i[!is.na(i)]
}

#' Due to changes in ProteomDiscover some mzML files have a different
#' `spectrumId` format that is not supported by proteowizard (and hence by
#' mzR), e.g. "scan=[acquisitionNum], file=[fileId]" vs the supported
#' "controllerType=0 controllerNumber=1 scan=[acquisitionNum]".
#' As fallback solution proteowizard just use 1:n as acquisitionNum which cause
#' mismatches to the ScanHeadsman output. See
#' https://github.com/sgibb/topdownr/issues/73 for details.
#'
#' @param x `character`
#' @return `double`
#' @noRd
.translateThermoIdToScanId <- function(x) {
    if(!is.character(x) || !length(x) || !all(nzchar(x)))
        stop("'x' has to be a non-empty character vector.")
    as.double(gsub("^.*scan=([0-9]+).*$", "\\1", x))
}

#' Create a valid filename
#'
#' According to
#' https://docs.microsoft.com/de-de/windows/desktop/FileIO/naming-a-file#naming_conventions
#' ":" etc are not allowed in file names on windows
#'
#' @param x `character`
#' @return `character`
#' @noRd
.validFilename <- function(x) {
    if(!is.character(x) || !length(x) || !all(nzchar(x)))
        stop("'x' has to be a non-empty character vector.")
    gsub("[:;\\<>|\\?\\*]", "-", gsub("\\.+", ".", x))
}

#' wrapper around vapply for FUN.VALUE=double(1L)
#' @noRd
.vapply1d <- function(X, FUN, ..., USE.NAMES=FALSE) {
    vapply(X=X, FUN=FUN, FUN.VALUE=NA_real_, ..., USE.NAMES=USE.NAMES)
}

#' wrapper around vapply for FUN.VALUE=logical(1L)
#' @noRd
.vapply1l <- function(X, FUN, ..., USE.NAMES=FALSE) {
    vapply(X=X, FUN=FUN, FUN.VALUE=NA, ..., USE.NAMES=USE.NAMES)
}
