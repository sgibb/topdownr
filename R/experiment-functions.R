#' Create fragment optimisation experiment
#'
#' This function is used to create a tree-like `list` of
#' all combinations of a user-given set of MS1 and TMS2 settings for an
#' fragment optimisation experiment. The list could be written to an
#' Orbitrap Fusion Lumos method xml file using [writeMethodXmls()].
#'
#' @param ms1 `data.frame`, MS1 settings.
#' @param ... further named arguments with `data.frame`s containing the TMS2
#' settings.
#' @param groupBy `character`, group experiments by columns in the TMS2
#' `data.frame`s. The columns have to be present in all `data.frame`s. Each
#' group will be written to its own XML file.
#' @param nMs2perMs1 `integer`, how many TMS2 scans should be run after a MS1
#' scan?
#' @param scanDuration `double`, if greater than zero (e.g. `scanDuration=0.5`)
#' the Start/EndTimeMin are overwritten with a duration of `scanDuration`. If
#' `scanDuration` is zero (default) Start/EndTimeMin are not overwritten.
#' @param replications `integer`, number of replications.
#' @param randomise `logical`, should the TMS2 scan settings randomised?
#' @return `list`, able to be written via [xml2::as_xml_document()]
#' @export
#' @seealso [writeMethodXmls()],
#' [`expandMs1Conditions()`][expandMs1Conditions],
#' [`expandTms2Conditions()`][expandMs1Conditions]
#' @examples
#' ## build experiments within R
#' ms1 <- expandMs1Conditions(
#'     FirstMass=400,
#'     LastMass=1200,
#'     Microscans=as.integer(10)
#' )
#'
#' targetMz <- cbind(mz=c(560.6, 700.5, 933.7), z=rep(1, 3))
#' common <- list(
#'     OrbitrapResolution="R120K",
#'     IsolationWindow=1,
#'     MaxITTimeInMS=200,
#'     Microscans=as.integer(40),
#'     AgcTarget=c(1e5, 5e5, 1e6)
#' )
#'
#' cid <- expandTms2Conditions(
#'     MassList=targetMz,
#'     common,
#'     ActivationType="CID",
#'     CIDCollisionEnergy=seq(7, 35, 7)
#' )
#' hcd <- expandTms2Conditions(
#'     MassList=targetMz,
#'     common,
#'     ActivationType="HCD",
#'     HCDCollisionEnergy=seq(7, 35, 7)
#' )
#' etd <- expandTms2Conditions(
#'     MassList=targetMz,
#'     common,
#'     ActivationType="ETD",
#'     ETDReactionTime= 1:2,
#' )
#' etcid <- expandTms2Conditions(
#'     MassList=targetMz,
#'     common,
#'     ActivationType="ETD",
#'     ETDReactionTime= 1:2,
#'     ETDSupplementalActivation="ETciD",
#'     ETDSupplementalActivationEnergy=1:2
#' )
#' uvpd <- expandTms2Conditions(
#'     MassList=targetMz,
#'     common,
#'     ActivationType="UVPD"
#' )
#'
#' exps <- createExperimentsFragmentOptimisation(
#'     ms1=ms1, cid, hcd, etd, etcid, uvpd,
#'     groupBy=c("AgcTarget", "replication"), nMs2perMs1=10, scanDuration=0.5,
#'     replications=2, randomise=TRUE
#' )
#'
#' ## use different settings for CID
#' cid560 <- expandTms2Conditions(
#'     MassList=cbind(560.6, 1),
#'     common,
#'     ActivationType="CID",
#'     CIDCollisionEnergy=seq(7, 21, 7)
#' )
#' cid700 <- expandTms2Conditions(
#'     MassList=cbind(700.5, 1),
#'     common,
#'     ActivationType="CID",
#'     CIDCollisionEnergy=seq(21, 35, 7)
#' )
#'
#' exps <- createExperimentsFragmentOptimisation(
#'     ms1=ms1, cid560, cid700,
#'     groupBy=c("AgcTarget", "replication"), nMs2perMs1=10, scanDuration=0.5,
#'     replications=2, randomise=TRUE
#' )
#'
#' ## use a CSV (or excel) file as input
#' myCsvContent <- "
#' ActivationType, ETDReactionTime, UVPDActivationTime
#' UVPD,,1000
#' ETD,1000,
#' "
#' myCsvSettings <- read.csv(text=myCsvContent, stringsAsFactors=FALSE)
#' myCsvSettings
#' #   ActivationType ETDReactionTime UVPDActivationTime
#' # 1           UVPD              NA               1000
#' # 2            ETD            1000                 NA
#'
#' exps <- createExperimentsFragmentOptimisation(
#'     ms1 = data.frame(FirstMass=500, LastMass=1000),
#'     ## TMS2
#'     myCsvSettings,
#'     ## other arguments
#'     groupBy="ActivationType"
#' )
createExperimentsFragmentOptimisation <-
    function(ms1, ...,
             groupBy=c("AgcTarget", "replication"),
             nMs2perMs1=10, scanDuration=0,
             replications=2, randomise=TRUE) {

    stopifnot(
        is.numeric(nMs2perMs1), length(nMs2perMs1) == 1,
        is.numeric(scanDuration), length(scanDuration) == 1,
        is.numeric(replications), length(replications) == 1,
        is.logical(randomise), length(randomise) == 1
    )

    ms2 <- do.call(.rbind, .flatten(list(...)))
    nr <- nrow(ms2)

    ms2 <- ms2[rep(seq_len(nr), replications),, drop=FALSE]
    ms2$replication <- rep(seq_len(replications), each=nr)
    ms2$ScanDescription <- .scanDescription(n=nr, replications=replications)

    if (length(ms2)) {
        ms2 <- .groupBy(ms2, groupBy)
    } else {
        ms2 <- list(ms2=ms2)
    }
    if (randomise) {
        ms2 <- lapply(ms2, .resample)
    }

    nrs <- .nrows(ms2)
    mnrs <- max(nrs)

    times <- .startEndTime(
        nMs2=mnrs,
        nMs2perMs1=nMs2perMs1,
        duration=scanDuration
    )

    l <- vector(mode="list", length(ms2))
    names(l) <- names(ms2)

    for (i in seq(along=l)) {
        nMs1 <- ceiling(nrs[i] * 1L/nMs2perMs1)
        nMs2 <- nrs[i]
        n <- nMs1 + nMs2
        nNodes <- n + nMs2 - 1L     # CopyAndAppendExperiment

        sbtimes <- times[seq_len(n),]
        sbtimes$ExpId <- seq_len(n) - 1L
        sbtimes$SrcId <- as.integer(sbtimes$Type == "MS2")

        l[[i]] <- list(MethodModifications=vector(mode="list", length=nNodes))
        names(l[[i]][[1L]]) <- rep("Modification", nNodes)

        l[[i]][[1L]][[1L]] <- structure(
            .ms1ConditionToTree(
                ms1,
                times=unlist(
                    sbtimes[1L, c("StartTimeMin", "EndTimeMin")],
                    use.names=FALSE
                )
            ),
            Order=1L
        )

        for (j in seq_len(nrow(sbtimes) - 2L) + 2L) {
            orderId <- j - 1L
            l[[i]][[1L]][[orderId]] <- structure(
                .copyAndAppendExperiment(
                    srcId=sbtimes$SrcId[j]
                ),
                Order=orderId
            )
        }

        sbtimes <- sbtimes[sbtimes$Type == "MS2",, drop=FALSE]

        for (j in seq_len(nrow(sbtimes))) {
            orderId <- n + j - 1L
            l[[i]][[1L]][[orderId]] <- structure(
                .tms2ConditionToTree(
                    ms2[[i]][j,],
                    id=sbtimes$ExpId[j],
                    times=unlist(
                        sbtimes[j, c("StartTimeMin", "EndTimeMin")],
                        use.names=FALSE
                    )
                ),
                Order=orderId
            )
        }

        attr(l[[i]][[1L]], "Version") <- 2L
        attr(l[[i]][[1L]], "Model") <- "OrbitrapFusionLumos"
        attr(l[[i]][[1L]], "Family") <- "Calcium"
        attr(l[[i]][[1L]], "Type") <- "SL"
    }
    l
}

#' Create a FullMSScan node in a nested list
#'
#' @param x `data.frame` row
#' @param id `integer`, experiment id
#' @param times `double(2)`, start/end time
#' @param \dots arguments passed to internal functions
#' @return nested `list`
#' @noRd
.ms1ConditionToTree <- function(x, id=0L, times=NA_real_, ...) {
    if (anyNA(times)) {
        l <- list(Experiment=list(FullMSScan=NULL))
    } else {
        l <- list(Experiment=vector(mode="list", length=3L))
        names(l$Experiment) <- c(
            "FullMSScan", "StartTimeMin", "EndTimeMin"
        )
        l$Experiment$StartTimeMin <- list(times[1L])
        l$Experiment$EndTimeMin <- list(times[2L])
    }
    l$Experiment$FullMSScan <- lapply(x[, !is.na(x)], as.list)
    attr(l$Experiment, "ExperimentIndex") <- id
    l
}

#' Create a CopyAndAppendExperiment node in a nested list
#'
#' @param srcId `integer`, source experiment id
#' @param \dots arguments passed to internal functions
#' @return nested `list`
#' @noRd
.copyAndAppendExperiment <- function(srcId=0L, ...) {
    list(
        CopyAndAppendExperiment=structure(list(), SourceExperimentIndex=srcId)
    )
}

#' Convert single data.frame row (condition) to nested list
#'
#' @param x `data.frame` row
#' @param id `integer`, experiment id
#' @param times `double(2)`, start/end time
#' @param \dots arguments passed to internal functions
#' @return nested `list`
#' @noRd
.tms2ConditionToTree <- function(x, id, times=NA_real_, ...) {
    if (anyNA(times)) {
        l <- list(Experiment=list(TMSnScan=NULL))
    } else {
        l <- list(Experiment=vector(mode="list", length=3L))
        names(l$Experiment) <- c(
            "TMSnScan", "StartTimeMin", "EndTimeMin"
        )
        l$Experiment$StartTimeMin <- list(times[1L])
        l$Experiment$EndTimeMin <- list(times[2L])
    }
    if (!is.null(x$MassList)) {
        l$Experiment$TMSnScan$MassList <- .massListToTree(x$MassList)
    }
    x[, c("MassList", "replication")] <- NA
    cn <- colnames(x)[!is.na(x)]
    if (length(cn)) {
        l$Experiment$TMSnScan[cn] <- lapply(x[, cn], as.list)
    }
    attr(l$Experiment, "ExperimentIndex") <- id
    l
}

#' Convert collapsed MassList into xml2 ready tree like list structure
#'
#' @param x `character`, collapsed MassList (output of `.collapseMassList`)
#' @return `list`
#' @noRd
.massListToTree <- function(x) {
    x <- unname(.expandMassList(x))
    l <- vector(mode="list", length=nrow(x))
    names(l) <- rep("MassListRecord", length(l))
    for (i in seq(along=l)) {
        l[[i]] <- list(MOverZ=list(x[i, 1L]), Z=list(x[i, 2L]))
    }
    l
}

#' Collapse MassList
#'
#' @param x `matrix`, 2 columns: mz, z
#' @return `character`
#' @noRd
.collapseMassList <- function(x) {
    stopifnot(is.matrix(x), ncol(x) == 2L)
    paste0(apply(x, 1L, paste0, collapse="/"), collapse=" ")
}

#' Expand MassList
#'
#' @param x `character` collapsed MassList
#' @return `matrix`
#' @noRd
.expandMassList <- function(x) {
    stopifnot(is.character(x))
    x <- strsplit(x, " |/")[[1L]]
    matrix(
        as.numeric(x),
        ncol=2, byrow=TRUE,
        dimnames=list(NULL, c("mz", "z"))
    )
}

#' Expand MS Conditions
#'
#' Create a `data.frame` of all possible combinations of the given arguments.
#' It ensures that just arguments are applied that yield a valid
#' MethodModification.xml file.
#'
#' @param ... further named arguments, used to create the combination of
#' conditions.
#' @param family `character`, currently just Calcium is supported
#' @param version `character`, currently 3.1, 3.2 (default), 3.3 are supported
#' @return `data.frame` with all possible combinations of conditions/settings.
#' @seealso [validMs1Settings()]
#' @rdname expandMsConditions
#' @seealso [validTms2Settings()], [expand.grid()]
#' @export
#' @examples
#' expandMs1Conditions(FirstMass=100, LastMass=400)
expandMs1Conditions <- function(..., family="Calcium", version="3.2") {
    settings <- .flatten(list(...))
    .validateMsSettings(type="MS1", settings, family=family, version=version)
    expand.grid(settings, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
}

#' @rdname expandMsConditions
#' @param ActivationType `character`, *ActivationType* for TMS2,
#' either CID, HCD, ETD, or UVPD.
#' @param MassList `matrix`, 2 columns (mass, z) for targeted mass list,
#' or `NULL` (default) to not overwrite targeted mass.
#' @export
#' @examples
#' expandTms2Conditions(
#'      ActivationType="CID",
#'      OrbitrapResolution="R120K",
#'      IsolationWindow=1,
#'      MaxITTimeInMS=200,
#'      Microscans=as.integer(40),
#'      AgcTarget=c(1e5, 5e5, 1e6),
#'      CIDCollisionEnergy=c(NA, seq(7, 35, 7)),
#'      MassList=cbind(mz=c(560.6, 700.5, 933.7), z=rep(1, 3))
#' )
expandTms2Conditions <- function(ActivationType=c("CID", "HCD", "ETD", "UVPD"),
                                 ...,
                                 MassList=NULL,
                                 family="Calcium", version="3.2") {
    ActivationType <- match.arg(ActivationType)
    settings <- .flatten(list(...))

    .validateMsSettings(
        type=ActivationType,
        settings,
        family=family,
        version=version
    )

    if (!is.null(MassList)) {
        MassList <- .collapseMassList(MassList)
    }

    expand.grid(
        c(
            MassList=MassList,
            ActivationType=ActivationType,
            settings
        ),
        KEEP.OUT.ATTRS=FALSE,
        stringsAsFactors=FALSE
    )
}

#' List valid MS settings.
#'
#' These functions list settings for MS1 or TMS2 that are supported by
#' *Thermo's XmlMethodChanger*.
#'
#' @param family `character`, currently just Calcium is supported
#' @param version `character`, currently 3.1, 3.2 (default), 3.3 are supported
#' @return `matrix` with three columns:
#'  - name: element name
#'  - class: expected R class of the value
#'  - type: MS/ActivationType, e.g. MS1/TMS2/ETD/...
#' @rdname validMsSettings
#' @export
#' @examples
#' validMs1Settings()
validMs1Settings <- function(family="Calcium", version="3.2") {
    .validMsSettings("MS1", family=family, version=version)
}

#' @rdname validMsSettings
#'
#' @param type `character`, type of activation.
#' @export
#' @examples
#' validTms2Settings()
#' validTms2Settings("TMS2")
#' validTms2Settings("ETD")
#' validTms2Settings(c("TMS2", "ETD"))
validTms2Settings <- function(type=c("All", "TMS2", "ETD", "CID", "HCD", "UVPD"),
                              family="Calcium", version="3.2") {
    type <- match.arg(type, several.ok=TRUE)
    if ("All" %in% type) {
        type <- c("TMS2", "ETD", "CID", "HCD", "UVPD")
    }
    .validMsSettings(type, family=family, version=version)
}

#' List valid MS settings
#'
#' @param type `character`, MS1/TMS2/Activation
#' @param family `character`, currently just Calcium is supported
#' @param version `character`, currently 3.1, 3.2 [default], 3.3 are supported
#' @return `matrix`
#' @noRd
.validMsSettings <- function(type, family="Calcium", version="3.2") {
    stopifnot(
        is.character(type),
        is.character(family) && family %in% names(.validMsSettingsXsd),
        is.character(version) &&
            version %in% names(.validMsSettingsXsd[[family]])
    )
    m <- .validMsSettingsXsd[[family]][[version]]
    m[m[, "type"] %in% type,, drop=FALSE]
}

#' Validate a single MS setting against internal .validMsSettings (derivated
#' from XSD)
#'
#' @param name `character`, element name
#' @param value any type, value of element
#' @param type `character`, type of setting
#' @param family `character`, currently just Calcium is supported
#' @param version `character`, currently 3.1, 3.2 [default], 3.3 are supported
#' @return `TRUE` if valid, else message
#' @noRd
.validateMsSetting <- function(name, value, type, family="Calcium",
                               version="3.2") {
    if (type %in% c("ETD", "CID", "HCD", "UVPD")) {
        type <- c("TMS2", type)
    }
    settings <- .validMsSettings(type, family=family, version=version)
    entry <- settings[settings[, "name"] == name,]

    if (!length(entry)) {
        return(
            paste0(
                name, " is not a valid element of type ",
                paste0("'", type, "'", collapse="/"),
                " and/or not defined in MethodModification.xsd.\n",
                "  Run `validMs", ifelse(all(type == "MS1"), "1", "2"),
                "Settings()` for a complete list of possible settings."
            )
        )
    }

    tcl <- entry["class"]
    cl <- typeof(value)

    if (isTRUE(cl == "character" && grepl(":", tcl))) {
        isValidValue <-
            .vapply1l(value, function(x)grepl(paste0("(^|:)", x, "(:|$)"), tcl))

        if (any(!isValidValue)) {
            return(
                paste0(
                    name, " could not be '", paste0(value[!isValidValue]),
                    "'. Should be one of '", gsub(":", ", ", tcl), "'."
                )
            )
        }
    } else if (cl != tcl) {
        return(
            paste0(name, " should be of class '", tcl, "' but is '", cl, "'.")
        )
    }
    TRUE
}

#' Validate MS settings
#'
#' @param type `character`, MS1/TMS2/ActivationType
#' @param settings `list`, named arguments used for validation
#' @param family `character`, currently just Calcium is supported
#' @param version `character`, currently 3.1, 3.2 [default], 3.3 are supported
#' @return `TRUE` if valid, else stops with an error
#' @noRd
.validateMsSettings <- function(type=c("MS1", "TMS2", "ETD", "CID", "HCD",
                                       "UVPD"),
                                settings, family="Calcium", version="3.2") {
    type <- match.arg(type)
    validation <- mapply(
        .validateMsSetting,
        name=names(settings),
        value=settings,
        type=type,
        family=family,
        version=version,
        SIMPLIFY=FALSE
    )
    isValid <- .vapply1l(validation, isTRUE)

    if (any(!isValid)) {
        stop(paste0(unlist(validation[!isValid]), collapse="\n  "))
    }
    TRUE
}
