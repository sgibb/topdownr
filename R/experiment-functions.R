#' Create fragment optimisation experiment
#'
#' @param ms1 `data.frame`, MS1 settings.
#' @param ... further named arguments with `data.frame`s containing the MS2
#' settings.
#' @param groupyBy `character`, group experiments by columns in the MS2
#' `data.frame`s. The columns have to be present in all `data.frame`s. Each
#' group will be written to its own XML file.
#' @param nMs2perMs1 `integer`,
#' @param duration `double`, how long should the scan be?
#' @param replications `integer`, number of replications.
#' @param randomise `logical`, should the MS2 scan settings randomised?
#' @param massLabeling `logical`,
#' how many MS2 scans should be run after a MS1 scan?
#' @return `list`, able to be written via [xml2::as_xml_document()]
#' @export
createExperimentsFragmentOptimisation <- function(
    ms1, ..., groupBy=c("AgcTarget", "replication"), nMs2perMs1=10, duration=0.5,
    replications=2, randomise=TRUE, massLabeling=TRUE) {

    ms2 <- do.call(.rbind, .flatten(list(...)))
    nr <- nrow(ms2)
    ms2 <- ms2[rep(seq_len(nr), replications),, drop=FALSE]
    ms2$replication <- rep(seq_len(replications), each=nr)
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

    if (mnrs > 999L) {
        stop(
            "Please choose a different 'groupBy' value or reduce the number ",
            "of combinations. We can't label more than 999 conditions and ",
            "you ask for ", mnrs, "."
        )
    }

    times <- .startEndTime(nMs2=mnrs, nMs2perMs1=nMs2perMs1, duration=duration)
    l <- vector(mode="list", length(ms2))
    names(l) <- names(ms2)

    for (i in seq(along=l)) {
        n <- ceiling(nrs[i] * 1L/nMs2perMs1) + nrs[i]
        sbtimes <- times[seq_len(n),]
        sbtimes$Id <- seq_len(n)
        l[[i]] <- list(MethodModifications=vector(mode="list", length=n))
        names(l[[i]][[1L]]) <- rep("Modification", n)

        l[[i]][[1L]][[1L]] <- .ms1ConditionToTree(ms1, times=sbtimes[1L, 2L:3L])
        attr(l[[i]][[1L]][[1L]], "Order") <- 1L
        idMs1 <- which(sbtimes$Type[-1L] == "MS1") + 1L
        idMs2 <- which(sbtimes$Type == "MS2")

        for (j in idMs1) {
            l[[i]][[1L]][[j]] <-
                .ms1CopyAndAppendExperiment(j - 1L, sbtimes[j, 2L:3L])
            attr(l[[i]][[1L]][[j]], "Order") <- j
        }
        for (j in seq(along=idMs2)) {
            l[[i]][[1L]][[idMs2[j]]] <- .ms2ConditionToTree(
                ms2[[i]][j,],
                expId=idMs2[j] - 1L,
                condId=j,
                times=sbtimes[idMs2[j], 2L:3L],
                massLabeling=massLabeling
            )
            attr(l[[i]][[1L]][[idMs2[j]]], "Order") <- idMs2[j]
        }

        attr(l[[i]][[1L]], "Version") <- 2L
        attr(l[[i]][[1L]], "Model") <- "OrbitrapFusionLumos"
        attr(l[[i]][[1L]], "Family") <- "Calcium"
        attr(l[[i]][[1L]], "Type") <- "SL"
    }
    l
}

.ms1ConditionToTree <- function(x, expId=0L, times, ...) {
    l <- list(Experiment=vector(mode="list", length=3L))
    names(l$Experiment) <- c(
        "FullMSScan", "StartTimeMin", "EndTimeMin"
    )
    l$Experiment$FullMSScan <- lapply(x[, !is.na(x)], as.list)
    l$Experiment[c("StartTimeMin", "EndTimeMin")] <- lapply(times, as.list)
    attr(l$Experiment, "ExperimentIndex") <- expId
    l
}

.ms1CopyAndAppendExperiment <- function(expId, times, ...) {
    l <- vector(mode="list", length=2L)
    names(l) <- c("CopyAndAppendExperiment", "Experiment")
    l$CopyAndAppendExperiment <- list()
    attr(l$CopyAndAppendExperiment, "SourceExperimentIndex") <- 0L
    l$Experiment <- list(
        StartTimeMin=list(times[1L]),
        EndTimeMin=list(times[2L])
    )
    attr(l$Experiment, "ExperimentIndex") <- expId
    l
}

#' Convert single data.frame row (condition) to nested list
#'
#' @param x `data.frame` row
#' @param expId `integer`, experiment id
#' @param condId `integer`, condition id
#' @param startTime `integer`, start time
#' @param endTime `integer`, end time
#' @param \dots arguments passed to internal functions
#' @return nested `list`
#' @noRd
.ms2ConditionToTree <- function(x, expId, condId, times, ...) {
    l <- list(Experiment=vector(mode="list", length=4L))
    names(l$Experiment) <- c(
        "TMSnScan", "MassListFilter",
        "StartTimeMin", "EndTimeMin"
    )
    l$Experiment$MassListFilter <-
        .massListToTree(x$MassList, id=condId, ...)
    attr(l$Experiment$MassListFilter, "MassListType") <-
        "TargetedMassInclusion"
    x[, c("MassList", "replication")] <- NA
    l$Experiment$TMSnScan <- lapply(x[, !is.na(x)], as.list)
    l$Experiment[c("StartTimeMin", "EndTimeMin")] <- lapply(times, as.list)
    attr(l$Experiment, "ExperimentIndex") <- expId
    l
}

.massListToTree <- function(x, id, massLabeling=TRUE) {
    x <- unname(.expandMassList(x))
    if (massLabeling) {
        x[, 1L] <- .massLabel(x[, 1L], id)
    }
    l <- vector(mode="list", length=nrow(x))
    names(l) <- rep("MassListRecord", length(l))
    for (i in seq(along=l)) {
        l[[i]] <- list(MOverZ=list(x[i, 1L]), Z=list(x[i, 2L]))
    }
    list(MassList=l)
}

#' Collapse MassList
#'
#' @param x `matrix`, 2 columns: mz, z
#' @return `character`
#' @noRd
.collapseMassList <- function(x) {
    stopifnot(is.matrix(x))
    stopifnot(ncol(x) == 2L)
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

#' Expand MS1 Conditions
#'
#' TODO:
#'
#' @param ... further named arguments, used to create the combination of
#' conditions.
#' @return `data.frame`
#' @seealso [validMs1Settings()]
#' @export
#' @examples
#' expandMs1Conditions(FirstMass=100, LastMass=400)
expandMs1Conditions <- function(...) {
    settings <- .flatten(list(...))
    .validateMsSettings(type="MS1", settings)
    expand.grid(settings, stringsAsFactors=FALSE)
}

#' Expand MS2 Conditions
#'
#' TODO:
#'
#' @param MassList `matrix`, 2 columns (mass, z).
#' @param ActivationType `character`, *ActivationType* for MS2, either CID, HCD, ETD, or
#' UVPD.
#' @param ... further named arguments, used to create the combination of
#' conditions.
#' @return `data.frame`
#' @seealso [validMs2Settings()]
#' @export
#' @examples
#' expandMs2Conditions(
#'      MassList=cbind(mz=c(560.6, 700.5, 933.7), z=rep(1, 3)),
#'      ActivationType="CID",
#'      OrbitrapResolution="R120K",
#'      IsolationWindow=1,
#'      MaxITTimeInMS=200,
#'      Microscans=as.integer(40),
#'      AgcTarget=c(1e5, 5e5, 1e6),
#'      CIDCollisionEnergy=c(NA, seq(7, 35, 7))
#' )
expandMs2Conditions <- function(MassList,
                                ActivationType=c("CID", "HCD", "ETD", "UVPD"),
                                ...) {
    ActivationType <- match.arg(ActivationType)
    settings <- .flatten(list(...))
    .validateMsSettings(type=ActivationType, settings)
    expand.grid(
        c(
            MassList=.collapseMassList(MassList),
            ActivationType=ActivationType,
            settings
        ),
        stringsAsFactors=FALSE
    )
}

#' List valid MS settings.
#'
#' These functions list settings for MS1 or MS2 that are supported by
#' *Thermo's XmlMethodChanger*.
#'
#' @return `matrix` with three columns:
#'  - name: element name
#'  - class: expected R class of the value
#'  - type: MS/ActivationType, e.g. MS1/MS2/ETD/...
#' @rdname validMsSettings
#' @export
#' @examples
#' validMs1Settings()
validMs1Settings <- function() {
    .validMsSettings("MS1")
}

#' @rdname validMsSettings
#'
#' @param type `character`, type of activation.
#' @export
#' @examples
#' validMs2Settings()
#' validMs2Settings("MS2")
#' validMs2Settings("ETD")
#' validMs2Settings(c("MS2", "ETD"))
validMs2Settings <- function(type=c("All", "MS2", "ETD", "CID", "HCD",
                                    "UVPD")) {
    type <- match.arg(type, several.ok=TRUE)
    if ("All" %in% type) {
        type <- c("MS2", "ETD", "CID", "HCD", "UVPD")
    }
    .validMsSettings(type)
}

#' List valid MS settings
#'
#' @param type `character`, MS1/MS2/Activation
#' @param version `character`, currently just Calcium3.1 supported
#' @return `matrix`
#' @noRd
.validMsSettings <- function(type, version="Calcium3.1") {
    stopifnot(is.character(type))
    stopifnot(is.character(version))
    m <- get(paste0(".validMsSettings", version))
    m[m[, "type"] %in% type,, drop=FALSE]
}

#' Validate a single MS setting against internal .validMsSettings (derivated
#' from XSD)
#'
#' @param name `character`, element name
#' @param value any type, value of element
#' @param type `character`, type of setting
#' @return `TRUE` if valid, else message
#' @noRd
.validateMsSetting <- function(name, value, type) {
    if (type %in% c("ETD", "CID", "HCD", "UVPD")) {
        type <- c("MS2", type)
    }
    settings <- .validMsSettings(type)
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
#' @param type `character`, MS1/MS2/ActivationType
#' @param settings `list`, named arguments used for validation
#' @return `TRUE` if valid, else stops with an error
#' @noRd
.validateMsSettings <- function(type=c("MS1", "MS2", "ETD", "CID", "HCD", "UVPD"),
                                settings) {
    type <- match.arg(type)
    validation <- mapply(
        .validateMsSetting,
        name=names(settings),
        value=settings,
        type=type,
        SIMPLIFY=FALSE
    )
    isValid <- .vapply1l(validation, isTRUE)

    if (any(!isValid)) {
        stop(paste0(unlist(validation[!isValid]), collapse="\n  "))
    }
    TRUE
}
