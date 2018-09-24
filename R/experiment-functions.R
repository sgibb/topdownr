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

    ms2 <- do.call(.rbind, list(...))
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

    if (mnrs > 999) {
        stop(
            "Please choose a different 'groupBy' value or reduce the number ",
            "of combinations. We can't label more than 999 conditions and ",
            "you ask for ", mnrs, "."
        )
    }

    times <- .startEndTime(nMs2=mnrs, nMs2perMs1=nMs2perMs1, duration=duration)
    l <- vector(mode="list", length(ms2))
    names(l) <- rep("MethodModifications", length(ms2))

    for (i in seq(along=l)) {
        n <- ceiling(nrs[i] * 1L/nMs2perMs1) + nrs[i]
        sbtimes <- times[seq_len(n),]
        sbtimes$Id <- seq_len(n)
        l[[i]] <- vector(mode="list", length=n)
        names(l[[i]]) <- rep("Modification", n)

        l[[i]][[1L]] <- .ms1ConditionToTree(ms1, times=sbtimes[1L, 2L:3L])
        attr(l[[i]][[1L]], "Order") <- 1L
        idMs1 <- which(sbtimes$Type[-1L] == "MS1") + 1L
        idMs2 <- which(sbtimes$Type == "MS2")

        for (j in idMs1) {
            l[[i]][[j]] <-
                .ms1CopyAndAppendExperiment(j - 1L, sbtimes[j, 2L:3L])
            attr(l[[i]][[j]], "Order") <- j
        }
        condId <- 0L
        for (j in idMs2) {
            condId <- condId + 1L
            l[[i]][[j]] <- .ms2ConditionToTree(
                ms2[[i]][condId,],
                expId=j - 1L,
                condId=condId,
                times=sbtimes[j, 2L:3L]
            )
            attr(l[[i]][[j]], "Order") <- j
        }

        attr(l[[i]], "Version") <- 2L
        attr(l[[i]], "Model") <- "OrbitrapFusionLumos"
        attr(l[[i]], "Family") <- "Calcium"
        attr(l[[i]], "Type") <- "SL"

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
        "TMSnScan", "TargetedInclusionMassListFilter",
        "StartTimeMin", "EndTimeMin"
    )
    l$Experiment$TargetedInclusionMassListFilter <-
        .massListToTree(x$MassList, id=condId, ...)
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
#' @export
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

#' Validate a single MS setting against internal .validMsSettings (derivated
#' from XSD)
#'
#' @param name `character`, element name
#' @param value any type, value of element
#' @param type `character`, type of setting
#' @return `TRUE` if valid, else message
#' @noRd
.validateMsSetting <- function(name, value, type) {
    isValidEntry <- .validMsSettings[, "name"] == name &
        (.validMsSettings[, "type"] == type |
         .validMsSettings[, "type"] == "MS2" &
         type != "MS1")

    if (!any(isValidEntry)) {
        return(
            paste0(
                name, " is not a valid element of type '", type, "'",
                " and/or not defined in MethodModification.xsd."
            )
        )
    }

    entry <- .validMsSettings[isValidEntry, "class"]
    cl <- typeof(value)

    if (isTRUE(cl == "character")) {
        isValidValue <-
            .vapply1l(value, function(x)grepl(paste0("(^|:)", x, "(:|$)"), entry))

        if (any(!isValidValue)) {
            return(
                paste0(
                    name, " could not be '", paste0(value[!isValidValue]),
                    "'. Should be one of '", gsub(":", ", ", entry), "'."
                )
            )
        }
    } else if (cl != entry) {
        return(
            paste0(name, " should be of class '", entry, "' but is '", cl, "'.")
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


