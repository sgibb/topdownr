#' Find best combination of conditions for highest coverage
#'
#' @param object `NCBSet`
#' @param n `integer`, max number of combinations/iterations.
#' @param minN `integer` stop if there are less than `minN` additional fragments
#' added.
#' @return `matrix`, first column: index of condition, second column: number of
#' newly covered bonds.
#' @export
#' @noRd
setMethod("bestConditions", "NCBSet",
          function(object, n=ncol(object), minN=0L) {
    m <- .bestCoverageCombination(object@assay, n=n, minN=minN)
    colnames(m) <- c("index", "bonds")
    m
})

#' Plot fragmentation map.
#'
#' @param object `NCBSet`
#' @param nCombinations `integer`, number of combinations to show (0 to avoid
#' plotting them at all)
#' @param cumCoverage `logical`, if `TRUE` (default) cumulative coverage of
#' combinations is shown.
#' @export
#' @noRd
setMethod("fragmentationMap", "NCBSet",
          function(object, nCombinations=10, cumCoverage=TRUE, ...) {
    d <- .dgcMatrix2data.frame(object@assay)
    d <- cbind(d,
               Condition=colnames(object)[d$col],
               Activation=object$Activation[d$col])

    if (nCombinations) {
        i <- bestConditions(object, nCombinations)[, "index"]
        combinations <- object[, i]

        if (cumCoverage) {
            cmb <- .dgcMatrix2data.frame(.cumComb(combinations@assay))
        } else {
            cmb <- .dgcMatrix2data.frame(combinations@assay)
        }
        cmb <- cbind(cmb,
                     Condition=paste(.formatNumbers(cmb$col),
                                     colnames(combinations)[cmb$col],
                                     sep=": "),
                     Activation="Cmb")
        cmb$Condition <- as.ordered(cmb$Condition)
        d <- rbind(d, cmb)
    }

    ## tell ggplot that we have an already ordered factor
    d$Condition <- as.ordered(d$Condition)
    d$Activation <- factor(d$Activation,
                           levels=c("CID", "HCD", "ETD", "ETcid", "EThcd",
                                    "Cmb"))

    ggplot(data=d) +
        geom_raster(aes(x=Condition, y=row, fill=as.factor(x))) +
        facet_grid(. ~ Activation, scales="free_x", space="free_x") +
        scale_fill_manual(name="Observed Fragments",
                          labels=c("N-terminal", "C-terminal", "Both"),
                          values=c("#1b9e77", "#d95f02", "#7570b3")) +
        scale_x_discrete(name="Condition", expand=c(0L, 0L)) +
        scale_y_discrete(name=expression(
                            paste("Bonds (", N-terminal %->% C-termnial, ")")),
                         expand=c(0L, 0L)) +
        geom_vline(xintercept=c(0L, seq_len(ncol(object))) + 0.5,
                   colour="#808080", size=0.1) +
        geom_hline(yintercept=c(0L, seq_len(nrow(object))) + 0.5,
                   colour="#808080", size=0.1) +
        ggtitle("fragmentation map") +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.4),
              plot.title=element_text(hjust=0.5, face="bold"),
              panel.grid.major=element_blank(),
              panel.border=element_blank(),
              panel.background=element_blank(),
              strip.background=element_rect(fill="#f0f0f0", colour="#ffffff"))
})

#' @param object `NCBSet`
#' @noRd
setMethod("show", "NCBSet", function(object) {
    callNextMethod()

    if (length(object@rowViews)) {
        cat("- - - Fragment data - - -\n")
        cat0(paste0("Number of ",
                    c("N-terminal", "C-terminal", "N- and C-terminal"),
                    " fragments: ", tabulate(object@assay@x, nbins=3L), "\n"))
    }

    if (nrow(object@colData)) {
        cat("- - - Condition data - - -\n")
        cat("Number of conditions:", length(unique(object$Sample)), "\n")
        cat("Number of scans:", nrow(object@colData), "\n")
        cat0("Condition variables (", ncol(object@colData), "): ",
             paste0(.hft(colnames(object@colData), n=2), collapse=", "), "\n")
    }

    if (all(dim(object))) {
        cat("- - - Assay data - - -\n")
        cat(sprintf("Size of array: %dx%d (%.2f%% != 0)\n",
                    nrow(object@assay), ncol(object@assay),
                    nnzero(object@assay) / length(object@assay) * 100L))
    }

    if (length(object@processing)) {
        cat("- - - Processing information - - -\n")
        cat(paste0(object@processing, collapse="\n"), "\n")
    }

    invisible(NULL)
})
