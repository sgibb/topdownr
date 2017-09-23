#' Find best combination of conditions for highest coverage
#'
#' @param object `NCBSet`
#' @param n `integer`, max number of combinations/iterations.
#' @param minN `integer` stop if there are less than `minN` additional fragments
#' added.
#' @return `matrix`, first column: index of condition, second column: number of
#' newly added fragments third column: number of newly covered bonds.
#' @export
#' @noRd
setMethod("bestConditions", "NCBSet",
          function(object, n=ncol(object), minN=0L) {
    m <- .bestNcbCoverageCombination(object@assay, n=n, minN=minN)
    rownames(m) <- colnames(object)[m[, "index"]]
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
          function(object, nCombinations=10, cumCoverage=TRUE,
                   labels=colnames(object), ...) {

    if (length(labels) != ncol(object)) {
        stop("The number of elements in 'labels' has to be the same as ",
             "the number of columns in 'object'.")
    }
    d <- .dgcMatrix2data.frame(object@assay)
    d$Activation <- as.character(object$Activation[d$col])

    if (nCombinations) {
        i <- bestConditions(object, nCombinations)[, "index"]
        combinations <- object[, i]
        labels <- c(labels, labels[i])

        if (cumCoverage) {
            cmb <- .dgcMatrix2data.frame(.cumComb(combinations@assay))
        } else {
            cmb <- .dgcMatrix2data.frame(combinations@assay)
        }
        cmb$col <- cmb$col + ncol(object)
        cmb$Activation <- "Cmb"
        d <- rbind(d, cmb)
    }

    ## tell ggplot that we have an already ordered factor
    d$col <- as.ordered(d$col)
    d$x <- as.factor(d$x)
    d$Activation <- factor(d$Activation,
                           levels=c("CID", "HCD", "ETD", "ETcid", "EThcd",
                                    "Cmb"),
                           ordered=TRUE)

    names(labels) <- seq_along(labels)

    ggplot(data=d) +
        geom_raster(aes_string(x="col", y="row", fill="x")) +
        facet_grid(. ~ Activation, scales="free_x", space="free_x") +
        scale_fill_manual(name="Observed Fragments",
                          labels=c("N-terminal", "C-terminal", "Both"),
                          values=c("#1b9e77", "#d95f02", "#7570b3")) +
        scale_x_discrete(name="Condition", expand=c(0L, 0L),
                         labels=labels) +
        scale_y_continuous(name=expression(
                            paste("Bonds (", N-terminal %->% C-termnial, ")")),
                           breaks=seq_len(nrow(object)),
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

#' @param object `NCBSet`
#' @param what `character`, summarise by "conditions" (col) or "fragments" (rows)
#' @return `data.frame`
#' @export
#' @noRd
setMethod("summary", "NCBSet",
          function(object, what=c("conditions", "bonds"), ...) {
    what <- if (match.arg(what) == "conditions") { "columns" } else { "rows" }
    callNextMethod(object=object, what=what)
})

