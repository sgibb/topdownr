#' @rdname FragmentViews-class
#' @param object FragmentViews
#' @aliases show,FragmentViews-method
#' @export
setMethod("show", "FragmentViews", function(object) {
    cat0(class(object), " on a ", length(subject(object)),
         "-letter sequence:\n")
    cat0("  ", .snippet(as.character(subject(object)),
                        getOption("width") - 2L), "\n")
    cat0("Views:\n")

    entries <- list(
        ## rownames
        c("", .hft(paste0("[", seq_len(length(object)), "]"), n=5L)),
        ## ranges
        start=c("start", .hft(start(object), n=5L)),
        end=c("end", .hft(end(object), n=5L)),
        width=c("width", .hft(width(object), n=5L)),
        ## metadata
        mass=c("mass", .hft(sprintf("%.2f", elementMetadata(object)$mass),
                            n=5L)),
        type=c("type", .hft(as.character(elementMetadata(object)$type), n=5L)),
        z=c("z", as.character(.hft(elementMetadata(object)$z, n=5L)))
    )

    cw <- .vapply1d(entries, function(e)max(nchar(e)))

    ## views
    mw <- sum(cw + 1L)
    ## n == 5 + 1 (to compensate for the fill ("...") below
    vw <- as(object[.hft(seq_len(length(object)), fill=NULL, n=6L)],
             "character")
    vw <- c("", .hft(.snippet(paste0("[", vw, "]"), getOption("width") - mw),
                     n=5L))

    ## format
    just <- c(rep("right", 5), "left", "right", "left")
    entries <- mapply(format, x=c(entries, list(vw)), justify=just)

    ## output
    out <- paste(apply(entries, 1L, paste, collapse=" "), collapse="\n")
    cat0(out, "\n")

    invisible(NULL)
})
