#' @rdname FragmentViews-class
#' @param object,x,y FragmentViews
#' @aliases combine,FragmentViews-method
#' @export
setMethod("combine", signature(x="FragmentViews", y="FragmentViews"),
          function(x, y) {
    if (!identical(x@subject, y@subject)) {
        stop("'subject' must be identical")
    }
    if (!identical(x@metadata$mass, y@metadata$mass)) {
        stop("'metadata(...)$mass' must be identical")
    }
    if (!identical(x@metadata$modifications, y@metadata$modifications)) {
        stop("'metadata(...)$modifications' must be identical")
    }
    un <- union(names(x), names(y))
    xn <- intersect(names(x), un)
    yn <- intersect(names(y), un)
    x@elementMetadata <- rbind(elementMetadata(x[xn]), elementMetadata(y[yn]))
    x@ranges <- c(x@ranges[xn], y@ranges[yn])
    x@metadata <- modifyList(metadata(x), metadata(y))
    x <- x[order(elementMetadata(x)$mass)]
    if (validObject(x)) {
        x
    }
})

#' @rdname FragmentViews-class
## @param object FragmentViews
#' @aliases show,FragmentViews-method
#' @export
setMethod("show", "FragmentViews", function(object) {
    cat0(class(object), " on a ",
        length(subject(object)), "-letter sequence:\n")
    cat0("  ",
        .snippet(
            as.character(subject(object)),
            getOption("width") - 2L
        ), "\n"
    )

    if (length(metadata(object)$mass)) {
        cat0("Mass:\n")
        cat0(paste0("  ", metadata(object)$mass, "\n"))
    }

    if (length(metadata(object)$modifications)) {
        cat0("Modifications:\n")
        cat0(paste0("  ", metadata(object)$modifications, "\n"))
    }

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
        name=c("name", .hft(names(object), n=5L)),
        type=c("type", .hft(as.character(elementMetadata(object)$type), n=5L)),
        z=c("z", as.character(.hft(elementMetadata(object)$z, n=5L)))
    )

    cw <- .vapply1d(entries, function(e)max(nchar(e)))

    ## views
    mw <- sum(cw + 1L)
    ## n == 5 + 1 (to compensate for the fill ("...") below
    vw <-
        as(object[.hft(seq_len(length(object)), fill=NULL, n=6L)], "character")
    vw <- c("",
        .hft(.snippet(paste0("[", vw, "]"), getOption("width") - mw), n=5L))

    ## format
    just <- c(rep("right", 5), "left", "left", "right", "left")
    entries <- mapply(format, x=c(entries, list(vw)), justify=just)

    ## output
    out <- paste(apply(entries, 1L, paste, collapse=" "), collapse="\n")
    cat0(out, "\n")

    invisible(NULL)
})

#' @rdname FragmentViews-class
#' @name coerce,FragmentViews,data.frame-method
#' @section Coercion:
#'
#' `as(object, "data.frame")`: Coerce an
#' [FragmentViews-class] object into a `data.frame`.
#' @examples
#' as(fv, "data.frame")
setAs("FragmentViews", "data.frame", function(from) {
    data.frame(
        fragment=as(from, "character"),
        start=start(from),
        end=end(from),
        width=width(from),
        name=names(from),
        type=elementMetadata(from)$type,
        mass=elementMetadata(from)$mass,
        z=elementMetadata(from)$z,
        row.names=names(from)
    )
})
