#' @rdname exampleESE
#' @name exampleESE
#' @aliases exampleESE
#'
#' @title Example EnrichmentSE object
#'
#' @description
#' Small sample signal from ENCODE ChIP-seq for H3K27ac, H3K4me3 and p300, 
#' around some p300 binding sites and TSS in mESC.
#'
#' @return a named character vector of length 1.
NULL

#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
.EnrichmentSE <- setClass("EnrichmentSE",
  contains="RangedSummarizedExperiment",
  representation(
    rowRanges="GenomicRanges"
  ),
  prototype(
    rowRanges=GRanges()
  )
)

EnrichmentSE <- function(assays, rowRanges=NULL, ...){
  if(is.null(rowRanges)){
    a <- SummarizedExperiment(assays, ...)
  }else{
    a <- SummarizedExperiment(assays, rowRanges=rowRanges, ...)
  }
  a@NAMES <- NULL
  epiwraps:::.EnrichmentSE(a)
}



#' @export
setMethod("[", c("EnrichmentSE", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE){
  if(missing(j)) j <- NULL
  if(missing(i)) i <- NULL
  .resizeESE(x, i, j, drop=drop)
})

#' @importFrom SummarizedExperiment assays assays<- rowRanges 
#' @importFrom SummarizedExperiment assayNames assayNames<-
setMethod("show", "EnrichmentSE", function(object){
  cat("class:", class(object), "\n")
  cat(ncol(object), "tracks across", nrow(object), "regions\n")
  
  ## assays()
  nms <- assayNames(object)
  if (is.null(nms))
    nms <- character(length(assays(object, withDimnames=FALSE)))
  coolcat("assays(%d): %s\n", nms)
  
  ## rownames()
  rownames <- rownames(object)
  if (!is.null(rownames)) coolcat("rownames(%d): %s\n", rownames)
  else cat("rownames: NULL\n")
  
  ## rowData`()
  coolcat("rowData names(%d): %s\n", names(rowData(object, use.names=FALSE)))
  
  ## colnames()
  colnames <- colnames(object)
  if (!is.null(colnames)) coolcat("colnames(%d): %s\n", colnames)
  else cat("colnames: NULL\n")
  
  ## colData()
  coolcat("colData names(%d): %s\n", names(colData(object)))
  
  ## metadata()
  expt <- names(metadata(object))
  if (is.null(expt))
    expt <- character(length(metadata(object)))
  coolcat("metadata(%d): %s\n", expt)
  
})

#' showTrackInfo
#' 
#' Provide some information about the relative signal ranges of each track.
#'
#' @param x A named list of signal matrices or an EnrichmentSE object as 
#'   produced by \code{\link{signal2Matrix}}
#' @param assay The assay to use, defaults to the input assay.
#' @param doPrint Logical; whether to print the information.
#'
#' @return An invisible list of captions.
#' @export
#' @examples
#' data(exampleESE)
#' showTrackInfo(exampleESE)
showTrackInfo <- function(x, assay="input", doPrint=TRUE){
  if(is(x, "EnrichmentSE")) x <- getSignalMatrices(x, assay=assay)
  stopifnot(is(x,"list") & all(vapply(x, FUN.VALUE=logical(1), 
                                      class2="normalizedMatrix", FUN=is)))
  out <- lapply(x, FUN=function(x){
    upstream_index = attr(x, "upstream_index")
    downstream_index = attr(x, "upstream_index")
    target_index_len = length(attr(x, "target_index"))
    extend = attr(x, "extend")
    out <- paste0("  -",extend[1],"/+",extend[2],"bp (")
    if(length(upstream_index)==length(downstream_index)){
      out <- paste0(out, length(upstream_index), " windows each")
    }else{
      out <- paste0(out, "respectively ", length(upstream_index), " and ",
                    length(downstream_index), " windows")
    }
    if(isTRUE(attr(x, "smooth"))) out <- paste0(out, ", smoothed")
    if(target_index_len>1){
      out <- paste0(out, ")\n  around given regions (", target_index_len, 
                    " windows)")
    }else{
      out <- paste0(out, ")\n  around the centers of given regions")
    }
    if(doPrint){
      cat(attr(x, "signal_name"),"(",paste(dim(x), collapse="x"),") :\n")
      cat(out, "\n")
    }
    out
  })
  invisible(out)
}

#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
.ml2assay <- function(ml){
  ml <- .comparableMatrices(ml)
  a <- DataFrame(row.names=row.names(ml[[1]]))
  for(f in names(ml)) a[[f]] <- ml[[f]]
  a
}

#' Creates an EnrichmentSE from a list of normalizedMatrix objects
#' 
#' @param ml A named list of normalizedMatrix objects with corresponding rows.
#' @param rowRanges An optional GRanges object corresponding to the rows of each
#'   object of `ml`.
#' @param assayName The name of the assay, defaults to 'input'
#' @param addScore Logical; whether to add an enriched_score assay.
#' @param ... Passed to `SummarizedExperiment()`
#'
#' @importFrom SummarizedExperiment assays<- assay<- assays
#' @export
#' @examples
#' # for an example we first need a list of signal matrices. To this end, 
#' # we first fetch the path to the example bigwig file:
#' bw <- system.file("extdata/example_atac.bw", package="epiwraps")
#' # we load example regions:
#' regions <- rtracklayer::import(system.file("extdata/example_peaks.bed", 
#'                                            package="epiwraps"))
#' # we obtain the matrix of the signal around the regions, indicating that we
#' # want the output as a list of signal matrices:
#' m <- signal2Matrix(bw, regions, ret="list")
#' # we can then transform this into an EnrichmentSE object:
#' m <- ml2ESE(m)
ml2ESE <- function(ml, rowRanges, assayName="input", addScore=FALSE, ...){
  a <- .ml2assay(ml)
  stopifnot(is(rowRanges,"GRanges"))
  if(length(rowRanges)!=nrow(a)){
    if(is.null(names(rowRanges)) || !all(row.names(a) %in% names(rowRanges)))
      stop("length(rowRanges) is different from the number of signal rows,",
           " and there are no matching names!")
    rowRanges <- rowRanges[row.names(a)]
  }
  fnames <- row.names(a)
  if(!is.null(rowRanges)){
    if(!is.null(names(rowRanges)) && !is.null(fnames)){
      if(!all(fnames %in% names(rowRanges))){
        if(length(fnames)==length(rowRanges)){
          warning("The names of the features do not match those of 'rowRanges'.",
                  "The names will be disregarded and their orders assumed to ",
                  "be the same")
        }else{
          stop("The matrices' features/rows do not match the 'rowRanges'.")
        }
      }else{
        rowRanges <- rowRanges[fnames]
      }
    }else{
      if(length(fnames)!=length(rowRanges))
        stop("The matrices' features/rows do not match the 'rowRanges'.")
    }
  }
  a <- EnrichmentSE(list(input=a), rowRanges=rowRanges, ...)
  if(addScore){
    es <- DataFrame(lapply(ml, enriched_score))
    dimnames(es) <- dimnames(a)
    suppressWarnings({assays(a)$enriched_score <- es})
  }
  suppressWarnings({assayNames(a)[1] <- assayName})
  a
}

.ese2ml <- function(x, assay=1L){
  if(is.numeric(assay) && length(setdiff(assayNames(x),"enriched_score"))>1 && 
     !is.null(assayNames(x)[assay]) && assayNames(x)[assay]!=""){
    message("Using assay ", assayNames(x)[assay])
  }
  x <- assays(x)[[assay]]
  stopifnot(is(x,"DFrame"))
  x <- as.list(x)
  stopifnot(all(unlist(lapply(x, class2="normalizedMatrix", is))))
  x
}

#' getSignalMatrices
#' 
#' Extracts a list of signal matrices from an EnrichmentSE object.
#'
#' @param x An object of class `EnrichmentSE`, as produced by 
#'   \code{\link{signal2Matrix}}.
#' @param assay The assay to extract (defaults to the first assay).
#'
#' @return A list of normalizedMatrix objects.
#' @export
#' @examples
#' # we first get an EnrichmentSE object:
#' data(exampleESE)
#' # then we can extract the list of signal matrices:
#' sm <- getSignalMatrices(exampleESE)
getSignalMatrices <- function(x, assay=1L){
  stopifnot(is(x,"EnrichmentSE"))
  .ese2ml(x, assay=assay)
}


.parseEseVar <- function(se, x, type=c("row","col")){
  if(is.null(x)) return(NULL)
  type <- match.arg(type)
  if(is.vector(x) & length(x)>1){
    if(length(x)!=ifelse(type=="row",nrow(se),ncol(se)))
      stop("The vector provided for '", deparse(substitute(x)), "' does not",
           " have the right length.")
    return(x)
  }
  stopifnot(is(se,"EnrichmentSE"))
  if(type=="row"){
    x <- rowData(se)[[x]]
  }else if(type=="col"){
    x <- colData(se)[[x]]
  }
  if(is.null(x)) stop("'",deparse(substitute(x)), "' not found in ",type,"Data")
  x
}

#' addAssayToESE
#' 
#' Adds an assay of signal matrices to an existing `EnrichmentSE` object.
#'
#' @param x An object of class `EnrichmentSE`, as produced by 
#'   \code{\link{signal2Matrix}}.
#' @param a The assay to add, e.g. a list of normalizedMatrix objects
#' @param name 
#' @param replace Logical, whether to replace any existing assay of the same 
#'   name (default TRUE). If FALSE and the assay already existed, the new assay 
#'   name is given a suffix.
#'
#' @return `x` with the added/updated assay.
#' @export
#' @examples
#' # we first get an EnrichmentSE object:
#' data(exampleESE)
#' # then we will create a new assay which is simply sqrt-transformed, and add 
#' # it back in the object
#' newAssay <- lapply(getSignalMatrices(x), sqrt)
#' exampleESE <- addAssayToESE(exampleESE, newAssay, named="sqrt")
addAssayToESE <- function(x, a, name="normalized", replace=TRUE){
  stopifnot(is(x,"EnrichmentSE"))
  if(is(a, "list")) a <- .ml2assay(a)
  if(ncol(a)!=ncol(x))
    stop("The new assay does not have as many tracks as the existing ones.")
  al <- list()
  if(!replace) name <- rev(make.unique(c(assayNames(x), name)))[1]
  al[[name]] <- a
  assays(x) <- c(al, as.list(assays(x)[setdiff(assayNames(x),name)]))
  x
}

# called by the method
.resizeESE <- function(x, i=NULL, j=NULL, drop=FALSE){
  stopifnot(is(x,"EnrichmentSE"))
  rr <- rowRanges(x)
  cd <- colData(x)
  al <- assays(x)
  if(!is.null(i)){
    al <- lapply(al, FUN=function(x){
      if(is(x[,1],"normalizedMatrix")){
        x <- .ml2assay(lapply(x, FUN=function(y) .resizeNmatrix(y,i)))
      }else{
        x <- x[i,,drop=FALSE]
      }
    })
    rr <- rr[i]
  }
  if(!is.null(j)){
    al <- lapply(al, FUN=function(x) x[,j,drop=FALSE])
    cd <- cd[j,]
  }
  EnrichmentSE(assays=al, rowRanges=rr, colData=cd, metadata=metadata(x))
}


#' @export
setMethod("score", "EnrichmentSE", function(x,...){
  sapply(getSignalMatrices(x), FUN=EnrichedHeatmap::enriched_score)
})
