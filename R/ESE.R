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
  .EnrichmentSE(a)
}



#' @export
setMethod("[", "EnrichmentSE", function(x, i, j, drop=TRUE){
  .resizeESE(x, i, j, drop=drop)
})

#' @importFrom SummarizedExperiment assays assays<- rowRanges
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
#' @param assayName The name of the assay, defaults to 'input'
#' @param rowRanges An optional GRanges object corresponding to the rows of each
#'   object of `ml`.
#' @param addScore Logical; whether to add an enriched_score assay.
#' @param ... Passed to `SummarizedExperiment()`
#'
#' @importFrom SummarizedExperiment assays<- assay<- assays
#' @export
ml2ESE <- function(ml, assayName="input", rowRanges=NULL, addScore=FALSE, ...){
  a <- .ml2assay(ml)
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
  al <- list()
  al[[assayName]] <- a
  a <- EnrichmentSE(a, rowRanges=rowRanges, ...)
  if(addScore){
    es <- DataFrame(lapply(ml, enriched_score))
    dimnames(es) <- dimnames(a)
    suppressWarnings({assays(a)$enriched_score <- es})
  }
  a
}

.ese2ml <- function(x, assay=1L){
  if(is.numeric(assay) && length(assays(x))>1 && 
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

.addAssayToESE <- function(ml, a, name="normalized", replace=FALSE){
  if(is(a, "list")) a <- .ml2assay(a)
  al <- list()
  if(replace) name <- rev(make.unique(c(assayNames(ml), name)))[1]
  al[[name]] <- a
  assays(ml) <- c(al, as.list(assays(ml)[setdiff(assayNames(ml),name)]))
  ml
}

.resizeESE <- function(x, i=NULL, j=NULL, drop=FALSE){
  stopifnot(is(x,"EnrichmentSE"))
  rr <- rowRanges(x)
  cd <- colData(x)
  al <- assays(x)
  if(!is.null(i)){
    al <- lapply(al, FUN=function(x){
      if(is(x[,1],"normalizedMatrix")){
        x <- DataFrame(lapply(x, FUN=function(y) .resizeNmatrix(y,i)))
      }else{
        x <- x[i,]
      }
    })
    rr <- rr[i]
  }
  if(!is.null(j)){
    al <- lapply(al, FUN=function(x) x[,j,drop=FALSE])
    cd <- cd[j,]
  }
  EnrichmentSE(SummarizedExperiment(assays=al, rowRanges=rr, colData=cd,
                                    metadata=metadata(x)))
}
