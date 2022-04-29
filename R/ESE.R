

.ml2assay <- function(ml){
  ml <- .comparableMatrices(ml)
  a <- DataFrame(row.names=row.names(ml[[1]]))
  for(f in names(ml)) a[[f]] <- ml[[f]]
  a
}

.ml2ese <- function(ml, assayName="density", rowRanges=NULL, ...){
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
  if(is.null(rowRanges)){
    a <- SummarizedExperiment(al, ...)
  }else{
    a <- SummarizedExperiment(al, rowRanges=rowRanges, ...)
  }
  a
}

.ese2ml <- function(x, assay=1){
  x <- assays(x)[[assay]]
  stopifnot(is(x,"DFrame"))
  x <- as.list(x)
  stopifnot(all(unlist(lapply(x, class2="normalizedMatrix", is))))
  x
}
