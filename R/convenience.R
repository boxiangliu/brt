#' @description Return the basename of a file name without path or extension.
basename.2 <- function(names, extension){
  base.extension = basename(names)
  base = stringr::str_replace(base.extension, extension, "")
  return(base)
}


#' @description Display correlation on pairs() function.
#' @title Display correlation on pairs() function.
#' @seealso \code{pairs()} function help where this function is scraped.

panel.cor <- function(x, y, digits = 4, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#' Reads all gene count files in a directory and creates a count data.frame from them. Each column of the data.frame represents a file; each row represent a gene. Empty rows are removed by default.
#' @name read.counts
#' @return A count matrix in class data.frame
#' @examples # If I want to generate a count matrix using count files with extension
#' # "Aligned.out.count" in directory "~/RNAseq" I would:
#' read.counts(path = '~/RNAseq', extension = 'Aligned.out.count')
#' @title Generate a Count Matrix
read.counts <- function(path = '.', extension = ".Aligned.out.count", remove.empty.rows = TRUE, log = FALSE){
  pattern = paste0("*", extension)
  filenames = list.files(path = path, pattern = pattern ,full.names = TRUE)
  files = lapply(filenames, read.table, header = FALSE)
  num_files = length(filenames)
  num_genes = dim(files[1][[1]])[1]
  rownames = files[1][[1]][,1]
  colnames = basename.2(filenames, extension)
  counts = matrix(0, nrow = num_genes, ncol = num_files, dimnames = list(rownames, colnames))
  for (i in seq(1,num_files)){
    if (dim(counts)[1] != dim(files[i][[1]])[1]) {stop("count matrix and count file dimensions do not agree")}
    counts[,i] = files[i][[1]][,2]
  }
  if (remove.empty.rows) {counts = counts[rowSums(counts) != 0, ]} # remove empty rows
  counts = counts[-((nrow(counts)-2):nrow(counts)),] # remove extra lines
  colnames(counts) = make.names(colnames(counts))
  counts = as.data.frame(counts)
  if (log) {counts = log10(counts + 1)}
  return(counts)
}

#' @title Merge Technical Duplicates in a Count Matrix
#' @description Merge the technical duplicate that has the same sample name but different extensions
#' @seealso \code{read.counts()}
#' @examples # merging the 4 lanes from NextSeq RNAseq data
#' counts.merged = merge.counts(count.matrix, method = 'average', extension = '_L00[1234]')

merge.counts <-function(count.matrix, method = 'average', extension = '_L00[1234]', log = FALSE){
  sample.names = unique(basename.2(colnames(counts), '_L00[1234]'))
  if (!is.data.frame(counts)) {counts = as.data.frame(counts)}
  counts.merged = as.data.frame(matrix(nrow = nrow(counts)))
  if (method == 'average') {
    for (sample.name in sample.names){
      counts.merged[,sample.name] = rowMeans(counts[,which(stringr::str_match(colnames(counts), sample.name) == sample.name)])
    }
  } else if (method == 'sum') {
    for (sample.name in sample.names){
      counts.merged[,sample.name] = rowSums(counts[,which(stringr::str_match(colnames(counts), sample.name) == sample.name)])
    }
  } else {
    stop("Unknown merge method. Use 'average' or 'sum'.")
  }

  counts.merged = counts.merged[which(colnames(counts.merged) != 'V1')]
  rownames(counts.merged) = rownames(counts)
  if(log) {counts.merged = log10(counts.merged + 1)}
  return(counts.merged)
}
