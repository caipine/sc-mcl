
get.data.dist <- function( object = object,
  pb.method = 'average',
  assays = "RNA",
  features = pth1,
  return.seurat = FALSE,
  group.by = "Finalname" , #'ident',
  add.ident = NULL,
  slot = 'data',
  verbose = TRUE
  ) {
  
  if (length(x = slot) == 1) {
    slot <- rep_len(x = slot, length.out = length(x = assays))
  } else if (length(x = slot) != length(x = assays)) {
    stop("Number of slots provided does not match number of assays")
  }
  data <- FetchData(object = object, vars = rev(x = group.by))
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
  if (nrow(x = data) < ncol(x = object)) {
    message("Removing cells with NA for 1 or more grouping variables")
    object <- subset(x = object, cells = rownames(x = data))
  }
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(x = data),
    FUN = function(i) {
      length(x = levels(x = data[, i]))
    }
  )
  if (any(num.levels == 1)) {
    message(paste0("The following grouping variables have 1 value and will be ignored: ",
                   paste0(colnames(x = data)[which(num.levels <= 1)], collapse = ", ")))
    group.by <- colnames(x = data)[which(num.levels > 1)]
    data <- data[, which(num.levels > 1), drop = F]
  }
  if (ncol(x = data) == 0) {
    message("All grouping variables have 1 value only. Computing across all cells.")
    category.matrix <- matrix(
      data = 1,
      nrow = ncol(x = object),
      dimnames = list(Cells(x = object), 'all')
    )
    if (pb.method == 'average') {
      category.matrix <- category.matrix / sum(category.matrix)
    }
  } else {
    category.matrix <- sparse.model.matrix(object = as.formula(
      object = paste0(
        '~0+',
        paste0(
          "data[,",
          1:length(x = group.by),
          "]",
          collapse = ":"
        )
      )
    ))
    colsums <- colSums(x = category.matrix)
    category.matrix <- category.matrix[, colsums > 0]
    colsums <- colsums[colsums > 0]
    if (pb.method == 'average') {
      category.matrix <- Sweep(
        x = category.matrix,
        MARGIN = 2,
        STATS = colsums,
        FUN = "/")
    }
    colnames(x = category.matrix) <- sapply(
      X = colnames(x = category.matrix),
      FUN = function(name) {
        name <- gsub(pattern = "data\\[, [1-9]*\\]", replacement = "", x = name)
        return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))), collapse = "_"))
      })
  }
  data.return <- list()
#  for (i in 1:length(x = assays)) {
i= 1
    data.use <- GetAssayData(
      object = object,
      assay = assays[i],
      slot = slot[i]
    )
    features.to.avg <- features %||% rownames(x = data.use)
    if (inherits(x = features, what = "list")) {
      features.to.avg <- features[i]
    }
    if (IsMatrixEmpty(x = data.use)) {
      warning(
        "The ", slot[i], " slot for the ", assays[i],
        " assay is empty. Skipping assay.", immediate. = TRUE, call. = FALSE)
      next
    }
    bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = bad.features) > 0) {
      warning(
        "The following ", length(x = bad.features),
        " features were not found in the ", assays[i], " assay: ",
        paste(bad.features, collapse = ", "), call. = FALSE, immediate. = TRUE)
    }
    features.assay <-   features.to.avg   ############################            intersect(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = features.assay) > 0) {
      data.use <- data.use[features.assay, ]
    } else {
      warning("None of the features specified were found in the ", assays[i],
              " assay.", call. = FALSE, immediate. = TRUE)
      next
    }
    if (slot[i] == 'data') {
      data.use <- expm1(x = data.use)
      if (any(data.use == Inf)) {
        warning("Exponentiation yielded infinite values. `data` may not be log-normed.")
      }
    }
data.return[[i]] <- as.matrix(x = (data.use %*% category.matrix))
data.return <- data.return[[i]]
return(data.return)
}

