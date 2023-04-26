`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}
FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  slots <- na.omit(object = Filter(
    f = function(x) {
      return(class(x = slot(object = object, name = x)) == 'list')
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(class(x = object[[i]]))
    }
  )
  object.classes <- object.classes[object.classes %in% classes.keep]
  return(names(x = object.classes))
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

Sweep <- function(x, MARGIN, STATS, FUN = '-', check.margin = TRUE, ...) {
  if (any(grepl(pattern = 'X', x = names(x = formals(fun = sweep))))) {
    return(sweep(
      X = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  } else {
    return(sweep(
      x = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  }
}


PseudobulkExpression <- function(
  object,
  pb.method = 'average',
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = 'ident',
  add.ident = NULL,
  slot = 'data',
  verbose = TRUE,
  ...
) {
  CheckDots(..., fxns = 'CreateSeuratObject')
  if (!is.null(x = add.ident)) {
    .Deprecated(msg = "'add.ident' is a deprecated argument, please use the 'group.by' argument instead")
    group.by <- c('ident', add.ident)
  }
  if (!(pb.method %in% c('average', 'aggregate'))) {
    stop("'pb.method' must be either 'average' or 'aggregate'")
  }
  object.assays <- FilterObjects(object = object, classes.keep = 'Assay')
  assays <- assays %||% object.assays
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(x = assays) == 0) {
      stop("None of the requested assays are present in the object")
    } else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (length(x = slot) == 1) {
    slot <- rep_len(x = slot, length.out = length(x = assays))
  } else if (length(x = slot) != length(x = assays)) {
    stop("Number of slots provided does not match number of assays")
  }
  
  ## data: Group cell variable
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
  
  ## category.matrix: factor size of averaging expression, it's 1/cell_cluster_size
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
    category.matrix <- Matrix::sparse.model.matrix(object = as.formula(
      object = paste0(
        '~0+',
        paste0(
          "data[,",
          1:length(x = group.by),
          "]",
          collapse = ":"
        )
      )
    ))#
    colsums <- Matrix::colSums(x = category.matrix)
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
  for (i in 1:length(x = assays)) {
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
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
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
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    if (slot[1] == 'scale.data') {
      na.matrix <- data.return[[1]]
      na.matrix[1:length(x = na.matrix)] <- NA
      # TODO: restore once check.matrix is in SeuratObject
      # toRet <- CreateSeuratObject(
      #   counts = na.matrix,
      #   project = if (pb.method == "average") "Average" else "Aggregate",
      #   assay = names(x = data.return)[1],
      #   check.matrix = FALSE,
      #   ...
      # )
      toRet <- CreateSeuratObject(
        counts = na.matrix,
        project = if (pb.method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1],
        ...
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "counts",
        new.data = matrix()
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "data",
        new.data = na.matrix
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "scale.data",
        new.data = data.return[[1]]
      )
    } else {
      # TODO: restore once check.matrix is in SeuratObject
      # toRet <- CreateSeuratObject(
      #   counts = data.return[[1]],
      #   project = if (pb.method == "average") "Average" else "Aggregate",
      #   assay = names(x = data.return)[1],
      #   check.matrix = FALSE,
      #   ...
      # )
      toRet <- CreateSeuratObject(
        counts = data.return[[1]],
        project = if (pb.method == "average") "Average" else "Aggregate",
        assay = names(x = data.return)[1]
      )
      toRet <- SetAssayData(
        object = toRet,
        assay = names(x = data.return)[1],
        slot = "data",
        new.data = log1p(x = as.matrix(x = data.return[[1]]))
      )
    }
    #for multimodal data
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        if (slot[i] == 'scale.data') {
          na.matrix <- data.return[[i]]
          na.matrix[1:length(x = na.matrix)] <- NA
          # TODO: restore once check.matrix is in SeuratObject
          # toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = na.matrix, check.matrix = FALSE)
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = na.matrix)
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "counts",
            new.data = matrix()
          )
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "data",
            new.data = na.matrix
          )
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "scale.data",
            new.data = as.matrix(x = data.return[[i]])
          )
        } else {
          # TODO: restore once check.matrix is in SeuratObject
          # toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]], check.matrix = FALSE)
          toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
          toRet <- SetAssayData(
            object = toRet,
            assay = names(x = data.return)[i],
            slot = "data",
            new.data = log1p(x = as.matrix(x = data.return[[i]]))
          )
        }
        
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
      if (slot[which(DefaultAssay(object = object) %in% names(x = data.return))[1]] != 'scale.data') {
        toRet <- ScaleData(object = toRet, verbose = verbose)
      }
    }
    if ('ident' %in% group.by) {
      first.cells <- c()
      for (i in 1:ncol(x = category.matrix)) {
        first.cells <- c(first.cells, Position(x = category.matrix[,i], f = function(x) {x > 0}))
      }
      Idents(object = toRet) <- Idents(object = object)[first.cells]
    }
    return(toRet)
  } else {
    return(data.return)
  }
}

AverageExpression_v4 <- function(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = 'ident',
  add.ident = NULL,
  slot = 'data',
  verbose = TRUE,
  ...
) {
  return(
    PseudobulkExpression(
      object = object,
      pb.method = 'average',
      assays = assays,
      features = features,
      return.seurat = return.seurat,
      group.by = group.by,
      add.ident = add.ident,
      slot = slot,
      verbose = verbose,
      ...
    )
  )
}

CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

AverageExpression_v3 <- function (object, assays = NULL, features = NULL, return.seurat = FALSE, 
          add.ident = NULL, slot = "data", use.scale = FALSE, use.counts = FALSE, 
          verbose = TRUE, ...) 
{
  CheckDots(..., fxns = "CreateSeuratObject")
  if (use.scale) {
    .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
    slot <- "scale.data"
  }
  if (use.counts) {
    .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
    if (use.scale) {
      warning("Both 'use.scale' and 'use.counts' were set; using counts", 
              call. = FALSE, immediate. = TRUE)
    }
    slot <- "counts"
  }
  fxn.average <- switch(EXPR = slot, data = function(x) {
    return(mean(x = expm1(x = x)))
  }, mean)
  object.assays <- FilterObjects(object = object, classes.keep = "Assay")
  assays <- assays %||% object.assays
  ident.orig <- Idents(object = object)
  orig.levels <- levels(x = Idents(object = object))
  ident.new <- c()
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(assays) == 0) {
      stop("None of the requested assays are present in the object")
    }
    else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (!is.null(x = add.ident)) {
    new.data <- FetchData(object = object, vars = add.ident)
    new.ident <- paste(Idents(object)[rownames(x = new.data)], 
                       new.data[, 1], sep = "_")
    Idents(object, cells = rownames(new.data)) <- new.ident
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(object = object, assay = assays[i], 
                             slot = slot)
    features.assay <- features
    if (length(x = intersect(x = features, y = rownames(x = data.use))) < 
        1) {
      features.assay <- rownames(x = data.use)
    }
    data.all <- data.frame(row.names = features.assay)
    for (j in levels(x = Idents(object))) {
      temp.cells <- WhichCells(object = object, idents = j)
      features.assay <- unique(x = intersect(x = features.assay, 
                                             y = rownames(x = data.use)))
      if (length(x = temp.cells) == 1) {
        data.temp <- (data.use[features.assay, temp.cells])
        if (slot == "data") {
          data.temp <- expm1(x = data.temp)
        }
      }
      if (length(x = temp.cells) > 1) {
        data.temp <- apply(X = data.use[features.assay, 
                                        temp.cells, drop = FALSE], MARGIN = 1, FUN = fxn.average)
      }
      data.all <- cbind(data.all, data.temp)
      colnames(x = data.all)[ncol(x = data.all)] <- j
      if (verbose) {
        message(paste("Finished averaging", assays[i], 
                      "for cluster", j))
      }
      if (i == 1) {
        ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
      }
    }
    names(x = ident.new) <- levels(x = Idents(object))
    data.return[[i]] <- data.all
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(counts = data.return[[1]], 
                                project = "Average", assay = names(x = data.return)[1])
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
    }
    Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
    Idents(object = toRet) <- factor(x = Idents(object = toRet), 
                                     levels = as.character(x = orig.levels), ordered = TRUE)
    toRet <- NormalizeData(object = toRet, verbose = verbose)
    toRet <- ScaleData(object = toRet, verbose = verbose)
    return(toRet)
  }
  else {
    return(data.return)
  }
}