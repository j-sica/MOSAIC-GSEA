options(warn = 1)  # Print each warning immediately so it is preserved in the stage log.
# Fit model and build contrasts

counts_file     <- "all_counts_final.txt"
metadata_file   <- "metadata_ethanol_inter.tsv"
formula_str     <- "~ strain * treatment * oxygen"
contrast_factor <- "treatment"
reference_level <- "YPD"
output_dir      <- "results/"

filter_column       <- NULL
filter_value        <- NULL
variable_types_file <- NULL

interaction_factors <- NULL
sig_threshold <- 0.05
cores <- 4L


suppressPackageStartupMessages({
  library(limma); library(edgeR); library(tools)
  library(variancePartition); library(BiocParallel)
})
# Respect small machines while keeping four workers as the default.
available_cores <- suppressWarnings(parallel::detectCores(logical=TRUE))
if (length(available_cores) == 1L && is.finite(available_cores) && available_cores >= 1L && cores > available_cores) {
  warning(sprintf("Requested %d workers, but only %d logical core(s) detected; using %d worker(s).",
                  cores, available_cores, available_cores), call.=FALSE)
  cores <- as.integer(available_cores)
}
has_random <- grepl("\\|", formula_str)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_auto <- function(path) {
  ext <- tolower(file_ext(path))
  if (ext == "csv") {
    read.csv(path, header=TRUE, check.names=FALSE,
             stringsAsFactors=FALSE, comment.char="")
  } else {
    read.table(path, header=TRUE, sep="\t", check.names=FALSE,
               stringsAsFactors=FALSE, quote="", comment.char="")
  }
}

read_metadata <- function(path) {
  ext <- tolower(file_ext(path))
  if (ext == "csv") {
    read.csv(path, header=TRUE, check.names=FALSE, colClasses="character",
             stringsAsFactors=FALSE, comment.char="", na.strings=c("", "NA", "NaN"))
  } else {
    read.table(path, header=TRUE, sep="\t", check.names=FALSE, colClasses="character",
               stringsAsFactors=FALSE, quote="", comment.char="", na.strings=c("", "NA", "NaN"))
  }
}

read_variable_types <- function(path, columns) {
  if (is.null(path) || !nzchar(path)) return(character(0))
  lines <- trimws(readLines(path, warn=FALSE, encoding="UTF-8"))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) stop("The variable-types file is empty.")
  parts <- strsplit(lines, "[\\t ]+", perl=TRUE)
  if (any(lengths(parts) < 2L))
    stop("Each variable-types row must contain a column name and factor or number.")
  tab <- data.frame(
    variable=vapply(parts, `[`, character(1), 1),
    type=tolower(vapply(parts, `[`, character(1), 2)),
    stringsAsFactors=FALSE
  )
  if (tolower(tab$variable[1]) %in% c("variable", "column", "name") &&
      tab$type[1] %in% c("type", "class")) tab <- tab[-1, , drop=FALSE]
  tab$type[tab$type == "numeric"] <- "number"
  if (any(!tab$type %in% c("factor", "number")))
    stop("Variable types must be factor or number.")
  if (anyDuplicated(tab$variable))
    stop("The variable-types file contains duplicate column names.")
  unknown <- setdiff(tab$variable, columns)
  if (length(unknown))
    stop(sprintf("Unknown column(s) in variable-types file: %s", paste(unknown, collapse=", ")))
  setNames(tab$type, tab$variable)
}

infer_variable_type <- function(x) {
  z <- trimws(as.character(x[!is.na(x)]))
  if (!length(z)) return("factor")
  numeric_text <- grepl("^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$", z, perl=TRUE)
  if (all(numeric_text) && any(grepl("\\.", z))) "number" else "factor"
}

apply_variable_types <- function(meta, path=NULL) {
  overrides <- read_variable_types(path, colnames(meta))
  rows <- vector("list", length(meta))
  for (i in seq_along(meta)) {
    nm <- colnames(meta)[i]
    type <- if (nm %in% names(overrides)) overrides[[nm]] else infer_variable_type(meta[[nm]])
    source <- if (nm %in% names(overrides)) "file" else "auto"
    if (type == "number") {
      old_na <- is.na(meta[[nm]])
      value <- suppressWarnings(as.numeric(meta[[nm]]))
      bad <- !old_na & is.na(value)
      if (any(bad))
        stop(sprintf("Column '%s' was set to number but contains non-numeric values.", nm))
      meta[[nm]] <- value
    } else {
      vals <- as.character(meta[[nm]])
      lev <- unique(vals[!is.na(vals)])
      meta[[nm]] <- factor(vals, levels=lev)
    }
    rows[[i]] <- data.frame(variable=nm, type=type, source=source, stringsAsFactors=FALSE)
  }
  list(data=meta, types=do.call(rbind, rows))
}

expand_formula_display <- function(fstr) {
  f <- gsub("\\s+", "", fstr)
  random_part <- regmatches(f, gregexpr("\\+?\\([^)]*\\|[^)]*\\)", f, perl=TRUE))[[1]]
  main_f <- Reduce(function(s, r) gsub(r, "", s, fixed=TRUE), random_part, f)
  main_f <- gsub("~", "", main_f); main_f <- trimws(gsub("\\++", "+", main_f))
  
  has_star <- grepl("\\*", main_f)
  if (!has_star) return(fstr)
  
  terms <- strsplit(main_f, "\\+")[[1]]
  expanded <- character(0)
  for (term in terms) {
    if (!grepl("\\*", term)) {
      expanded <- c(expanded, term)
    } else {
      parts <- strsplit(term, "\\*")[[1]]
      n <- length(parts)
      for (k in 1:n) {
        for (combo in combn(parts, k, simplify=FALSE)) {
          new_term <- paste(combo, collapse=":")
          if (!new_term %in% expanded) expanded <- c(expanded, new_term)
        }
      }
    }
  }
  random_str <- if (length(random_part) > 0) paste0("+", paste(gsub("^\\+","",random_part), collapse="+")) else ""
  paste0("~ ", paste(expanded, collapse=" + "), random_str)
}

expanded_formula_str <- expand_formula_display(formula_str)
cat(sprintf("Model: %s\n", formula_str))

auto_detect_interaction_factors <- function(fstr, cf) {
  f <- gsub("\\s+", "\\s*", fstr)
  f_clean <- gsub("\\s+", "", fstr)
  f_clean <- gsub("\\([^)]*\\|[^)]*\\)", "", f_clean)
  f_clean <- gsub("~", "", f_clean)
  
  found <- character(0)
  
  star_groups <- regmatches(f_clean, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*(?:\\*[A-Za-z_.][A-Za-z0-9_.]*)+",
                                              f_clean, perl=TRUE))[[1]]
  for (grp in star_groups) {
    parts <- strsplit(grp, "\\*")[[1]]
    if (cf %in% parts) found <- union(found, setdiff(parts, cf))
  }
  
  colon_terms <- regmatches(f_clean, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*(?::[A-Za-z_.][A-Za-z0-9_.]*)+",
                                              f_clean, perl=TRUE))[[1]]
  for (term in colon_terms) {
    parts <- strsplit(term, ":")[[1]]
    if (cf %in% parts) found <- union(found, setdiff(parts, cf))
  }
  
  found
}

formula_has_interaction <- function(fstr, vars) {
  fixed <- trimws(gsub("\\+?\\s*\\([^)]*\\|[^)]*\\)", "", fstr))
  labels <- attr(terms(as.formula(fixed)), "term.labels")
  if (!length(labels)) return(FALSE)
  any(vapply(strsplit(labels, ":", fixed=TRUE), function(parts) setequal(parts, vars), logical(1)))
}

counts_raw <- read_auto(counts_file)
gene_col   <- colnames(counts_raw)[1]

# Keep duplicate input rows distinct for DE and volcano plots; annotation uses original_gene.
original_gene_ids <- trimws(as.character(counts_raw[[gene_col]]))
if (anyNA(original_gene_ids) || any(!nzchar(original_gene_ids)))
  stop("The first count-matrix column contains missing or empty gene identifiers.")

gene_group_n <- ave(seq_along(original_gene_ids), original_gene_ids,
                    FUN=function(ii) length(ii))
gene_group_i <- ave(seq_along(original_gene_ids), original_gene_ids,
                    FUN=seq_along)
internal_gene_ids <- ifelse(
  gene_group_n > 1L,
  paste0(original_gene_ids, " (", gene_group_i, ")"),
  original_gene_ids
)
if (anyDuplicated(internal_gene_ids)) {
  stop(paste0(
    "Gene identifiers could not be made unambiguous with the ' (n)' suffix convention. ",
    "This can happen when an input identifier already has the same text as a generated duplicate label; ",
    "use unique biological identifiers in the count matrix."
  ))
}
original_gene_by_internal <- setNames(original_gene_ids, internal_gene_ids)

gene_identifier_map <- data.frame(
  gene=internal_gene_ids,
  original_gene=original_gene_ids,
  duplicate_input_id=(gene_group_n > 1L),
  duplicate_index=ifelse(gene_group_n > 1L, gene_group_i, NA_integer_),
  stringsAsFactors=FALSE,
  check.names=FALSE
)
write.table(gene_identifier_map, file.path(output_dir, "gene_identifier_map.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
if (any(gene_group_n > 1L)) {
  cat(sprintf(
    "Duplicate input gene labels: %d row(s) across %d duplicated identifier(s); retained separately as 'GENE (1)', 'GENE (2)', ... for DE/reporting.\n\n",
    sum(gene_group_n > 1L), length(unique(original_gene_ids[gene_group_n > 1L]))
  ))
}

rownames(counts_raw) <- internal_gene_ids
counts_raw <- counts_raw[, colnames(counts_raw) != gene_col, drop=FALSE]
counts_mat <- as.matrix(counts_raw); storage.mode(counts_mat) <- "numeric"
cat(sprintf("Input: %d samples | %d genes\n", ncol(counts_mat), nrow(counts_mat)))

meta <- read_metadata(metadata_file)
rownames(meta) <- meta[[1]]; meta <- meta[,-1,drop=FALSE]

if (anyDuplicated(rownames(counts_mat)))
  stop("Internal gene identifiers are unexpectedly duplicated after duplicate-label resolution.")
if (anyDuplicated(colnames(counts_mat)))
  stop("Duplicate sample identifiers found in the count matrix.")
if (anyDuplicated(rownames(meta)))
  stop("Duplicate sample identifiers found in metadata.")
if (anyNA(counts_mat) || any(!is.finite(counts_mat)))
  stop("Count matrix contains missing or non-finite values.")
if (any(counts_mat < 0)) stop("Count matrix contains negative values.")
if (any(abs(counts_mat - round(counts_mat)) > 1e-8))
  stop("Count matrix must contain integer counts.")

missing_meta <- setdiff(colnames(counts_mat), rownames(meta))
missing_counts <- setdiff(rownames(meta), colnames(counts_mat))
if (length(missing_meta) || length(missing_counts)) {
  details <- c(
    if (length(missing_meta))
      paste0("Samples present in the count matrix but absent from metadata: ",
             paste(missing_meta, collapse=", "), ".") else NULL,
    if (length(missing_counts))
      paste0("Samples present in metadata but absent from the count matrix: ",
             paste(missing_counts, collapse=", "), ".") else NULL
  )
  guidance <- if (length(missing_meta) && !length(missing_counts)) {
    "The metadata appear to describe a subset of the count matrix. Subset the count matrix to the intended samples, or add metadata rows for the additional samples, before rerunning the analysis."
  } else if (length(missing_counts) && !length(missing_meta)) {
    "Remove metadata rows without count data, or add the corresponding count-matrix columns, before rerunning the analysis."
  } else {
    "Ensure that the two files contain exactly the same sample identifiers, then rerun the analysis."
  }
  stop(paste(
    "Count-matrix and metadata sample identifiers must match exactly.",
    paste(details, collapse=" "),
    guidance
  ))
}
counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"
shared <- colnames(counts_mat)
counts_mat <- counts_mat[, shared, drop=FALSE]
meta <- meta[shared, , drop=FALSE]

if (!is.null(filter_column) && !is.null(filter_value)) {
  if (!filter_column %in% colnames(meta))
    stop(sprintf("Filter column '%s' not found. Available: %s",
                 filter_column, paste(colnames(meta), collapse=", ")))
  keep_samples <- !is.na(meta[[filter_column]]) & as.character(meta[[filter_column]]) == filter_value
  if (!any(keep_samples)) stop(sprintf("No samples left after filter '%s' == '%s'", filter_column, filter_value))
  meta <- meta[keep_samples, , drop=FALSE]
  counts_mat <- counts_mat[, rownames(meta), drop=FALSE]
  cat(sprintf("Sample filter: %s == '%s'  (%d samples)\n\n", filter_column, filter_value, nrow(meta)))
}

typed <- apply_variable_types(meta, variable_types_file)
meta <- typed$data

# Model factors must retain at least two observed levels.
model_vars_precheck <- unique(all.vars(as.formula(formula_str)))
for (v in intersect(model_vars_precheck, colnames(meta))) {
  if (is.factor(meta[[v]]) || is.character(meta[[v]])) {
    observed_levels <- unique(as.character(meta[[v]][!is.na(meta[[v]])]))
    if (length(observed_levels) < 2L) {
      stop(sprintf("Categorical model variable '%s' has fewer than two observed levels after filtering.", v))
    }
  }
}

# Reject missing model metadata to preserve sample alignment.
model_vars <- unique(all.vars(as.formula(formula_str)))
missing_model_vars <- setdiff(model_vars, colnames(meta))
if (length(missing_model_vars))
  stop(sprintf("Model variable(s) not found in metadata: %s", paste(missing_model_vars, collapse=", ")))
na_model_vars <- model_vars[vapply(model_vars, function(v) anyNA(meta[[v]]), logical(1))]
if (length(na_model_vars)) {
  details <- vapply(na_model_vars, function(v) {
    ids <- rownames(meta)[is.na(meta[[v]])]
    sprintf("%s [%s]", v, paste(head(ids, 10L), collapse=", "))
  }, character(1))
  stop(sprintf(
    "Missing values detected in model variable(s): %s. MOSAIC-GSEA does not silently drop samples; complete the metadata or filter the samples explicitly.",
    paste(details, collapse="; ")
  ))
}

validate_mosaic_identifier <- function(x, what) {
  vals <- unique(as.character(x[!is.na(x)]))
  reserved_tokens <- c("_vs_", "__", "_to_", "×")
  bad_reserved <- vals[vapply(vals, function(v) any(vapply(reserved_tokens, function(tok) grepl(tok, v, fixed=TRUE), logical(1))), logical(1))]
  bad_path <- vals[grepl('[<>:"/\\\\|?*]|[[:cntrl:]]', vals, perl=TRUE)]
  bad_space <- vals[!nzchar(trimws(vals)) | trimws(vals) != vals]
  bad <- unique(c(bad_reserved, bad_path, bad_space))
  if (length(bad)) {
    stop(sprintf("Unsupported %s value(s): %s. Factor levels used by MOSAIC-GSEA must be non-empty, have no leading/trailing whitespace, and must not contain filename characters or reserved tokens (_vs_, __, _to_, ×).",
                 what, paste(head(bad, 10L), collapse=", ")))
  }
}

factor_cols <- names(meta)[vapply(meta, is.factor, logical(1))]
for (nm in factor_cols) validate_mosaic_identifier(meta[[nm]], sprintf("factor level in '%s'", nm))
validate_mosaic_identifier(colnames(meta), "metadata column name")

write.table(typed$types, file.path(output_dir, "metadata_variable_types.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
if (any(typed$types$source == "file")) {
  overridden <- typed$types[typed$types$source == "file", , drop=FALSE]
  cat(sprintf("Variable types overridden: %s\n",
              paste(sprintf("%s=%s", overridden$variable, overridden$type), collapse=", ")))
}

if (!contrast_factor %in% colnames(meta))
  stop(sprintf("contrast_factor '%s' not found. Available: %s",
               contrast_factor, paste(colnames(meta), collapse=", ")))
if (is.numeric(meta[[contrast_factor]]))
  stop(paste0("contrast_factor '", contrast_factor, "' is numeric. ",
              "The primary contrast factor must be categorical. ",
              "Continuous variables can be included as numeric covariates in the model formula."))

original_levels <- levels(meta[[contrast_factor]])
if (is.null(reference_level)) reference_level <- original_levels[1]
if (!reference_level %in% original_levels)
  stop(sprintf("reference_level '%s' not in '%s'. Available: %s",
               reference_level, contrast_factor, paste(original_levels, collapse=", ")))

meta[[contrast_factor]] <- relevel(meta[[contrast_factor]], ref=reference_level)

factor_levels  <- levels(meta[[contrast_factor]])
non_ref_levels <- setdiff(factor_levels, reference_level)

if (is.null(interaction_factors)) {
  interaction_factors <- auto_detect_interaction_factors(formula_str, contrast_factor)
  if (length(interaction_factors) > 0) invisible(interaction_factors)
} else if (length(interaction_factors) == 1 && !nzchar(interaction_factors)) {
  interaction_factors <- character(0)
}

int_levels     <- list()
int_refs       <- list()
int_nonref     <- list()

if (length(interaction_factors) > 0) {
  for (i in seq_along(interaction_factors)) {
    ifac <- interaction_factors[i]
    if (!ifac %in% colnames(meta))
      stop(sprintf("interaction_factor '%s' not found in metadata.", ifac))
    
    if (is.numeric(meta[[ifac]]))
      stop(paste0("interaction_factor '", ifac, "' is numeric. ",
                  "Automatic condition-specific contrasts currently require categorical interaction factors. ",
                  "Use continuous variables as additive covariates, or set the column to factor with --variable-types."))
    lvls <- levels(meta[[ifac]])
    ref <- lvls[1]
    
    meta[[ifac]] <- relevel(meta[[ifac]], ref=ref)
    int_levels[[ifac]] <- levels(meta[[ifac]])
    int_refs[[ifac]]   <- ref
    int_nonref[[ifac]] <- setdiff(lvls, ref)
  }
}

n_intfac <- length(interaction_factors)
if (n_intfac > 2L)
  stop("Interaction models with more than three factors are not supported. Use at most two interaction factors in addition to the contrast factor.")
if (n_intfac == 2L && !formula_has_interaction(formula_str, c(contrast_factor, interaction_factors)))
  stop(sprintf("Two interaction factors were detected for '%s', but the model formula does not contain the corresponding three-factor interaction term %s:%s:%s.",
               contrast_factor, contrast_factor, interaction_factors[1], interaction_factors[2]))
mode_str <- switch(as.character(n_intfac),
                   "0" = "STANDARD",
                   "1" = "2-FACTOR INTERACTION",
                   "2" = "3-FACTOR INTERACTION")

cat(sprintf("Analysis: %s | contrast: %s (reference: %s)\n", mode_str, contrast_factor, reference_level))
if (n_intfac > 0) {
  for (ifac in interaction_factors)
    cat(sprintf("  %s — reference: %s  |  non-ref: %s\n",
                ifac, int_refs[[ifac]], paste(int_nonref[[ifac]], collapse=", ")))
}
cat("\n")

if (n_intfac == 2L) {
  model_factors <- c(contrast_factor, interaction_factors)
  nlev <- vapply(model_factors, function(f) nlevels(factor(meta[[f]])), integer(1))
  if (any(nlev != 2L)) {
    stop(sprintf(
      "Three-factor interaction models are restricted to complete 2 x 2 x 2 designs. Levels found: %s",
      paste(sprintf("%s=%d", model_factors, nlev), collapse=", ")
    ))
  }
  observed <- do.call(table, lapply(meta[model_factors], factor))
  if (length(observed) != 8L || any(observed == 0L)) {
    stop("Three-factor interaction models require all eight factor combinations to be represented by at least one sample.")
  }
}

fixed_f       <- trimws(gsub("\\+?\\s*\\([^)]*\\|[^)]*\\)", "", formula_str))
fixed_formula <- as.formula(fixed_f)
if (attr(terms(fixed_formula), "intercept") != 1L) {
  stop(paste0(
    "No-intercept formulas (~ 0 + ... or ~ -1 + ...) are not supported. ",
    "MOSAIC-GSEA requires an intercept because its automatic reference-based contrasts ",
    "and variance-report interpretation assume an intercept-containing model."
  ))
}
fixed_design  <- model.matrix(fixed_formula, data=meta)
fixed_rank <- qr(fixed_design)$rank
if (fixed_rank < ncol(fixed_design)) {
  aliased <- tryCatch(limma::nonEstimable(fixed_design), error=function(e) character(0))
  stop(paste0("The model design matrix is not full rank; requested effects are not uniquely estimable.",
              if (length(aliased)) paste0(" Non-estimable coefficients: ", paste(aliased, collapse=", ")) else ""))
}
residual_df_design <- nrow(fixed_design) - fixed_rank
if (residual_df_design <= 0L) {
  stop(sprintf(
    "The fixed-effect design is saturated: %d samples, rank %d, residual degrees of freedom = %d. Add biological replication or simplify the model before running MOSAIC-GSEA.",
    nrow(fixed_design), fixed_rank, residual_df_design
  ))
}
if (residual_df_design < 2L) {
  warning(sprintf(
    "The fitted design has only %d residual degree(s) of freedom; inference may be unstable. More biological replication is strongly recommended.",
    residual_df_design
  ), call.=FALSE)
}

lib_sizes_raw <- colSums(counts_mat)
if (any(!is.finite(lib_sizes_raw) | lib_sizes_raw <= 0)) {
  bad_libs <- names(lib_sizes_raw)[!is.finite(lib_sizes_raw) | lib_sizes_raw <= 0]
  stop(sprintf("Count matrix contains zero/non-finite library size for sample(s): %s", paste(bad_libs, collapse=", ")))
}

dge           <- DGEList(counts=counts_mat)
keep          <- filterByExpr(dge, design=fixed_design)
if (!any(keep)) stop("filterByExpr removed all genes; no differential-expression analysis can be fitted.")
dge           <- dge[keep, , keep.lib.sizes=FALSE]
dge           <- calcNormFactors(dge, method="TMM")
cat(sprintf("filterByExpr: %d/%d genes retained\n", sum(keep), nrow(counts_mat)))

vp_bpparam <- if (.Platform$OS.type=="windows") SnowParam(cores,"SOCK") else MulticoreParam(cores)
if (has_random) {
  register(vp_bpparam)
  form <- as.formula(formula_str)
  vobj <- voomWithDreamWeights(dge, form, meta)
  fit <- NULL
  all_coefs <- colnames(fixed_design)
} else {
  cat("Fitting limma-voom...\n")
  design <- fixed_design
  vobj   <- voom(dge, design)
  fit    <- lmFit(vobj, design); fit <- eBayes(fit, robust=TRUE)
  all_coefs <- colnames(coef(fit))
}

vp_output_files <- file.path(output_dir, c("variance_partition_gene_wise.tsv",
                                           "variance_partition_summary.tsv",
                                           "variance_partition_status.tsv"))
unlink(vp_output_files[file.exists(vp_output_files)])
cat("Estimating gene-wise fractions of variance explained...\n")
vp_formula <- if (has_random) form else fixed_formula
vp_fixed_vars <- all.vars(fixed_formula)
vp_categorical_fixed <- vp_fixed_vars[vapply(vp_fixed_vars, function(v) {
  v %in% colnames(meta) && (is.factor(meta[[v]]) || is.character(meta[[v]]))
}, logical(1))]

write_vp_status <- function(status, method, note="", n_genes=0L,
                            max_dev=NA_real_) {
  vp_status <- data.frame(
    key=c("status", "method", "formula", "weights", "REML", "n_genes",
          "max_abs_component_sum_minus_one", "contrast_coding", "note"),
    value=c(status, method, formula_str, "voom precision weights", "FALSE",
            as.character(n_genes),
            if (is.finite(max_dev)) format(max_dev, scientific=TRUE) else "NA",
            paste(getOption("contrasts"), collapse=","), note),
    stringsAsFactors=FALSE
  )
  write.table(vp_status, file.path(output_dir, "variance_partition_status.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

write_vp_tables <- function(vp_df, vp_gene, method, note) {
  if (is.null(vp_df) || !nrow(vp_df)) {
    msg <- "Variance partitioning returned no complete gene-wise estimates."
    warning(msg, call.=FALSE)
    write_vp_status("FAILED", method, msg)
    return(invisible(FALSE))
  }
  vp_numeric <- vapply(vp_df, is.numeric, logical(1))
  if (!all(vp_numeric)) {
    bad <- names(vp_df)[!vp_numeric]
    msg <- sprintf("Variance partitioning returned non-numeric component column(s): %s",
                   paste(bad, collapse=", "))
    warning(msg, call.=FALSE)
    write_vp_status("FAILED", method, msg)
    return(invisible(FALSE))
  }
  if (!"Residuals" %in% colnames(vp_df)) {
    msg <- "Variance partitioning did not return the expected Residuals component."
    warning(msg, call.=FALSE)
    write_vp_status("FAILED", method, msg)
    return(invisible(FALSE))
  }
  complete <- apply(vp_df, 1, function(z) all(is.finite(z)))
  if (any(!complete)) {
    warning(sprintf("Variance partitioning returned incomplete estimates for %d gene(s); those genes were omitted from the variance report.", sum(!complete)), call.=FALSE)
    vp_df <- vp_df[complete, , drop=FALSE]
    vp_gene <- vp_gene[complete]
  }
  if (!nrow(vp_df) || length(vp_gene) != nrow(vp_df)) {
    msg <- "Variance partitioning gene identifiers could not be aligned to the expression matrix."
    warning(msg, call.=FALSE)
    write_vp_status("FAILED", method, msg)
    return(invisible(FALSE))
  }
  vp_components <- colnames(vp_df)
  vp_nonmodel <- intersect(c("Residuals", "Measurement.error"), vp_components)
  vp_model_components <- setdiff(vp_components, vp_nonmodel)
  vp_df$Model_explained <- if (length(vp_model_components)) rowSums(vp_df[, vp_model_components, drop=FALSE]) else 0
  vp_out <- data.frame(gene=vp_gene, vp_df, check.names=FALSE, stringsAsFactors=FALSE)
  write.table(vp_out, file.path(output_dir, "variance_partition_gene_wise.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

  vp_summary_components <- c(vp_components, "Model_explained")
  vp_summary <- do.call(rbind, lapply(vp_summary_components, function(nm) {
    x <- as.numeric(vp_df[[nm]])
    data.frame(component=nm, mean=mean(x), median=median(x),
               q25=as.numeric(quantile(x, 0.25, names=FALSE)),
               q75=as.numeric(quantile(x, 0.75, names=FALSE)),
               min=min(x), max=max(x), stringsAsFactors=FALSE, check.names=FALSE)
  }))
  write.table(vp_summary, file.path(output_dir, "variance_partition_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

  vp_component_sums <- rowSums(vp_df[, vp_components, drop=FALSE])
  vp_max_dev <- max(abs(vp_component_sums - 1))
  if (!is.finite(vp_max_dev) || vp_max_dev > 1e-6) {
    warning(sprintf("Variance fractions do not sum to one within numerical tolerance (maximum absolute deviation %.3g).", vp_max_dev), call.=FALSE)
  }
  write_vp_status("SUCCESS", method, note, nrow(vp_df), vp_max_dev)
  cat(sprintf("Variance partition: %d genes, %d modeled component(s) + residual variance.\n\n", nrow(vp_df), length(vp_model_components)))
  invisible(TRUE)
}

if (has_random && length(vp_categorical_fixed)) {
  # For categorical fixed effects in mixed models, report aggregate Nakagawa R2 fractions.
  method <- "Nakagawa R2 via performance::r2_nakagawa"
  if (!requireNamespace("performance", quietly=TRUE)) {
    msg <- "Mixed-model aggregate variance report requires the R package 'performance'."
    warning(msg, call.=FALSE)
    write_vp_status("SKIPPED", method, msg)
  } else {
    vp_fits <- tryCatch(
      suppressWarnings(variancePartition::fitVarPartModel(
        vobj, vp_formula, meta, REML=FALSE, useWeights=TRUE,
        showWarnings=FALSE, BPPARAM=vp_bpparam
      )),
      error=function(e) {
        msg <- sprintf("Mixed-model variance fitting failed for formula %s: %s", formula_str, conditionMessage(e))
        warning(msg, call.=FALSE)
        write_vp_status("FAILED", method, msg)
        NULL
      }
    )
    if (!is.null(vp_fits)) {
      fit_list <- tryCatch(as.list(vp_fits), error=function(e) vp_fits)
      genes <- names(fit_list)
      if (is.null(genes) || length(genes) != length(fit_list) || any(!nzchar(genes))) genes <- rownames(vobj$E)[seq_along(fit_list)]
      r2_tol <- 1e-8
      failure_reason <- rep("", length(fit_list))
      vals <- lapply(seq_along(fit_list), function(ii) {
        fit_i <- fit_list[[ii]]
        invalid <- function(reason) {
          failure_reason[ii] <<- reason
          c(Fixed_effects=NA_real_, Random_effects=NA_real_, Residuals=NA_real_)
        }
        if (is.null(fit_i)) return(invalid("model fit unavailable"))
        rr <- tryCatch(suppressWarnings(performance::r2_nakagawa(fit_i, verbose=FALSE)), error=function(e) NULL)
        if (is.null(rr)) return(invalid("Nakagawa R2 could not be estimated"))
        rm <- suppressWarnings(as.numeric(rr$R2_marginal[1]))
        rc <- suppressWarnings(as.numeric(rr$R2_conditional[1]))
        if (!is.finite(rm) || !is.finite(rc)) return(invalid("non-finite marginal/conditional R2"))
        if (rm < -r2_tol || rm > 1 + r2_tol || rc < -r2_tol || rc > 1 + r2_tol)
          return(invalid(sprintf("R2 outside [0,1]: marginal=%.6g, conditional=%.6g", rm, rc)))
        if (rc + r2_tol < rm)
          return(invalid(sprintf("conditional R2 below marginal R2: marginal=%.6g, conditional=%.6g", rm, rc)))

        # Correct only floating-point excursions within tolerance.
        rm <- min(1, max(0, rm))
        rc <- min(1, max(0, rc))
        if (rc < rm && rm - rc <= r2_tol) rc <- rm
        c(Fixed_effects=rm, Random_effects=rc-rm, Residuals=1-rc)
      })
      vp_df <- as.data.frame(do.call(rbind, vals), check.names=FALSE)
      failed_idx <- which(nzchar(failure_reason))
      failed_file <- file.path(output_dir, "variance_partition_failed_genes.tsv")
      if (length(failed_idx)) {
        failed_tab <- data.frame(
          gene=genes[failed_idx], reason=failure_reason[failed_idx],
          stringsAsFactors=FALSE, check.names=FALSE
        )
        write.table(failed_tab, failed_file, sep="\t", row.names=FALSE, quote=FALSE)
        cat(sprintf("Variance partition: excluded %d gene(s) with invalid/unavailable Nakagawa R2; see %s.\n",
                    length(failed_idx), basename(failed_file)))
      } else if (file.exists(failed_file)) {
        unlink(failed_file)
      }
      note <- paste0(
        "The differential-expression model contains categorical fixed effects and random effects. ",
        "Term-wise variancePartition fractions are not reported because variancePartition warns that categorical fixed effects in an LMM can give misleading term-specific fractions. ",
        "Instead, the report uses Nakagawa marginal R2 as the aggregate fraction explained by all fixed effects, conditional-minus-marginal R2 as the aggregate random-effect contribution, and one-minus-conditional R2 as residual variation."
      )
      write_vp_tables(vp_df, genes, method, note)
    }
  }
} else if (!has_random) {
  # Hierarchical factor-level Shapley decomposition of the weighted full-model R2.
  # Each model term joins a coalition only when ALL variables in that term are
  # present. Thus interactions are allocated across their participating factors
  # without choosing an arbitrary reference level or an order of entry.
  method <- "Hierarchical Shapley (voom-weighted R2)"
  shapley_result <- tryCatch({
    players <- unique(all.vars(fixed_formula))
    n_players <- length(players)
    if (!n_players || n_players > 10L)
      stop(sprintf("Hierarchical Shapley variance partition requires 1-10 modeled fixed variables (found %d); exhaustive coalitions would be impractical.", n_players))

    term_labels <- attr(terms(fixed_formula), "term.labels")
    term_players <- strsplit(term_labels, ":", fixed=TRUE)
    if (any(!vapply(term_players, function(x) length(x) > 0L && all(x %in% players), logical(1))))
      stop("Hierarchical Shapley requires plain fixed-variable names and their interactions; transformed formula terms are not supported.")

    column_term <- attr(fixed_design, "assign")
    if (length(column_term) != ncol(fixed_design) || !any(column_term == 0L))
      stop("Hierarchical Shapley could not identify the fitted design's intercept and model terms.")
    num_coalitions <- bitwShiftL(1L, n_players)
    masks <- 0:(num_coalitions - 1L)
    memberships <- vapply(seq_len(n_players), function(j) bitwAnd(masks, bitwShiftL(1L, j - 1L)) != 0L,
                          logical(num_coalitions))
    if (n_players == 1L) memberships <- matrix(memberships, ncol=1L)
    term_masks <- vapply(term_players, function(vars) {
      mask <- 0L
      for (v in vars) mask <- bitwOr(mask, bitwShiftL(1L, match(v, players) - 1L))
      mask
    }, integer(1))
    # Subset the ORIGINAL model matrix rather than refitting with new references.
    # Fixed basis and voom weights are identical across all nested models.
    subset_columns <- lapply(masks, function(mask) {
      selected_terms <- which(vapply(term_masks, function(tm) bitwAnd(mask, tm) == tm, logical(1)))
      which(column_term == 0L | column_term %in% selected_terms)
    })
    subset_by_player <- lapply(seq_len(n_players), function(j) {
      bit <- bitwShiftL(1L, j - 1L)
      idx <- which(bitwAnd(masks, bit) == 0L)
      sizes <- rowSums(memberships[idx, , drop=FALSE])
      list(base=idx, added=bitwOr(masks[idx], bit) + 1L,
           weight=exp(lgamma(sizes + 1) + lgamma(n_players - sizes) - lgamma(n_players + 1)))
    })

    expression <- as.matrix(vobj$E)
    precision <- vobj$weights
    if (is.null(precision) || !identical(dim(precision), dim(expression)))
      stop("Hierarchical Shapley needs a gene-by-sample voom precision-weight matrix.")
    if (ncol(expression) != nrow(fixed_design))
      stop("Hierarchical Shapley expression samples are not aligned to the fixed-effect design.")
    if (!identical(colnames(expression), rownames(fixed_design)))
      stop("Hierarchical Shapley sample identifiers differ between voom expression and the design.")
    if (any(!is.finite(expression)) || any(!is.finite(precision)) || any(precision <= 0))
      stop("Hierarchical Shapley requires finite expression values and strictly positive voom weights.")

    # Each R2 uses the SAME weighted response and the SAME gene-specific voom
    # weights; the intercept-only model supplies the weighted total sum of squares.
    one_gene <- function(i) {
      y <- expression[i, ]
      sqrt_w <- sqrt(precision[i, ])
      y_w <- y * sqrt_w
      y_mean <- sum(y * sqrt_w^2) / sum(sqrt_w^2)
      tss <- sum(((y - y_mean) * sqrt_w)^2)
      if (!is.finite(tss) || tss <= 1e-12 * max(1, sum(y_w^2)))
        return(rep(NA_real_, n_players + 1L))
      x_w <- fixed_design * sqrt_w
      r2 <- vapply(subset_columns, function(cols) {
        subfit <- qr(x_w[, cols, drop=FALSE])
        residual <- qr.resid(subfit, y_w)
        1 - sum(residual^2) / tss
      }, numeric(1))
      if (any(!is.finite(r2)) || any(r2 < -1e-8) || any(r2 > 1 + 1e-8))
        return(rep(NA_real_, n_players + 1L))
      r2 <- pmin(1, pmax(0, r2))
      allocations <- vapply(subset_by_player, function(sub) {
        gains <- r2[sub$added] - r2[sub$base]
        if (any(gains < -1e-8)) return(NA_real_)
        sum(sub$weight * pmax(0, gains))
      }, numeric(1))
      if (any(!is.finite(allocations)) ||
          abs(sum(allocations) - r2[num_coalitions]) > 1e-6)
        return(rep(NA_real_, n_players + 1L))
      c(allocations, 1 - r2[num_coalitions])
    }
    fractions <- t(vapply(seq_len(nrow(expression)), one_gene,
                          numeric(n_players + 1L)))
    colnames(fractions) <- c(players, "Residuals")
    rownames(fractions) <- rownames(expression)
    fractions
  }, error=function(e) {
    msg <- sprintf("Hierarchical Shapley variance partition failed for formula %s: %s", formula_str, conditionMessage(e))
    warning(msg, call.=FALSE)
    write_vp_status("FAILED", method, msg)
    NULL
  })
  if (!is.null(shapley_result)) {
    note <- paste0(
      "Gene-wise hierarchical Shapley allocation of voom-weighted full-model R2. ",
      "Each fixed variable is a player (multi-level categorical factors are grouped); a coalition includes a model term only if all its variables are present. ",
      "Interaction contributions are allocated to participating variables, rather than displayed as independent bars. ",
      "Nested models use the same samples, original design-matrix coding and per-gene voom weights. ",
      "Model_explained equals weighted full-model R2; Residuals equals one minus that R2. ",
      "Fractions measure model-dependent explanatory association, not causal effects."
    )
    write_vp_tables(as.data.frame(shapley_result, check.names=FALSE),
                    rownames(shapley_result), method, note)
  }
} else {
  # Unchanged mixed-model route for designs without categorical fixed effects.
  method <- "variancePartition::fitExtractVarPartModel"
  vp_result <- tryCatch(
    variancePartition::fitExtractVarPartModel(
      vobj, vp_formula, meta, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,
      BPPARAM=vp_bpparam
    ),
    error=function(e) {
      msg <- sprintf("Variance partitioning failed for formula %s: %s", formula_str, conditionMessage(e))
      warning(msg, call.=FALSE)
      write_vp_status("FAILED", method, msg)
      NULL
    }
  )
  if (!is.null(vp_result)) {
    vp_df <- as.data.frame(vp_result, check.names=FALSE)
    genes <- rownames(vp_df)
    if (is.null(genes) || any(!nzchar(genes))) genes <- rownames(vobj$E)[seq_len(nrow(vp_df))]
    note <- "Mixed-model variance fractions are reported term by term because the mixed model contains no categorical fixed effects; random effects and continuous fixed effects are estimated jointly."
    write_vp_tables(vp_df, genes, method, note)
  }
}

is_estimable_contrast <- function(L_named, design, tolerance=1e-8) {
  q <- qr(t(design), tol=tolerance)
  projected <- qr.fitted(q, as.numeric(L_named))
  max(abs(as.numeric(L_named) - projected)) <= tolerance * max(1, max(abs(L_named)))
}

dream_contrasts <- list()
contrast_vectors <- list()  # Exact fitted-model coefficient weights for audit.

format_contrast_result <- function(res, label) {
  res$gene <- rownames(res)
  res$original_gene <- unname(original_gene_by_internal[res$gene])
  missing_original <- is.na(res$original_gene) | !nzchar(res$original_gene)
  res$original_gene[missing_original] <- res$gene[missing_original]
  res$contrast <- label
  res <- res[!is.na(res$adj.P.Val), ]
  res <- res[match(rownames(dge$counts), res$gene), ]
  res <- res[!is.na(res$gene), ]; rownames(res) <- NULL
  keep_cols <- c("gene","original_gene","AveExpr","logFC","t")
  if ("z.std" %in% colnames(res)) keep_cols <- c(keep_cols, "z.std")
  keep_cols <- c(keep_cols, "P.Value","adj.P.Val","contrast")
  res[, keep_cols, drop=FALSE]
}

extract_contrast <- function(L_named, label) {
  if (length(L_named) != ncol(fixed_design) || !identical(names(L_named), colnames(fixed_design)))
    stop(sprintf("Contrast '%s' is not aligned to the fitted design matrix.", label))
  if (!is_estimable_contrast(L_named, fixed_design))
    stop(sprintf("Contrast '%s' is not estimable from the observed design.", label))
  if (sum(abs(L_named)) < 1e-12)
    stop(sprintf("Contrast '%s' is identically zero.", label))
  contrast_vectors[[label]] <<- L_named
  if (has_random) {
    dream_contrasts[[label]] <<- L_named
    return(NULL)
  }
  L   <- matrix(L_named, ncol=1, dimnames=list(names(L_named), label))
  tmp <- contrasts.fit(fit, L); tmp <- eBayes(tmp, robust=TRUE)
  res <- topTable(tmp, number=Inf, sort.by="none")
  format_contrast_result(res, label)
}

design_columns   <- colnames(fixed_design)
design_contrasts <- attr(fixed_design, "contrasts")
design_terms     <- terms(fixed_formula, data=meta)
design_xlevels   <- .getXlevels(design_terms, model.frame(fixed_formula, data=meta))

baseline_values <- lapply(meta, function(x) {
  if (is.factor(x)) levels(x)[1]
  else if (is.numeric(x)) mean(x, na.rm=TRUE)
  else x[which(!is.na(x))[1]]
})

build_design_row <- function(overrides=list()) {
  unknown <- setdiff(names(overrides), colnames(meta))
  if (length(unknown) > 0L)
    stop("Unknown metadata variable(s) in contrast condition: ",
         paste(unknown, collapse=", "))

  condition <- as.data.frame(baseline_values, stringsAsFactors=FALSE,
                             check.names=FALSE)
  for (nm in names(overrides)) condition[[nm]] <- overrides[[nm]]

  for (nm in colnames(meta)) {
    if (is.factor(meta[[nm]])) {
      raw_value <- as.character(condition[[nm]][1])
      if (!raw_value %in% levels(meta[[nm]]))
        stop(sprintf("Unknown level '%s' for factor '%s'. Available: %s",
                     raw_value, nm, paste(levels(meta[[nm]]), collapse=", ")))
      condition[[nm]] <- factor(raw_value, levels=levels(meta[[nm]]))
    } else if (is.numeric(meta[[nm]])) {
      condition[[nm]] <- as.numeric(condition[[nm]])
    }
  }

  row <- model.matrix(fixed_formula, data=condition,
                      contrasts.arg=design_contrasts, xlev=design_xlevels)

  extra <- setdiff(colnames(row), design_columns)
  if (length(extra) > 0L)
    stop("Unexpected design columns while constructing contrast: ",
         paste(extra, collapse=", "))

  missing <- setdiff(design_columns, colnames(row))
  if (length(missing) > 0L) {
    pad <- matrix(0, nrow=nrow(row), ncol=length(missing),
                  dimnames=list(NULL, missing))
    row <- cbind(row, pad)
  }
  row <- row[, design_columns, drop=FALSE]

  absent_in_fit <- setdiff(design_columns, all_coefs)
  absent_in_design <- setdiff(all_coefs, design_columns)
  if (length(absent_in_fit) > 0L || length(absent_in_design) > 0L) {
    stop(paste0(
      "Fixed-effect design columns and fitted coefficients do not match. ",
      if (length(absent_in_fit)) paste0("Only in design: ", paste(absent_in_fit, collapse=", "), ". ") else "",
      if (length(absent_in_design)) paste0("Only in fit: ", paste(absent_in_design, collapse=", "), ".") else ""
    ))
  }
  out <- as.numeric(row[1, all_coefs, drop=TRUE])
  names(out) <- all_coefs
  out
}

build_condition_L <- function(cf_level, target_combo=list()) {
  num <- target_combo
  den <- target_combo
  num[[contrast_factor]] <- cf_level
  den[[contrast_factor]] <- reference_level
  build_design_row(num) - build_design_row(den)
}

build_interaction_L <- function(cf_level, int_combo) {
  included <- c(contrast_factor, names(int_combo))
  target <- c(setNames(list(cf_level), contrast_factor), int_combo)
  refs <- c(setNames(list(reference_level), contrast_factor), int_refs[names(int_combo)])
  k <- length(included)
  bit_grid <- expand.grid(rep(list(c(0L, 1L)), k), KEEP.OUT.ATTRS=FALSE,
                          stringsAsFactors=FALSE)
  names(bit_grid) <- included
  L <- setNames(numeric(length(all_coefs)), all_coefs)

  for (rr in seq_len(nrow(bit_grid))) {
    bits <- as.integer(bit_grid[rr, , drop=TRUE])
    values <- lapply(seq_along(included), function(j) {
      nm <- included[j]
      if (bits[j] == 1L) target[[nm]] else refs[[nm]]
    })
    names(values) <- included

    for (other in setdiff(interaction_factors, names(int_combo)))
      values[[other]] <- int_refs[[other]]

    weight <- (-1)^(k - sum(bits))
    L <- L + weight * build_design_row(values)
  }
  L
}

save_result <- function(df, label) {
  if (is.null(df)) return(invisible(NULL))
  out_file <- file.path(output_dir, paste0(label, ".txt"))
  write.table(df, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  out_file
}

results  <- list()
manifest <- list()

if (n_intfac == 0) {
  cat("Contrasts:\n")
  for (lvl in non_ref_levels) {
    label <- paste0(lvl, "_vs_", reference_level)
    L <- build_condition_L(lvl, list())
    res <- extract_contrast(L, label); if (!is.null(res)) results[[label]] <- res
    save_result(res, label)
    manifest[[label]] <- data.frame(label=label, condition_label=label,
                                    role="vs_ref", header1=lvl,
                                    stringsAsFactors=FALSE)
  }
  if (length(non_ref_levels) > 1) {
    cat("Pairwise contrasts:\n")
    for (p in combn(non_ref_levels, 2, simplify=FALSE)) {
      lvl_a <- p[1]; lvl_b <- p[2]
      label  <- paste0(lvl_a, "_vs_", lvl_b)
      L <- build_condition_L(lvl_a, list()) - build_condition_L(lvl_b, list())
      res <- extract_contrast(L, label); if (!is.null(res)) results[[label]] <- res
      save_result(res, label)
      manifest[[label]] <- data.frame(label=label, condition_label=label,
                                      role="pairwise", header1=lvl_a,
                                      stringsAsFactors=FALSE)
    }
  }
}

if (n_intfac >= 1) {
  
  f_clean_order <- gsub("\\s+", "", gsub("\\([^)]*\\|[^)]*\\)", "", formula_str))
  f_clean_order <- gsub("~","", f_clean_order)
  positions <- vapply(c(contrast_factor, interaction_factors), function(fac) {
    m <- regexpr(paste0("(?<![A-Za-z0-9_.])", fac, "(?![A-Za-z0-9_.])"), f_clean_order, perl=TRUE)
    if (m[1] < 0) Inf else as.numeric(m[1])
  }, numeric(1))
  
  int_order <- interaction_factors[order(positions[interaction_factors])]
  
  header_factors <- int_order
  
  combo_grid <- expand.grid(lapply(interaction_factors, function(ifac) int_levels[[ifac]]),
                            stringsAsFactors=FALSE)
  colnames(combo_grid) <- interaction_factors
  
  for (cf_lvl in non_ref_levels) {
    cat(sprintf("Contrasts for %s:\n", cf_lvl))
    
    for (r in seq_len(nrow(combo_grid))) {
      combo <- as.list(combo_grid[r, , drop=FALSE])
      
      lhs_parts <- c(cf_lvl,  unname(unlist(combo)))
      rhs_parts <- c(reference_level, unname(unlist(combo)))
      label <- paste0(paste(lhs_parts, collapse="×"), "_vs_",
                      paste(rhs_parts, collapse="×"))
      
      L   <- build_condition_L(cf_lvl, combo)
      if (sum(abs(L)) < 1e-12) { cat(sprintf("  WARNING: empty L for %s\n", label)); next }
      res <- extract_contrast(L, label); if (!is.null(res)) results[[label]] <- res
      save_result(res, label)
      
      hdr <- vapply(header_factors, function(hf) unlist(combo)[[hf]], character(1))
      manifest[[label]] <- data.frame(
        label          = label,
        condition_label = paste(c(cf_lvl, unname(unlist(combo))), collapse="×"),
        role           = "simple_effect",
        contrast_level = cf_lvl,
        stringsAsFactors = FALSE
      )
      for (hf in header_factors) manifest[[label]][[paste0("header_", hf)]] <- unlist(combo)[[hf]]
    }
    
    if (identical(cf_lvl, non_ref_levels[1]) && length(non_ref_levels) > 1L) {
      cat("Pairwise simple effects within strata:\n")
      for (r_pw in seq_len(nrow(combo_grid))) {
        combo_pw <- as.list(combo_grid[r_pw, , drop=FALSE])
        combo_vals <- unname(unlist(combo_pw))
        for (pw in combn(non_ref_levels, 2, simplify=FALSE)) {
          lvl_a <- pw[1]; lvl_b <- pw[2]
          lhs_parts <- c(lvl_a, combo_vals)
          rhs_parts <- c(lvl_b, combo_vals)
          label_pw <- paste0(paste(lhs_parts, collapse="×"), "_vs_",
                             paste(rhs_parts, collapse="×"))
          L_a <- build_condition_L(lvl_a, combo_pw)
          L_b <- build_condition_L(lvl_b, combo_pw)
          L_pw <- L_a - L_b
          if (sum(abs(L_pw)) == 0) {
            cat(sprintf("  WARNING: empty pairwise L for %s\n", label_pw)); next
          }
          res <- extract_contrast(L_pw, label_pw); if (!is.null(res)) results[[label_pw]] <- res
          save_result(res, label_pw)
          manifest[[label_pw]] <- data.frame(
            label=label_pw,
            condition_label=paste(c(lvl_a, combo_vals), collapse="×"),
            role="pairwise_simple",
            contrast_level=lvl_a,
            stringsAsFactors=FALSE
          )
        }
      }
    }

    if (identical(cf_lvl, non_ref_levels[1]) && length(non_ref_levels) > 1L) {
      cat("Marginal pairwise contrasts across strata:\n")
      mean_condition_L <- function(cf_level) {
        Ls <- lapply(seq_len(nrow(combo_grid)), function(rr) {
          build_condition_L(cf_level, as.list(combo_grid[rr, , drop=FALSE]))
        })
        Reduce(`+`, Ls) / length(Ls)
      }
      marginal_L <- setNames(lapply(non_ref_levels, mean_condition_L), non_ref_levels)
      for (pw in combn(non_ref_levels, 2, simplify=FALSE)) {
        lvl_a <- pw[1]; lvl_b <- pw[2]
        label_pw <- paste0("marginal__", lvl_a, "_vs_", lvl_b)
        L_pw <- marginal_L[[lvl_a]] - marginal_L[[lvl_b]]
        if (sum(abs(L_pw)) == 0) {
          cat(sprintf("  WARNING: empty marginal pairwise L for %s\n", label_pw)); next
        }
        res <- extract_contrast(L_pw, label_pw); if (!is.null(res)) results[[label_pw]] <- res
        save_result(res, label_pw)
        manifest[[label_pw]] <- data.frame(
          label=label_pw, condition_label=lvl_a,
          role="pairwise_marginal", contrast_level=lvl_a,
          stringsAsFactors=FALSE
        )
      }
    }

    for (ifac in interaction_factors) {
      ref_lv <- int_refs[[ifac]]
      nonref_lv <- int_nonref[[ifac]]
      level_pairs <- lapply(nonref_lv, function(nr_lv) c(nr_lv, ref_lv))
      if (length(nonref_lv) > 1L)
        level_pairs <- c(level_pairs, combn(nonref_lv, 2, simplify=FALSE))

      other_factors <- setdiff(interaction_factors, ifac)
      fixed_grid <- if (length(other_factors)) {
        z <- expand.grid(lapply(other_factors, function(f) int_levels[[f]]),
                         stringsAsFactors=FALSE)
        colnames(z) <- other_factors
        z
      } else data.frame(.dummy=1L)

      for (pair in level_pairs) {
        lv_a <- pair[1]; lv_b <- pair[2]
        for (rr in seq_len(nrow(fixed_grid))) {
          fixed_combo <- if (length(other_factors))
            as.list(fixed_grid[rr, other_factors, drop=FALSE]) else list()
          combo_a <- c(setNames(list(lv_a), ifac), fixed_combo)
          combo_b <- c(setNames(list(lv_b), ifac), fixed_combo)
          cond_tag <- if (length(other_factors)) paste0("__at__",
            paste(vapply(other_factors, function(f)
              paste0(f, "=", fixed_combo[[f]]), character(1)), collapse="__")) else ""
          label_int <- paste0(cf_lvl, "×", ifac, "__", lv_a, "_vs_", lv_b, cond_tag)
          L <- build_condition_L(cf_lvl, combo_a) - build_condition_L(cf_lvl, combo_b)
          if (sum(abs(L)) < 1e-12) { cat(sprintf("  WARNING: empty L for %s\n", label_int)); next }
          res <- extract_contrast(L, label_int); if (!is.null(res)) results[[label_int]] <- res
          save_result(res, label_int)
          manifest[[label_int]] <- data.frame(
            label=label_int, condition_label=paste0("x_",ifac,"_",lv_a,"_vs_",lv_b),
            role=paste0("interaction_x_", ifac),
            contrast_level=cf_lvl, interaction_factor=ifac,
            tested_level=lv_a, tested_reference=lv_b,
            conditioned_factor=if (length(other_factors)) paste(other_factors, collapse=";") else NA_character_,
            conditioned_level=if (length(other_factors)) paste(unname(unlist(fixed_combo)), collapse=";") else NA_character_,
            stringsAsFactors=FALSE)
        }
      }
    }
    
    if (n_intfac >= 2) {
      nr_combos <- expand.grid(lapply(interaction_factors, function(ifac) int_nonref[[ifac]]),
                               stringsAsFactors=FALSE)
      colnames(nr_combos) <- interaction_factors
      for (r in seq_len(nrow(nr_combos))) {
        int_combo  <- as.list(nr_combos[r, , drop=FALSE])
        level_tag <- paste(vapply(names(int_combo), function(f)
          paste0(f, "=", int_combo[[f]], "_vs_", int_refs[[f]]), character(1)), collapse="__")
        label_int  <- paste0(cf_lvl, "×", paste(names(int_combo), collapse="×"), "__", level_tag)
        L   <- build_interaction_L(cf_lvl, int_combo)
        if (sum(abs(L)) < 1e-12) { cat(sprintf("  WARNING: empty L for %s\n", label_int)); next }
        res <- extract_contrast(L, label_int); if (!is.null(res)) results[[label_int]] <- res
        save_result(res, label_int)
        manifest[[label_int]] <- data.frame(
          label=label_int, condition_label=paste0("x_", paste(unname(unlist(int_combo)), collapse="_x_")),
          role=paste0("interaction_", n_intfac+1, "way"),
          contrast_level=cf_lvl,
          interaction_factor=paste(names(int_combo), collapse=";"),
          tested_level=paste(unname(unlist(int_combo)), collapse=";"),
          tested_reference=paste(unname(unlist(int_refs[names(int_combo)])), collapse=";"),
          stringsAsFactors=FALSE)
      }
    }
    cat("\n")
  }
}

if (has_random) {
  if (!length(dream_contrasts)) stop("No estimable contrasts were generated for the mixed model.")
  L_all <- do.call(cbind, dream_contrasts)
  colnames(L_all) <- names(dream_contrasts)
  cat(sprintf("Fitting dream (mixed model) with %d contrast(s)...\n", ncol(L_all)))
  fit <- dream(vobj, form, meta, L=L_all)
  fit <- eBayes(fit, robust=TRUE)
  for (label in colnames(L_all)) {
    res <- topTable(fit, coef=label, number=Inf, sort.by="none")
    if (!"z.std" %in% colnames(res))
      stop(sprintf("dream did not return z.std for contrast '%s'.", label))
    res <- format_contrast_result(res, label)
    results[[label]] <- res
    save_result(res, label)
  }
}

analysis_info <- data.frame(
  key   = c("analysis_mode","contrast_factor","reference_level",
            "interaction_factors","n_interaction_factors",
            "formula_str","formula_expanded"),
  value = c(
    if (n_intfac==0) "standard" else paste0("interaction_", n_intfac+1, "way"),
    contrast_factor, reference_level,
    paste(interaction_factors, collapse=";"),
    as.character(n_intfac),
    formula_str, expanded_formula_str
  ),
  stringsAsFactors = FALSE
)
if (n_intfac > 0) {
  for (ifac in interaction_factors) {
    analysis_info <- rbind(analysis_info, data.frame(
      key   = c(paste0("int_factor_", ifac),
                paste0("int_ref_", ifac),
                paste0("int_nonref_", ifac)),
      value = c(ifac, int_refs[[ifac]], paste(int_nonref[[ifac]], collapse=";")),
      stringsAsFactors=FALSE))
  }
  analysis_info <- rbind(analysis_info, data.frame(
    key=c("header_factor_order", "condition_component_order"),
    value=c(paste(c(contrast_factor, int_order), collapse=";"),
            paste(c(contrast_factor, interaction_factors), collapse=";")),
    stringsAsFactors=FALSE))
}
write.table(analysis_info, file.path(output_dir, "analysis_manifest.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

if (length(manifest) > 0) {
  all_cols <- unique(unlist(lapply(manifest, colnames)))
  cm <- do.call(rbind, lapply(manifest, function(r) {
    missing <- setdiff(all_cols, colnames(r))
    for (col in missing) r[[col]] <- NA_character_
    r[, all_cols, drop=FALSE]
  }))
  rownames(cm) <- NULL
  write.table(cm, file.path(output_dir, "contrast_manifest.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

# A long-format audit trail: each label and each design-matrix coefficient weight.
# Positive and negative weights give the exact direction of the fitted contrast.
if (length(contrast_vectors)) {
  coefficient_table <- do.call(rbind, lapply(names(contrast_vectors), function(lbl) {
    v <- contrast_vectors[[lbl]]
    data.frame(contrast = rep(lbl, length(v)),
               model_coefficient = names(v),
               weight = as.numeric(v),
               stringsAsFactors = FALSE, check.names = FALSE)
  }))
  write.table(coefficient_table,
              file.path(output_dir, "contrast_coefficients.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("Contrasts completed: %d\n", length(results)))
cat(sprintf("\n  %-50s  %6s  %5s  %5s\n", "Contrast", "sig", "up", "dn"))
for (label in names(results)) {
  df <- results[[label]]
  n_sig <- sum(df$adj.P.Val < sig_threshold, na.rm=TRUE)
  cat(sprintf("  %-50s  %6d  %5d  %5d\n", label, n_sig,
              sum(df$adj.P.Val < sig_threshold & df$logFC > 0, na.rm=TRUE),
              sum(df$adj.P.Val < sig_threshold & df$logFC < 0, na.rm=TRUE)))
}
