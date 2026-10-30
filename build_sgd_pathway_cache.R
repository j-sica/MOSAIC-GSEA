#!/usr/bin/env Rscript

# Build the S. cerevisiae pathway cache from SGD biochemical_pathways.tab.
# Usage: Rscript build_sgd_pathway_cache.R [output.rds] [table-or-URL]

args <- commandArgs(trailingOnly = TRUE)

out_file <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  file.path(
    "saccharomyces_cerevisiae_GO_and_KEGG_data",
    "sgd_pathway_cache_sce.rds"
  )
}

source_arg <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else NA_character_

primary_url <- "https://downloads.yeastgenome.org/curation/literature/biochemical_pathways.tab"
fallback_urls <- c(
  primary_url,
  "https://sgd-archive.yeastgenome.org/curation/literature/biochemical_pathways.tab"
)

options(timeout = max(180, getOption("timeout", 60)))

required <- c("AnnotationDbi", "org.Sc.sgd.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    paste0(
      "Missing required package(s): ", paste(missing, collapse = ", "), "\n",
      "Install them with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install(c('AnnotationDbi','org.Sc.sgd.db'))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(AnnotationDbi))
orgdb_obj <- getExportedValue("org.Sc.sgd.db", "org.Sc.sgd.db")

if (!("ENSEMBL" %in% AnnotationDbi::keytypes(orgdb_obj))) {
  stop("org.Sc.sgd.db does not expose the ENSEMBL key type required by MOSAIC-GSEA.")
}

all_orfs <- sort(unique(as.character(
  AnnotationDbi::keys(orgdb_obj, keytype = "ENSEMBL")
)))
all_orfs <- all_orfs[!is.na(all_orfs) & nzchar(all_orfs)]
if (!length(all_orfs)) stop("No systematic yeast identifiers were obtained from org.Sc.sgd.db.")

# Map SGD aliases to systematic ORFs.
org_cols <- intersect(
  c("SYMBOL", "GENENAME", "ORF", "SGD"),
  AnnotationDbi::columns(orgdb_obj)
)

alias_to_orf <- setNames(all_orfs, toupper(all_orfs))
if (length(org_cols)) {
  ann <- suppressMessages(AnnotationDbi::select(
    orgdb_obj,
    keys = all_orfs,
    keytype = "ENSEMBL",
    columns = org_cols
  ))
  if (nrow(ann)) {
    for (col in org_cols) {
      vals <- as.character(ann[[col]])
      keep <- !is.na(vals) & nzchar(trimws(vals)) & !is.na(ann$ENSEMBL)
      if (any(keep)) {
        keys <- toupper(trimws(vals[keep]))
        vals_orf <- as.character(ann$ENSEMBL[keep])
        # Keep unambiguous aliases.
        spl <- split(vals_orf, keys)
        spl <- spl[vapply(spl, function(x) length(unique(x)) == 1L, logical(1))]
        if (length(spl)) {
          add <- vapply(spl, function(x) unique(x)[1], character(1))
          alias_to_orf[names(add)] <- add
        }
      }
    }
  }
}

is_url <- function(x) grepl("^https?://", x, ignore.case = TRUE)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

fetch_bulk <- function() {
  if (!is.na(source_arg) && !is_url(source_arg)) {
    if (!file.exists(source_arg)) stop("Local biochemical_pathways.tab not found: ", source_arg)
    cat("Using local SGD bulk file:\n  ", normalizePath(source_arg), "\n", sep = "")
    return(list(path = normalizePath(source_arg), source = normalizePath(source_arg), downloaded = FALSE))
  }

  urls <- if (!is.na(source_arg) && is_url(source_arg)) source_arg else fallback_urls
  dest <- file.path(dirname(out_file), "biochemical_pathways.tab")
  errors <- character(0)

  for (u in urls) {
    cat("Downloading SGD biochemical pathway table:\n  ", u, "\n", sep = "")
    ok <- tryCatch({
      suppressWarnings(utils::download.file(u, destfile = dest, mode = "wb", quiet = FALSE, method = "libcurl"))
      file.exists(dest) && is.finite(file.info(dest)$size) && file.info(dest)$size > 1000
    }, error = function(e) {
      errors <<- c(errors, paste0(u, " -> ", conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) return(list(path = dest, source = u, downloaded = TRUE))
  }

  stop(
    paste0(
      "Could not download biochemical_pathways.tab from SGD.\n",
      if (length(errors)) paste(errors, collapse = "\n") else "No download endpoint succeeded.",
      "\n\nYou can download the file manually and rerun:\n",
      "  Rscript build_sgd_pathway_cache.R <output.rds> biochemical_pathways.tab"
    ),
    call. = FALSE
  )
}

bulk <- fetch_bulk()
raw_lines <- readLines(bulk$path, warn = FALSE, encoding = "UTF-8")
raw_lines <- sub("^\\ufeff", "", raw_lines)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
raw_lines <- raw_lines[!grepl("^\\s*#", raw_lines)]
if (length(raw_lines) < 20L) stop("Downloaded biochemical_pathways.tab is unexpectedly small or empty.")

rows <- strsplit(raw_lines, "\\t", fixed = FALSE)
ncols <- vapply(rows, length, integer(1))
modal_ncol <- as.integer(names(sort(table(ncols), decreasing = TRUE))[1])
rows <- rows[ncols == modal_ncol]
if (length(rows) < 20L || modal_ncol < 2L) {
  stop("Could not parse biochemical_pathways.tab as a regular tab-delimited table.")
}

mat <- do.call(rbind, rows)
mat <- apply(mat, 2, trimws)
if (is.null(dim(mat))) mat <- matrix(mat, ncol = modal_ncol)

# Remove a header row when present.
header_words <- tolower(paste(mat[1, ], collapse = " "))
if (grepl("pathway|gene|enzyme|reference|ec number|ec_number", header_words)) {
  mat <- mat[-1, , drop = FALSE]
}

# Identify columns that map reliably to yeast genes.
field_tokens <- function(x) {
  x <- gsub("[(){}\\[\\],;|/]", " ", x)
  x <- gsub("[:=]", " ", x)
  unique(toupper(unlist(strsplit(x, "[[:space:]]+", perl = TRUE))))
}

map_field <- function(x) {
  toks <- field_tokens(x)
  toks <- toks[nzchar(toks)]
  out <- unname(alias_to_orf[toks])
  unique(out[!is.na(out) & nzchar(out)])
}

mapped_counts <- vapply(seq_len(ncol(mat)), function(j) {
  sum(vapply(mat[, j], function(x) length(map_field(x)) > 0L, logical(1)))
}, integer(1))

# SGD schema: pathway, enzyme, EC number, gene, references.
reference_field <- function(x) {
  z <- gsub("[[:space:]]+", "", trimws(x))
  grepl(
    "^(?:(?:PMID|SGD_REF):[^|]+)(?:\\|(?:PMID|SGD_REF):[^|]+)*$",
    z, ignore.case = TRUE, perl = TRUE
  )
}

min_gene_rows <- max(10L, floor(0.05 * nrow(mat)))
canonical_schema <- ncol(mat) >= 5L &&
  mapped_counts[4] >= min_gene_rows &&
  mean(reference_field(mat[, 1]), na.rm = TRUE) < 0.10

if (isTRUE(canonical_schema)) {
  path_col <- 1L
  gene_cols <- 4L
  cat("Detected canonical SGD biochemical_pathways.tab schema: pathway=column 1, gene=column 4.\n")
} else {
  # Fallback gene-column detection.
  gene_cols <- which(mapped_counts >= min_gene_rows)
  if (!length(gene_cols)) {
    stop(
      paste0(
        "Could not identify a gene/ORF column in biochemical_pathways.tab.\n",
        "Rows: ", nrow(mat), "; columns: ", ncol(mat), ".\n",
        "The SGD file format may have changed. No cache was written."
      ),
      call. = FALSE
    )
  }

  # Infer the pathway-name column when the canonical schema is unavailable.
  pathway_score <- rep(-Inf, ncol(mat))
  for (j in seq_len(ncol(mat))) {
    if (j %in% gene_cols) next
    vals <- trimws(mat[, j])
    vals <- vals[nzchar(vals)]
    if (!length(vals)) next
    u <- unique(vals)
    nu <- length(u)
    rep_frac <- 1 - nu / length(vals)
    alpha_frac <- mean(grepl("[A-Za-z]", vals))
    ec_frac <- mean(grepl("^[0-9-]+(\\.[0-9-]+){2,3}$", vals))
    ref_frac <- mean(reference_field(vals))
    ref_token_frac <- mean(grepl("(?:PMID|SGD_REF):", vals, ignore.case = TRUE, perl = TRUE))
    plausible_n <- nu >= 10L && nu <= min(1000L, floor(0.9 * nrow(mat)))
    if (plausible_n && alpha_frac > 0.5 && ec_frac < 0.5 && ref_frac < 0.10 && ref_token_frac < 0.10) {
      # Prefer descriptive text and earlier columns.
      pathway_score[j] <- rep_frac + 0.20 * alpha_frac - 0.002 * j
    }
  }

  path_col <- which.max(pathway_score)
  if (!length(path_col) || !is.finite(pathway_score[path_col])) {
    stop(
      paste0(
        "Could not identify the pathway-name column in biochemical_pathways.tab.\n",
        "The SGD file format may have changed. No cache was written."
      ),
      call. = FALSE
    )
  }
}

path_names <- trimws(mat[, path_col])
valid_rows <- nzchar(path_names) & !reference_field(path_names)
mat <- mat[valid_rows, , drop = FALSE]
path_names <- path_names[valid_rows]

# Use Pathway Tools IDs when present.
id_pattern <- "^(?:[A-Z0-9]+-)?PWY-[A-Z0-9-]+$|^PWY[0-9A-Z-]+$"
id_col <- NA_integer_
for (j in setdiff(seq_len(ncol(mat)), c(gene_cols, path_col))) {
  vals <- trimws(mat[, j])
  hit <- mean(grepl(id_pattern, vals, perl = TRUE))
  if (is.finite(hit) && hit >= 0.25) {
    id_col <- j
    break
  }
}

# Resolve pathway names to YeastPathways IDs when available.
html_decode_basic <- function(x) {
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x <- gsub("&quot;", '"', x, fixed = TRUE)
  x <- gsub("&#39;", "'", x, fixed = TRUE)
  x <- gsub("&apos;", "'", x, fixed = TRUE)
  x <- gsub("&nbsp;", " ", x, fixed = TRUE)
  x <- gsub("&gamma;", "gamma", x, fixed = TRUE)
  x <- gsub("&alpha;", "alpha", x, fixed = TRUE)
  x <- gsub("&beta;", "beta", x, fixed = TRUE)
  x
}

norm_path_name <- function(x) {
  x <- html_decode_basic(x)
  x <- gsub("<[^>]+>", "", x)
  x <- gsub("[[:space:]]+", " ", x)
  tolower(trimws(x))
}

resolve_yeastpathway_catalog <- function(pathway_names) {
  listing_url <- "https://pathway.yeastgenome.org/YEAST/hierarchical.html"
  out <- data.frame(
    name = pathway_names,
    object_id = NA_character_,
    url = rep(listing_url, length(pathway_names)),
    stringsAsFactors = FALSE
  )
  if (!length(pathway_names)) return(out)

  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({
    suppressWarnings(utils::download.file(listing_url, destfile = tmp, mode = "wb", quiet = TRUE, method = "libcurl"))
    file.exists(tmp) && is.finite(file.info(tmp)$size) && file.info(tmp)$size > 1000
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    cat("  YeastPathways direct-link catalog unavailable; using hierarchical-list links.\n")
    return(out)
  }

  html <- paste(readLines(tmp, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Extract pathway links regardless of HTML attribute order.
  anchors <- regmatches(
    html,
    gregexpr("<a\\b[^>]*href=[\\\"'][^\\\"']*NEW-IMAGE\\?[^\\\"']*[\\\"'][^>]*>.*?</a>",
             html, perl = TRUE, ignore.case = TRUE)
  )[[1]]
  if (!length(anchors) || identical(anchors, "")) return(out)

  hrefs <- sub(".*href=[\\\"']([^\\\"']+)[\\\"'].*", "\\1", anchors,
               perl = TRUE, ignore.case = TRUE)
  texts <- sub("^.*?>", "", anchors, perl = TRUE)
  texts <- sub("</a>.*$", "", texts, perl = TRUE, ignore.case = TRUE)
  texts <- html_decode_basic(gsub("<[^>]+>", "", texts))
  hrefs <- gsub("&amp;", "&", hrefs, fixed = TRUE)
  is_path <- grepl("(?:[?&])type=PATHWAY(?:&|$)", hrefs, ignore.case = TRUE, perl = TRUE)
  hrefs <- hrefs[is_path]
  texts <- texts[is_path]
  if (!length(hrefs)) return(out)

  object_ids <- sub("^.*(?:[?&])object=([^&]+).*$", "\\1", hrefs, perl = TRUE)
  valid_obj <- nzchar(object_ids) & object_ids != hrefs
  hrefs <- hrefs[valid_obj]
  texts <- texts[valid_obj]
  object_ids <- object_ids[valid_obj]
  if (!length(object_ids)) return(out)

  keys <- norm_path_name(texts)
  # Keep the first duplicate pathway label.
  keep <- !duplicated(keys) & nzchar(keys)
  keys <- keys[keep]
  object_ids <- object_ids[keep]
  hrefs <- hrefs[keep]

  idx <- match(norm_path_name(pathway_names), keys)
  hit <- !is.na(idx)
  if (any(hit)) {
    out$object_id[hit] <- object_ids[idx[hit]]
    out$url[hit] <- paste0(
      "https://pathway.yeastgenome.org/YEAST/NEW-IMAGE?type=PATHWAY&object=",
      utils::URLencode(out$object_id[hit], reserved = TRUE)
    )
  }
  cat(sprintf("  YeastPathways direct links resolved: %d/%d pathway names\n", sum(hit), length(pathway_names)))
  out
}

# Build pathway gene membership.
row_genes <- lapply(seq_len(nrow(mat)), function(i) {
  unique(unlist(lapply(gene_cols, function(j) map_field(mat[i, j])), use.names = FALSE))
})

keep <- vapply(row_genes, length, integer(1)) > 0L
if (!any(keep)) stop("No pathway-associated genes could be mapped to systematic ORFs.")
path_names2 <- path_names[keep]
row_genes2 <- row_genes[keep]
mat2 <- mat[keep, , drop = FALSE]

# Normalize pathway display names.
clean_name <- function(x) {
  x <- gsub("<[^>]+>", "", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}
path_names2 <- vapply(path_names2, clean_name, character(1))

# Prefer Pathway Tools IDs; otherwise use stable local IDs.
unique_names <- sort(unique(path_names2))
local_id_map <- setNames(sprintf("SGD_BIOCHEM_%03d", seq_along(unique_names)), unique_names)
path_catalog <- resolve_yeastpathway_catalog(unique_names)
catalog_id <- setNames(path_catalog$object_id, path_catalog$name)
catalog_url <- setNames(path_catalog$url, path_catalog$name)

pathways <- setNames(vector("list", length(all_orfs)), all_orfs)
status <- setNames(as.list(rep("no_annotation", length(all_orfs))), all_orfs)
for (orf in all_orfs) pathways[[orf]] <- list()

for (i in seq_along(row_genes2)) {
  pname <- path_names2[[i]]
  pid <- if (!is.na(id_col)) trimws(mat2[i, id_col]) else ""
  if (!nzchar(pid) || !grepl(id_pattern, pid, perl = TRUE)) {
    resolved_id <- unname(catalog_id[pname])
    if (length(resolved_id) && !is.na(resolved_id) && nzchar(resolved_id)) pid <- resolved_id
    else pid <- local_id_map[[pname]]
  }

  purl <- unname(catalog_url[pname])
  if (!length(purl) || is.na(purl) || !nzchar(purl)) {
    purl <- if (grepl(id_pattern, pid, perl = TRUE)) {
      paste0(
        "https://pathway.yeastgenome.org/YEAST/NEW-IMAGE?type=PATHWAY&object=",
        utils::URLencode(pid, reserved = TRUE)
      )
    } else {
      "https://pathway.yeastgenome.org/YEAST/hierarchical.html"
    }
  }

  entry <- list(id = pid, name = pname, url = purl)
  for (orf in row_genes2[[i]]) {
    if (!(orf %in% all_orfs)) next
    if (is.null(pathways[[orf]][[pid]])) pathways[[orf]][[pid]] <- entry
    status[[orf]] <- "annotated"
  }
}

annotated_orfs <- names(status)[vapply(status, identical, logical(1), "annotated")]
noann_orfs <- setdiff(all_orfs, annotated_orfs)

n_pathways <- length(unique(unlist(lapply(pathways[annotated_orfs], names), use.names = FALSE)))

# Abort if the SGD schema cannot be parsed safely.
if (length(annotated_orfs) < 100L) {
  stop(
    "Only ", length(annotated_orfs),
    " ORFs were mapped from biochemical_pathways.tab; this is implausibly low. No cache was written."
  )
}
if (n_pathways < 20L) {
  stop(
    "Only ", n_pathways,
    " pathways were parsed from biochemical_pathways.tab; this is implausibly low. No cache was written."
  )
}

payload <- list(
  cache_format = "MOSAIC_SGD_3",
  pathways = pathways,
  status = status,
  source = list(
    type = "SGD_bulk_biochemical_pathways",
    schema = if (isTRUE(canonical_schema)) "canonical_5col_pathway_enzyme_ec_gene_reference" else "heuristic",
    location = bulk$source,
    downloaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    rows_parsed = nrow(mat),
    pathway_column = path_col,
    gene_columns = gene_cols,
    pathway_id_column = if (is.na(id_col)) NA_integer_ else id_col
  )
)

saveRDS(payload, out_file)

status_file <- sub("\\.rds$", "_status.tsv", out_file, ignore.case = TRUE)
if (identical(status_file, out_file)) status_file <- paste0(out_file, "_status.tsv")
status_tab <- data.frame(
  SGD_ID = all_orfs,
  status = vapply(status, as.character, character(1)),
  n_pathways = vapply(pathways, length, integer(1)),
  stringsAsFactors = FALSE
)
write.table(status_tab, status_file, sep = "\t", quote = FALSE, row.names = FALSE)

membership_file <- sub("\\.rds$", "_membership.tsv", out_file, ignore.case = TRUE)
if (identical(membership_file, out_file)) membership_file <- paste0(out_file, "_membership.tsv")
membership <- do.call(rbind, lapply(annotated_orfs, function(orf) {
  pws <- pathways[[orf]]
  if (!length(pws)) return(NULL)
  data.frame(
    SGD_ID = orf,
    pathway_id = vapply(pws, `[[`, character(1), "id"),
    pathway_name = vapply(pws, `[[`, character(1), "name"),
    pathway_url = vapply(pws, `[[`, character(1), "url"),
    stringsAsFactors = FALSE
  )
}))
if (!is.null(membership) && nrow(membership)) {
  write.table(membership, membership_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\nBulk SGD cache created successfully\n")
cat("-----------------------------------\n")
cat(sprintf("source rows parsed : %d\n", nrow(mat)))
cat(sprintf("pathways parsed    : %d\n", n_pathways))
cat(sprintf("annotated ORFs     : %d\n", length(annotated_orfs)))
cat(sprintf("no annotation ORFs : %d\n", length(noann_orfs)))
cat(sprintf("request failed     : 0\n"))
cat(sprintf("cache              : %s\n", normalizePath(out_file, mustWork = FALSE)))
cat(sprintf("status table       : %s\n", normalizePath(status_file, mustWork = FALSE)))
if (!is.null(membership) && nrow(membership)) {
  cat(sprintf("membership table   : %s\n", normalizePath(membership_file, mustWork = FALSE)))
}
cat("\nNo per-gene YeastPathways requests were used.\n")
