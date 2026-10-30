options(warn = 1)  # Print each warning immediately so it is preserved in the stage log.
# Build reports

gsea_dir          <- "GSEA_results/"
results_dir       <- "results/"
output_dir        <- "heatmaps/"

organism_config_file <- "organism_config.R"
if (!file.exists(organism_config_file))
  stop(sprintf("Missing organism configuration file: %s", organism_config_file))
source(organism_config_file, local = FALSE)

if (!exists("map_gene_labels", inherits = FALSE)) {
  map_gene_labels <- !identical(symbol_column, key_type)
}
if (!exists("enable_sgd_gene_links", inherits = FALSE)) {
  enable_sgd_gene_links <- identical(orgdb, "org.Sc.sgd.db")
}
if (!exists("enable_sgd_pathway_lookup", inherits = FALSE)) {
  enable_sgd_pathway_lookup <- identical(orgdb, "org.Sc.sgd.db")
}

reference_level   <- NULL
contrast_factor  <- "treatment"

sig_threshold     <- 0.05
do_similarity_reduction <- TRUE
similarity_threshold    <- 0.375

similarity_metric       <- "jaccard_overlap"
max_nes_scale     <- NULL
sgd_cache_file    <- "saccharomyces_cerevisiae_GO_and_KEGG_data/sgd_pathway_cache_sce.rds"

filter_mode <- "interpretation"

winner_score <- "score_weighted"


suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(tools)
  library(AnnotationDbi)
  library(orgdb, character.only = TRUE)
})

if (!requireNamespace(orgdb, quietly = TRUE))
  stop(sprintf("OrgDb package '%s' is not installed.", orgdb))
orgdb_obj <- getExportedValue(orgdb, orgdb)
available_keytypes <- AnnotationDbi::keytypes(orgdb_obj)
available_columns  <- AnnotationDbi::columns(orgdb_obj)
if (!key_type %in% available_keytypes)
  stop(sprintf("key_type '%s' is not supported by %s.", key_type, orgdb))
if (!symbol_column %in% available_columns) symbol_column <- key_type
is_sgd_yeast <- identical(orgdb, "org.Sc.sgd.db")
enable_sgd_gene_links <- isTRUE(enable_sgd_gene_links) && is_sgd_yeast
enable_sgd_pathway_lookup <- isTRUE(enable_sgd_pathway_lookup) && is_sgd_yeast

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

build_variance_partition_report <- function() {
  stale_outputs <- file.path(output_dir, c("Variance_partition.png",
                                           "variance_partition_gene_wise.tsv",
                                           "variance_partition_summary.tsv",
                                           "variance_partition_status.tsv"))
  unlink(stale_outputs[file.exists(stale_outputs)])
  gene_file <- file.path(results_dir, "variance_partition_gene_wise.tsv")
  summary_file <- file.path(results_dir, "variance_partition_summary.tsv")
  status_file <- file.path(results_dir, "variance_partition_status.tsv")

  if (!file.exists(gene_file)) {
    if (file.exists(status_file)) {
      file.copy(status_file, file.path(output_dir, basename(status_file)), overwrite=TRUE)
      st <- tryCatch(read.table(status_file, header=TRUE, sep="\t", quote="", comment.char="",
                                stringsAsFactors=FALSE, check.names=FALSE), error=function(e) NULL)
      if (!is.null(st) && all(c("key", "value") %in% colnames(st))) {
        sm <- setNames(st$value, st$key)
        cat(sprintf("Variance partition report: %s%s\n",
                    if ("status" %in% names(sm)) sm[["status"]] else "not available",
                    if ("note" %in% names(sm) && nzchar(sm[["note"]])) paste0(" - ", sm[["note"]]) else ""))
      } else {
        cat("Variance partition report: gene-wise table not found; skipping.\n")
      }
    } else {
      cat("Variance partition report: gene-wise table not found; skipping.\n")
    }
    return(invisible(NULL))
  }

  vp <- read.table(gene_file, header=TRUE, sep="\t", quote="", comment.char="",
                   stringsAsFactors=FALSE, check.names=FALSE)
  if (!nrow(vp) || !"gene" %in% colnames(vp)) {
    warning("Variance partition gene-wise table is empty or malformed; plot was not generated.")
    return(invisible(NULL))
  }
  component_cols <- setdiff(colnames(vp), c("gene", "Model_explained"))
  numeric_cols <- component_cols[vapply(vp[, component_cols, drop=FALSE], is.numeric, logical(1))]
  if (!length(numeric_cols)) {
    warning("Variance partition table contains no numeric variance components; plot was not generated.")
    return(invisible(NULL))
  }

  complete_rows <- apply(vp[, numeric_cols, drop=FALSE], 1, function(z) all(is.finite(z)))
  if (any(!complete_rows)) {
    warning(sprintf("Omitting %d gene(s) with incomplete variance fractions from the variance plot.",
                    sum(!complete_rows)))
    vp <- vp[complete_rows, , drop=FALSE]
  }
  if (!nrow(vp)) {
    warning("No complete gene-wise variance fractions are available; plot was not generated.")
    return(invisible(NULL))
  }

  for (src in c(gene_file, summary_file, status_file)) {
    if (file.exists(src)) file.copy(src, file.path(output_dir, basename(src)), overwrite=TRUE)
  }

  means <- vapply(numeric_cols, function(nm) mean(vp[[nm]]), numeric(1))
  med <- vapply(numeric_cols, function(nm) median(vp[[nm]]), numeric(1))
  q25 <- vapply(numeric_cols, function(nm) quantile(vp[[nm]], 0.25, names=FALSE), numeric(1))
  q75 <- vapply(numeric_cols, function(nm) quantile(vp[[nm]], 0.75, names=FALSE), numeric(1))

  last_components <- intersect(c("Measurement.error", "Residuals"), numeric_cols)
  modeled_components <- setdiff(numeric_cols, last_components)
  modeled_components <- modeled_components[order(means[modeled_components], decreasing=TRUE)]
  numeric_cols <- c(modeled_components, last_components)
  means <- means[numeric_cols]
  med <- med[numeric_cols]
  q25 <- q25[numeric_cols]
  q75 <- q75[numeric_cols]

  mean_sum <- sum(means)
  if (abs(mean_sum - 1) > 1e-6) {
    warning(sprintf("Mean variance fractions sum to %.8f instead of one.", mean_sum))
  }

  base_palette <- c("#8AAED6", "#E5A45B", "#7FB77E", "#9C8BC4",
                    "#70B7B3", "#D2B45C", "#A8866B", "#78909C")
  fills <- rep(base_palette, length.out=length(numeric_cols))
  fills[numeric_cols == "Residuals"] <- "#B8B8B8"
  fills[numeric_cols == "Measurement.error"] <- "#D4D4D4"

  pretty_label <- function(x) {
    x <- gsub("\\.", " ", x)
    x <- gsub(":", " × ", x, fixed=TRUE)
    x
  }
  labs <- vapply(numeric_cols, pretty_label, character(1))
  ncomp <- length(numeric_cols)
  width_cm <- max(17, 2.2*ncomp + 6)
  height_cm <- 13

  caption_text <- paste(
    "Bars = mean across genes (component means sum to 1);",
    "points = individual genes; diamonds = median;",
    "central whiskers = interquartile range (25th-75th percentile)."
  )
  # Wrap the variance caption to the device width.
  caption_width_chars <- max(68L, floor(width_cm * 4.5))
  caption_lines <- strwrap(caption_text, width=caption_width_chars)
  caption_base_line <- max(5.1, 0.42*max(nchar(labs)) + 2.3)
  bottom_margin <- max(
    7,
    0.42*max(nchar(labs)) + 4,
    caption_base_line + 0.85*(length(caption_lines) - 1L) + 1.4
  )

  png_path <- file.path(output_dir, "Variance_partition.png")
  png(png_path, width=width_cm, height=height_cm, units="cm", res=600, type="cairo-png")
  oldpar <- par(no.readonly=TRUE)
  on.exit({par(oldpar); dev.off()}, add=TRUE)
  par(mar=c(bottom_margin, 5.2, 3.2, 1.2), xpd=FALSE)

  # Draw the grid before bars and overlays.
  bar_width <- 0.82
  bar_space <- 0.35
  bp <- barplot(means, names.arg=labs, las=2, ylim=c(0, 1),
                width=bar_width, space=bar_space,
                col=NA, border=NA, ylab="Fraction of expression variance",
                main="Gene-wise variance partition", cex.names=0.82, cex.axis=0.9,
                cex.lab=0.95, cex.main=1.05)
  abline(h=seq(0, 1, by=0.1), col="#E8E8E8", lwd=0.7)
  rect(bp - bar_width/2, 0, bp + bar_width/2, means,
       col=fills, border="#666666", lwd=0.8)
  set.seed(104729)
  for (i in seq_along(numeric_cols)) {
    vals <- vp[[numeric_cols[i]]]
    vals <- vals[is.finite(vals)]
    if (!length(vals)) next
    xx <- bp[i] + runif(length(vals), -0.27, 0.27)
    points(xx, vals, pch=16, cex=0.23, col=grDevices::adjustcolor("#222222", alpha.f=0.12))
    segments(bp[i], q25[i], bp[i], q75[i], lwd=2.2, col="#333333")
    segments(bp[i]-0.11, q25[i], bp[i]+0.11, q25[i], lwd=1.5, col="#333333")
    segments(bp[i]-0.11, q75[i], bp[i]+0.11, q75[i], lwd=1.5, col="#333333")
    points(bp[i], med[i], pch=18, cex=0.7, col="#111111")
  }
  box(col="#666666")
  for (j in seq_along(caption_lines)) {
    mtext(caption_lines[j], side=1,
          line=caption_base_line + 0.85*(j - 1L),
          cex=0.68, col="#555555")
  }
  dev.off()
  on.exit(NULL, add=FALSE)
  cat(sprintf("Variance partition: %s\n", basename(png_path)))
  invisible(NULL)
}

.manifest_path <- file.path(results_dir, "analysis_manifest.tsv")
pipeline_info  <- if (file.exists(.manifest_path)) {
  mi <- tryCatch(read.table(.manifest_path, header=TRUE, sep="\t",
                            stringsAsFactors=FALSE, quote="", comment.char="",
                            check.names=FALSE), error=function(e) NULL)
  if (!is.null(mi)) setNames(mi$value, mi$key) else character(0)
} else character(0)

pipeline_interaction_mode <- if ("analysis_mode" %in% names(pipeline_info)) {
  grepl("interaction", pipeline_info[["analysis_mode"]], ignore.case=TRUE)
} else {
  FALSE
}

int_factors_from_manifest <- if ("interaction_factors" %in% names(pipeline_info) &&
                                 nzchar(pipeline_info[["interaction_factors"]])) {
  strsplit(pipeline_info[["interaction_factors"]], ";")[[1]]
} else {
  character(0)
}


load_tsv <- function(path) {
  tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               quote = "", comment.char = "", check.names = FALSE),
    error = function(e) NULL
  )
}

detect_reference <- function(nms) {
  rhs <- vapply(nms, function(n) {
    p <- strsplit(n, "_vs_")[[1]]
    if (length(p) < 2) NA_character_ else paste(p[-1], collapse = "_vs_")
  }, character(1))
  rhs <- rhs[!is.na(rhs)]
  names(which.max(table(rhs)))
}

htmlesc <- function(s) {
  s <- gsub("&", "&amp;",  s, fixed = TRUE)
  s <- gsub("<", "&lt;",   s, fixed = TRUE)
  s <- gsub(">", "&gt;",   s, fixed = TRUE)
  s <- gsub('"', "&quot;", s, fixed = TRUE)
  gsub("'", "&#39;", s, fixed = TRUE)
}

jesc <- function(s) {
  if (is.null(s) || length(s) == 0) return("")
  s <- as.character(s); s[is.na(s)] <- ""
  s <- chartr("\\", "/", s)
  s <- gsub('"',  '\\"',  s, fixed = TRUE)
  s <- gsub("</", "<\\/", s, fixed = TRUE)
  s <- gsub("\n", "\\n",  s, fixed = TRUE)
  gsub("\r", "\\r", s, fixed = TRUE)
}

display_condition <- function(t) as.character(t)

split_condition_components <- function(x) {
  if (grepl("×", x, fixed=TRUE)) strsplit(x, "×", fixed=TRUE)[[1]]
  else as.character(x)
}

condition_component_order <- if ("condition_component_order" %in% names(pipeline_info) &&
                                   nzchar(pipeline_info[["condition_component_order"]])) {
  strsplit(pipeline_info[["condition_component_order"]], ";", fixed=TRUE)[[1]]
} else c(contrast_factor, int_factors_from_manifest)

manifest_header_factor_order <- if ("header_factor_order" %in% names(pipeline_info) &&
                                      nzchar(pipeline_info[["header_factor_order"]])) {
  strsplit(pipeline_info[["header_factor_order"]], ";", fixed=TRUE)[[1]]
} else condition_component_order

background_header_order <- setdiff(manifest_header_factor_order, contrast_factor)
missing_background <- setdiff(condition_component_order,
                              c(contrast_factor, background_header_order))
header_factor_order <- intersect(
  unique(c(contrast_factor, background_header_order, missing_background)),
  condition_component_order
)

order_parts_for_header <- function(parts) {
  vals <- setNames(rep("", length(condition_component_order)), condition_component_order)
  n <- min(length(parts), length(condition_component_order))
  if (n > 0L) vals[condition_component_order[seq_len(n)]] <- parts[seq_len(n)]
  unname(vals[header_factor_order])
}

jstr <- function(s) {
  if (is.null(s) || length(s) == 0 || (length(s)==1 && is.na(s))) return("null")
  paste0('"', jesc(as.character(s)), '"')
}

nes_to_hex <- function(nes, max_abs) {
  if (is.na(nes)) return("#D0D0D0")
  t <- max(-1, min(1, nes / max_abs))
  anchors <- c(-1, -0.55, 0, 0.55, 1)
  cols <- c("#3F6F9F", "#86ACCF", "#F7F7F7", "#E49482", "#B34F55")
  rgbm <- grDevices::col2rgb(cols)
  k <- if (t <= anchors[2]) 1L else if (t <= anchors[3]) 2L else if (t <= anchors[4]) 3L else 4L
  u <- (t - anchors[k]) / (anchors[k+1] - anchors[k])
  rgbv <- round(rgbm[,k] + u * (rgbm[,k+1] - rgbm[,k]))
  sprintf("#%02X%02X%02X", rgbv[1], rgbv[2], rgbv[3])
}

cfn_svg <- function(nes, sig, max_abs) {
  if (is.na(nes)) return("#E0E0E0")
  if (!is.na(sig) && sig == "ns") return("#D0D0D0")
  nes_to_hex(nes, max_abs)
}

metric_label <- function(metric) {
  switch(metric,
         jaccard_overlap = "Jaccard+overlap",
         combined = "Jaccard+overlap",
         jaccard  = "Jaccard similarity",
         overlap  = "overlap coefficient",
         dice     = "Sørensen-Dice similarity",
         stop("Unsupported similarity_metric: ", metric)
  )
}

compute_set_similarity <- function(a, b, metric = similarity_metric) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(NA_real_)
  if (!length(a) || !length(b)) return(0)
  ni <- length(intersect(a, b))
  na <- length(a); nb <- length(b)
  jac <- ni / length(union(a, b))
  ovl <- ni / min(na, nb)
  dice <- 2 * ni / (na + nb)
  switch(metric,
         jaccard_overlap = mean(c(jac, ovl)),
         combined = mean(c(jac, ovl)),
         jaccard = jac,
         overlap = ovl,
         dice = dice,
         stop("Unsupported similarity_metric: ", metric)
  )
}

cluster_from_similarity <- function(sim, selected_ids) {
  sim <- as.matrix(sim)
  if (is.null(rownames(sim)) || is.null(colnames(sim)))
    stop("similarity matrix has no term identifiers")
  if (length(selected_ids) < 2L) stop("too few selected terms for clustering")

  # Clustering orders terms but never changes the selected set.
  common <- selected_ids[selected_ids %in% rownames(sim) & selected_ids %in% colnames(sim)]
  if (!length(common)) stop("similarity matrix contains none of the selected term identifiers")

  full_sim <- matrix(0, length(selected_ids), length(selected_ids),
                     dimnames=list(selected_ids, selected_ids))
  full_sim[common, common] <- sim[common, common, drop=FALSE]
  full_sim[!is.finite(full_sim)] <- 0
  full_sim[full_sim < 0] <- 0
  full_sim[full_sim > 1] <- 1
  full_sim <- (full_sim + t(full_sim)) / 2
  diag(full_sim) <- 1

  distance_matrix <- 1 - full_sim
  diag(distance_matrix) <- 0
  dimnames(distance_matrix) <- list(selected_ids, selected_ids)
  list(
    hc = hclust(as.dist(distance_matrix), method="average"),
    sel = selected_ids,
    missing_from_similarity = setdiff(selected_ids, common)
  )
}



.clean_pw_name <- function(raw) {
  nm <- trimws(raw)
  if (grepl(">", nm, fixed=TRUE)) nm <- trimws(strsplit(nm,">",fixed=TRUE)[[1]][[ length(strsplit(nm,">",fixed=TRUE)[[1]]) ]])
  nm <- trimws(gsub("\\|\\s*$","",nm))
  if (!nchar(nm)) nm <- raw; nm
}

if (!exists(".SGD_CACHE")) {
  .SGD_CACHE <- list()
  .SGD_STATUS <- list()
  .SGD_CACHE_AVAILABLE <- FALSE
  .SGD_CACHE_SOURCE <- NULL

  if (enable_sgd_pathway_lookup && file.exists(sgd_cache_file)) {
    cached_sgd <- readRDS(sgd_cache_file)
    cache_fmt_ok <- is.list(cached_sgd) && !is.null(cached_sgd$cache_format) &&
      cached_sgd$cache_format %in% c("MOSAIC_SGD_2", "MOSAIC_SGD_3")
    cache_names <- if (cache_fmt_ok && is.list(cached_sgd$pathways)) unique(unlist(lapply(cached_sgd$pathways, function(g) {
      if (is.null(g) || !length(g)) return(character(0))
      vapply(g, function(e) if (is.list(e) && !is.null(e$name)) as.character(e$name[[1]]) else "", character(1))
    }), use.names=FALSE)) else character(0)
    cache_names <- trimws(cache_names[nzchar(trimws(cache_names))])
    cache_not_reference_column <- length(cache_names) >= 10L &&
      mean(grepl("(?:PMID|SGD_REF):", cache_names, ignore.case=TRUE, perl=TRUE)) < 0.10
    if (cache_fmt_ok && cache_not_reference_column) {
      .SGD_CACHE <- cached_sgd$pathways
      .SGD_STATUS <- cached_sgd$status
      .SGD_CACHE_SOURCE <- cached_sgd$source
      .SGD_CACHE_AVAILABLE <- TRUE
    } else {
      warning("The SGD pathway cache is unsupported or malformed (e.g. reference IDs stored as pathway names); SGD pathway annotations will be skipped.", call.=FALSE)
    }
    if (.SGD_CACHE_AVAILABLE) {
    }
  } else if (enable_sgd_pathway_lookup) {
    cat("No local SGD biochemical pathway cache found; SGD pathway annotations will be skipped.\n")
  }
}
if (!exists(".SGD_STATUS")) .SGD_STATUS <- list()
if (!exists(".SGD_CACHE_AVAILABLE")) .SGD_CACHE_AVAILABLE <- FALSE

fetch_sgd_pathways <- function(orf) {
  if (!enable_sgd_pathway_lookup || !isTRUE(.SGD_CACHE_AVAILABLE)) return(list())
  pw <- .SGD_CACHE[[orf]]
  if (is.null(pw)) {
    if (is.null(.SGD_STATUS[[orf]])) .SGD_STATUS[[orf]] <<- "not_in_cache"
    return(list())
  }
  pw
}

.anc_cache_global <- NULL

build_anc_cache <- function() {
  tryCatch({
    suppressPackageStartupMessages(library(GO.db))
    as.list(GOBPANCESTOR)
  }, error=function(e) NULL)
}

find_lca <- function(go_ids, anc_cache) {
  go_ids <- unique(go_ids[nchar(go_ids)>0])
  if (!length(go_ids)) return(NA_character_)
  go_roots <- c("GO:0008150", "GO:0003674", "GO:0005575",
                "GO:0009987", "GO:0065007", "GO:0044699")
  if (length(go_ids)==1 && !go_ids[1] %in% go_roots) return(go_ids[1])
  anc_list <- lapply(go_ids, function(g) {
    a <- anc_cache[[g]]
    if (is.null(a)) character(0) else setdiff(a[a!="all"], go_roots)
  })
  common <- Reduce(intersect, anc_list)
  common <- setdiff(common, go_roots)
  if (!length(common)) return(NA_character_)
  depths <- vapply(common, function(id) {
    a <- anc_cache[[id]]
    if (is.null(a)) 0L else length(setdiff(a[a!="all"], go_roots))
  }, integer(1))
  common[which.max(depths)]
}

go_term_label <- function(go_id) {
  go_roots <- c("GO:0008150","GO:0003674","GO:0005575",
                "GO:0009987","GO:0065007","GO:0044699")
  if (go_id %in% go_roots) return("")
  tryCatch({
    suppressPackageStartupMessages(library(GO.db))
    t <- Term(GOTERM[[go_id]])
    if (is.null(t)||is.na(t)) go_id else paste0(t," (",go_id,")")
  }, error=function(e) go_id)
}

find_cluster_category <- function(go_ids, anc_cache) {
  go_ids <- unique(go_ids[nzchar(go_ids)])
  if (!length(go_ids) || is.null(anc_cache)) return(NA_character_)
  roots <- c("GO:0008150", "GO:0003674", "GO:0005575",
             "GO:0009987", "GO:0065007", "GO:0044699", "all")
  anc_sets <- lapply(go_ids, function(g) {
    a <- anc_cache[[g]]
    unique(setdiff(c(g, if (is.null(a)) character(0) else a), roots))
  })
  candidates <- unique(unlist(anc_sets, use.names=FALSE))
  if (!length(candidates)) return(NA_character_)
  coverage <- vapply(candidates, function(id)
    sum(vapply(anc_sets, function(z) id %in% z, logical(1))), integer(1))
  depth <- vapply(candidates, function(id) {
    a <- anc_cache[[id]]
    if (is.null(a)) 0L else length(setdiff(a, roots))
  }, integer(1))
  candidates[order(-coverage, -depth, candidates)][1]
}


format_heatmap_value <- function(nes, sig, cld = "") {
  cld <- if (is.na(cld)) "" else trimws(as.character(cld))
  if (!is.finite(nes)) return(if (nzchar(cld)) paste0("ns ", cld) else "ns")
  if (is.na(sig) || sig != "*")
    return(paste0("ns (", sprintf("%.3f", nes), ")", if (nzchar(cld)) paste0(" ", cld) else ""))
  paste0(sprintf("%.3f", nes), if (nzchar(cld)) paste0(" ", cld) else "")
}

write_analysis_workbook <- function(path, sheets) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path, format_headers = TRUE)
    return(TRUE)
  }
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    header_style <- openxlsx::createStyle(
      fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold",
      halign = "center", valign = "center", border = "Bottom"
    )
    for (nm in names(sheets)) {
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, sheets[[nm]], withFilter = TRUE)
      if (ncol(sheets[[nm]]) > 0) {
        openxlsx::addStyle(wb, nm, header_style, rows = 1,
                           cols = seq_len(ncol(sheets[[nm]])), gridExpand = TRUE)
        openxlsx::freezePane(wb, nm, firstRow = TRUE)
        openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
      }
    }
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    return(TRUE)
  }
  warning("Excel export skipped: install either 'writexl' or 'openxlsx'.")
  FALSE
}


process_type <- function(type_suffix) {
  use_cld <- FALSE
  
  cat(sprintf("%s: ", type_suffix))
  
  summ_file <- file.path(gsea_dir, paste0(type_suffix,"_treatment_summary.tsv"))
  if (!file.exists(summ_file)) { cat("  Not found:", basename(summ_file),"\n"); return(invisible(NULL)) }
  summ <- load_tsv(summ_file)
  if (is.null(summ)||nrow(summ)==0) { cat("  Empty summary\n"); return(invisible(NULL)) }
  
  nes_cols   <- grep("_NES$",        colnames(summ), value=TRUE)
  sig_cols   <- grep("_sig_vs_ref$", colnames(summ), value=TRUE)
  treatments <- sub("_NES$","",nes_cols)
  n_treat    <- length(treatments)
  
  cat(sprintf("%d term(s), %d condition(s)\n", nrow(summ), length(treatments)))
  use_cld    <- FALSE
  
  nes_raw <- as.matrix(summ[,nes_cols,drop=FALSE])
  nes_mat <- matrix(
    suppressWarnings(as.numeric(nes_raw)),
    nrow=nrow(summ), ncol=length(nes_cols),
    dimnames=list(summ$ID, treatments)
  )
  
  sig_mat <- as.matrix(summ[,sig_cols,drop=FALSE])
  rownames(sig_mat) <- summ$ID; colnames(sig_mat) <- treatments
  
  nes_mat_all <- nes_mat
  sig_mat_all <- sig_mat
  
  cld_mat <- NULL
  cld_mat_all <- NULL
  
  desc_lu      <- setNames(summ$Description, summ$ID)
  desc_lu_original <- desc_lu

  elimination_rows <- list()
  elimination_i <- 0L
  add_eliminated_terms <- function(ids, reason, representative = NA_character_,
                                   similarity = NA_real_) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) return(invisible(NULL))
    reps <- rep_len(as.character(representative), length(ids))
    sims <- rep_len(as.numeric(similarity), length(ids))
    for (ii in seq_along(ids)) {
      gid <- ids[ii]
      elimination_i <<- elimination_i + 1L
      d <- desc_lu_original[gid]
      if (!length(d) || is.na(d) || !nzchar(trimws(d))) d <- gid
      elimination_rows[[elimination_i]] <<- data.frame(
        ID = gid,
        Description = unname(d),
        Removal_reason = reason,
        Representative_ID = ifelse(is.na(reps[ii]), "", reps[ii]),
        Similarity_to_representative = sims[ii],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    invisible(NULL)
  }

  build_elimination_report <- function() {
    out <- if (length(elimination_rows)) do.call(rbind, elimination_rows) else
      data.frame(ID=character(0), Description=character(0),
                 Removal_reason=character(0), Representative_ID=character(0),
                 Similarity_to_representative=numeric(0),
                 stringsAsFactors=FALSE, check.names=FALSE)
    out <- out[!duplicated(out[, c("ID", "Removal_reason")]), , drop=FALSE]
    rownames(out) <- NULL
    if (nrow(out) > 0L) cat(sprintf("  Annotation terms removed before reporting: %d\n", nrow(out)))
    out
  }

  if (type_suffix == "GO_BP") {
    desc_lu <- as.character(desc_lu)
    names(desc_lu) <- summ$ID

    explicit_obsolete <- !is.na(desc_lu) &
      grepl("^obsolete\\b", trimws(desc_lu), ignore.case=TRUE)

    bad_desc <- is.na(desc_lu) | !nzchar(trimws(desc_lu)) |
      trimws(desc_lu) == names(desc_lu) |
      trimws(desc_lu) == paste0(names(desc_lu), " (", names(desc_lu), ")")

    if (any(bad_desc) && requireNamespace("GO.db", quietly=TRUE)) {
      ids_bad <- names(desc_lu)[bad_desc]
      recovered <- vapply(ids_bad, function(go_id) {
        tryCatch({
          obj <- GO.db::GOTERM[[go_id]]
          if (is.null(obj)) return(NA_character_)
          val <- AnnotationDbi::Term(obj)
          if (!length(val) || is.na(val[1]) || !nzchar(trimws(val[1]))) NA_character_
          else as.character(val[1])
        }, error=function(e) NA_character_)
      }, character(1))
      ok <- !is.na(recovered) & nzchar(trimws(recovered)) &
        !grepl("^obsolete\\b", trimws(recovered), ignore.case=TRUE)
      desc_lu[ids_bad[ok]] <- recovered[ok]
    }

    unresolved <- is.na(desc_lu) | !nzchar(trimws(desc_lu)) |
      trimws(desc_lu) == names(desc_lu) |
      trimws(desc_lu) == paste0(names(desc_lu), " (", names(desc_lu), ")")

    obsolete_go_ids <- names(desc_lu)[explicit_obsolete | unresolved]

    if (length(obsolete_go_ids)) {
      removed_now <- intersect(obsolete_go_ids, rownames(nes_mat))
      keep_go_ids <- setdiff(rownames(nes_mat), obsolete_go_ids)
      cat(sprintf("  Removed obsolete/unresolved GO terms: %d\n", length(removed_now)))
      if (length(removed_now)) {
        add_eliminated_terms(removed_now, "obsolete_or_unresolved_GO")
      }

      summ <- summ[summ$ID %in% keep_go_ids, , drop=FALSE]
      nes_mat <- nes_mat[keep_go_ids, , drop=FALSE]
      sig_mat <- sig_mat[keep_go_ids, , drop=FALSE]
      nes_mat_all <- nes_mat_all[keep_go_ids, , drop=FALSE]
      sig_mat_all <- sig_mat_all[keep_go_ids, , drop=FALSE]
      if (!is.null(cld_mat)) cld_mat <- cld_mat[keep_go_ids, , drop=FALSE]
      if (!is.null(cld_mat_all)) cld_mat_all <- cld_mat_all[keep_go_ids, , drop=FALSE]
      desc_lu <- desc_lu[keep_go_ids]
    }

    if (nrow(summ)) summ$Description <- unname(desc_lu[summ$ID])
  }

  if (type_suffix == "KEGG") {
    desc_lu <- as.character(desc_lu)
    names(desc_lu) <- summ$ID
    clean_desc <- trimws(desc_lu)
    dup_any <- duplicated(clean_desc) | duplicated(clean_desc, fromLast=TRUE)
    if (any(dup_any)) {
      # Keep distinct KEGG IDs even when descriptions match.
      ids_dup <- names(desc_lu)[dup_any]
      desc_lu[ids_dup] <- paste0(clean_desc[dup_any], " [", ids_dup, "]")
      cat(sprintf("  KEGG display labels disambiguated with pathway IDs: %d pathway(s)\n", length(ids_dup)))
    }
    summ$Description <- unname(desc_lu[summ$ID])
  }
  
  int_x_cols    <- grep("^interaction_x_", colnames(summ), value=TRUE)
  int_3way_cols <- grep("^interaction_3way", colnames(summ), value=TRUE)
  has_3way_col  <- length(int_3way_cols) > 0L
  n_int_cols    <- length(int_x_cols) + length(int_3way_cols)
  
  n_parts_in_treatments <- if (length(treatments) > 0)
    max(vapply(lapply(treatments, split_condition_components), length, integer(1)))
  else 1L
  
  interaction_mode <- pipeline_interaction_mode ||
    "interaction_sig" %in% colnames(summ) ||
    n_int_cols > 0L ||
    (n_parts_in_treatments >= 3)
  
  interaction_3way <- has_3way_col ||
    (interaction_mode && n_parts_in_treatments >= 3)

  header_parts <- lapply(treatments, function(t)
    order_parts_for_header(split_condition_components(t)))
  if (interaction_mode && length(header_parts) && length(header_factor_order) >= 2L) {
    ord_df <- as.data.frame(do.call(rbind, header_parts), stringsAsFactors=FALSE)
    ord_idx <- do.call(order, c(ord_df, list(na.last=TRUE)))
    treatments <- treatments[ord_idx]
    header_parts <- header_parts[ord_idx]
    nes_mat <- nes_mat[, treatments, drop=FALSE]
    sig_mat <- sig_mat[, treatments, drop=FALSE]
    nes_mat_all <- nes_mat_all[, treatments, drop=FALSE]
    sig_mat_all <- sig_mat_all[, treatments, drop=FALSE]
    n_treat <- length(treatments)
    }
  names(header_parts) <- treatments
  
  if (interaction_mode && n_int_cols == 0L)
    cat("  NOTE: interaction mode from treatment names but no interaction columns in summary.\n",
        "        Interaction significance columns are required in the treatment summary.\n")
  
  cld_cols <- paste0(treatments, "_CLD")
  use_cld <- (!interaction_mode && n_treat >= 2) ||
    (interaction_mode && all(cld_cols %in% colnames(summ)))
  if (use_cld) {
    missing_cld_cols <- setdiff(cld_cols, colnames(summ))
    if (length(missing_cld_cols) == 0L) {
      cld_mat <- as.matrix(summ[, cld_cols, drop = FALSE])
      rownames(cld_mat) <- summ$ID
      colnames(cld_mat) <- treatments
      cld_mat_all <- cld_mat
    } else {
      use_cld <- FALSE
      cld_mat <- NULL
      cld_mat_all <- NULL
      cat(sprintf(
        "  CLD columns missing: %s\n  The treatment summary does not contain the required CLD data.\n",
        paste(missing_cld_cols, collapse = ", ")
      ))
    }
  }
  
  int_label_prefix <- if (length(treatments) > 0L) {
    split_condition_components(treatments[1])[1]
  } else contrast_factor

  parse_int_col <- function(col) {
    if (grepl("^interaction_x_", col)) {
      z <- sub("^interaction_x_", "", col)
      pp <- strsplit(z, "__", fixed=TRUE)[[1]]
      list(kind="2way", factor=pp[1], level=if(length(pp)>=2) pp[2] else NA_character_,
           tested_level=if(length(pp)>=3) pp[3] else NA_character_,
           tested_reference=if(length(pp)>=4) pp[4] else NA_character_)
    } else if (grepl("^interaction_3way", col)) {
      z <- sub("^interaction_3way__?", "", col)
      pp <- strsplit(z, "__", fixed=TRUE)[[1]]
      list(kind="3way", level=if(length(pp)>=1) pp[1] else NA_character_,
           factors=if(length(pp)>=2) pp[-1] else sub("^interaction_x_", "", int_x_cols))
    } else list(kind="unknown")
  }

  condition_cf_levels <- unique(vapply(treatments, function(z)
    split_condition_components(z)[1], character(1)))
  order_interaction_cols <- function(cols) {
    if (!length(cols)) return(cols)
    lev <- vapply(cols, function(z) {
      q <- parse_int_col(z); if(is.null(q$level)) NA_character_ else q$level
    }, character(1))
    fac <- vapply(cols, function(z) {
      q <- parse_int_col(z)
      if(identical(q$kind,"2way")) q$factor else paste(q$factors, collapse="×")
    }, character(1))
    lev_rank <- match(lev, condition_cf_levels)
    lev_rank[is.na(lev_rank)] <- length(condition_cf_levels) + 1L
    background_order <- setdiff(header_factor_order, contrast_factor)
    fac_rank <- match(fac, background_order)
    fac_rank[is.na(fac_rank)] <- length(background_order) + seq_len(sum(is.na(fac_rank)))
    cols[order(lev_rank, fac_rank, na.last=TRUE)]
  }
  int_x_cols <- order_interaction_cols(int_x_cols)
  int_3way_cols <- order_interaction_cols(int_3way_cols)
  all_int_cols <- c(int_x_cols, int_3way_cols)

  int_col_labels_full <- setNames(vapply(all_int_cols, function(col) {
    q <- parse_int_col(col)
    if (identical(q$kind, "2way")) {
      pfx <- if(!is.na(q$level) && nzchar(q$level)) q$level else tolower(int_label_prefix)
      paste0(pfx, " × ", q$factor,
             if (!is.null(q$tested_level) && !is.na(q$tested_level) && nzchar(q$tested_level))
               paste0(" [", q$tested_level,
                      if(!is.null(q$tested_reference) && !is.na(q$tested_reference) && nzchar(q$tested_reference))
                        paste0(" vs ", q$tested_reference) else "", "]") else "")
    } else if (identical(q$kind, "3way")) {
      pfx <- if(!is.na(q$level) && nzchar(q$level)) q$level else tolower(int_label_prefix)
      paste0(pfx, " × ", paste(q$factors, collapse=" × "))
    } else col
  }, character(1)), all_int_cols)

  int_col_labels <- setNames(vapply(all_int_cols, function(col) {
    q <- parse_int_col(col)
    if (identical(q$kind, "2way")) {
      paste0(tolower(contrast_factor), " × ", q$factor)
    } else if (identical(q$kind, "3way")) {
      paste0(tolower(contrast_factor), " × ", paste(q$factors, collapse=" × "))
    } else col
  }, character(1)), all_int_cols)

  if (length(int_col_labels) && anyDuplicated(unname(int_col_labels))) {
    seen_labels <- integer(0)
    names(seen_labels) <- character(0)
    int_col_labels <- setNames(vapply(unname(int_col_labels), function(lbl) {
      n_seen <- if (lbl %in% names(seen_labels)) seen_labels[[lbl]] else 0L
      seen_labels[[lbl]] <<- n_seen + 1L
      if (n_seen == 0L) lbl else paste0(lbl, paste(rep("\u200B", n_seen), collapse=""))
    }, character(1)), names(int_col_labels))
  }

  multi_level_2way <- interaction_mode && !interaction_3way && length(all_int_cols) > 1L
  show_formal_int_cols <- interaction_3way || multi_level_2way

  render_int_header <- function(lbl, xc, h_total) {
    parts    <- trimws(strsplit(lbl, "\u00d7", fixed=TRUE)[[1]])
    all_ln   <- unlist(lapply(seq_along(parts), function(k)
      c(if(k>1) "\u00d7" else NULL, parts[k])))
    n_ln     <- length(all_ln)
    line_h   <- h_total / (n_ln + 1)
    paste(vapply(seq_len(n_ln), function(k) {
      is_sym <- all_ln[k] == "\u00d7"
      paste0("<text x='",xc,"' y='",round(line_h*k),
             "' text-anchor='middle' font-size='",if(is_sym)"9" else "6.5",
             "' fill='white' font-weight='",if(is_sym)"normal" else "bold","'>",
             htmlesc(all_ln[k]),"</text>")
    }, character(1)), collapse="")
  }
  
  format_factor_label <- function(x) {
    x <- trimws(gsub("_", " ", as.character(x), fixed=TRUE))
    if (!nzchar(x)) return("Condition")
    paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
  }
  difference_column_label <- paste0(format_factor_label(contrast_factor), " differences")

  annotation_label <- if (interaction_3way) "Interaction" else
    if (interaction_mode) "Interaction" else difference_column_label
  
  diff_col_vec <- if (length(all_int_cols) > 0) {
    int_mats <- lapply(all_int_cols, function(col) {
      setNames(summ[[col]], summ$ID)
    })
    setNames(vapply(summ$ID, function(gid) {
      if (any(vapply(int_mats, function(v) {
        !is.na(v[gid]) && v[gid] == "*"
      }, logical(1)))) "*" else ""
    }, character(1)), summ$ID)
  } else if ("treatment_differences" %in% colnames(summ)) {
    setNames(summ$treatment_differences, summ$ID)
  } else setNames(rep("", nrow(summ)), summ$ID)
  
  int_sig_vecs <- setNames(lapply(all_int_cols, function(col) {
    setNames(summ[[col]], summ$ID)
  }), all_int_cols)
  
  orig_files <- list.files(gsea_dir, pattern=paste0("_",type_suffix,"\\.tsv$"), full.names=TRUE)
  orig_files <- orig_files[!grepl("_treatment_summary\\.tsv$", orig_files)]
  orig_names <- gsub(paste0("_",type_suffix,"\\.tsv$"),"",basename(orig_files))
  
  if (interaction_mode) {
    ref_level_manifest <- if ("reference_level" %in% names(pipeline_info))
      pipeline_info[["reference_level"]] else NA_character_
    has_vs_idx <- grepl("_vs_", orig_names, fixed=TRUE)
    rhs_first <- vapply(orig_names, function(nm) {
      if (!grepl("_vs_", nm, fixed=TRUE)) return(NA_character_)
      zz <- strsplit(nm, "_vs_", fixed=TRUE)[[1]]
      split_condition_components(zz[2])[1]
    }, character(1))
    is_ref_idx <- has_vs_idx & (is.na(ref_level_manifest) | rhs_first == ref_level_manifest)
    ref_treats  <- vapply(orig_names[is_ref_idx], function(nm)
      strsplit(nm, "_vs_", fixed=TRUE)[[1]][1], character(1))
    is_ref      <- is_ref_idx
    ref_level   <- if (!is.na(ref_level_manifest)) ref_level_manifest else if (length(ref_treats) > 0) {
      split_condition_components(sub("^.*_vs_", "", orig_names[is_ref_idx][1]))[1]
    } else NA_character_
  } else {
    ref_level  <- detect_reference(orig_names)
    ref_sfx    <- paste0("_vs_",ref_level)
    is_ref     <- endsWith(orig_names,ref_sfx)
    ref_treats <- sub(paste0(ref_sfx,"$"),"",orig_names[is_ref])
  }
  
  normalize_term_id <- function(x) {
    x <- as.character(x)
    x <- sub("^\ufeff", "", x)
    trimws(x)
  }
  normalize_condition_key <- function(x) trimws(as.character(x))

  orig_data <- setNames(lapply(orig_files[is_ref], function(f) {
    df <- tryCatch(
      read.delim(f, header=TRUE, stringsAsFactors=FALSE, quote="",
                 comment.char="", check.names=FALSE),
      error=function(e) {
        cat(sprintf("  WARNING: cannot read original GSEA file %s: %s\n",
                    basename(f), conditionMessage(e)))
        NULL
      }
    )
    if (is.null(df)) return(NULL)
    colnames(df) <- sub("^\ufeff", "", colnames(df))
    if (!"ID" %in% colnames(df)) {
      cat(sprintf("  WARNING: no ID column in %s\n", basename(f)))
      return(NULL)
    }
    df$ID <- normalize_term_id(df$ID)
    df <- df[nzchar(df$ID), , drop=FALSE]
    df <- df[!duplicated(df$ID), , drop=FALSE]
    rownames(df) <- df$ID
    df
  }), normalize_condition_key(ref_treats))

  le_cache <- new.env(hash=TRUE, parent=emptyenv())
  for (t0 in names(orig_data)) {
    t <- normalize_condition_key(t0)
    df <- orig_data[[t0]]
    if (is.null(df) || !"core_enrichment" %in% colnames(df)) {
      cat(sprintf("  WARNING: core_enrichment is missing for condition %s.\n", t))
      next
    }
    n_nonempty <- 0L
    for (i in seq_len(nrow(df))) {
      gid <- normalize_term_id(df$ID[i])
      ce  <- trimws(as.character(df$core_enrichment[i]))
      genes <- if (!is.na(ce) && nzchar(ce)) {
        z <- trimws(strsplit(ce, "/", fixed=TRUE)[[1]])
        unique(z[!is.na(z) & nzchar(z)])
      } else character(0)
      if (length(genes)) n_nonempty <- n_nonempty + 1L
      le_cache[[paste0(gid, "|||", t)]] <- genes
    }
    if (n_nonempty < nrow(df)) {
      n_empty <- nrow(df) - n_nonempty
      cat(sprintf(
        "  WARNING: %d of %d terms have an empty leading edge in condition %s.\n",
        n_empty, nrow(df), t
      ))
    }
  }
  get_le <- function(gid, t) {
    key <- paste0(normalize_term_id(gid), "|||", normalize_condition_key(t))
    v <- le_cache[[key]]
    if (is.null(v)) character(0) else v
  }
  
  get_padj_orig <- function(gid,t) {
    df <- orig_data[[t]]
    if (is.null(df)||!gid %in% rownames(df)) return(NA_real_)
    v  <- df[gid,"p.adjust"]; if(is.null(v)||is.na(v)) NA_real_ else as.numeric(v)
  }
  
  
  cld_are_different <- function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    
    a <- gsub("[^[:alpha:]]", "", trimws(as.character(a)))
    b <- gsub("[^[:alpha:]]", "", trimws(as.character(b)))
    
    if (!nzchar(a) || !nzchar(b)) return(FALSE)
    
    letters_a <- strsplit(a, "")[[1]]
    letters_b <- strsplit(b, "")[[1]]
    
    length(intersect(letters_a, letters_b)) == 0L
  }
  
  is_not_significant_anywhere <- function(gid) {
    if (!gid %in% rownames(sig_mat_all)) return(TRUE)
    vals <- as.character(sig_mat_all[gid, treatments, drop=TRUE])
    !any(!is.na(vals) & vals == "*")
  }

  is_interpretable_term <- function(gid) {
    nes_vals <- as.numeric(nes_mat[gid, ])
    sig_vals <- as.character(sig_mat[gid, ])
    valid_nes <- is.finite(nes_vals)
    valid_sig <- !is.na(sig_vals) & nzchar(trimws(sig_vals))
    if (!all(valid_nes) || !all(valid_sig)) return(FALSE)
    
    if (interaction_mode) {
      any_sig <- any(sig_vals == "*")
      if (!any_sig) return(FALSE)

      if (n_int_cols == 0L) return(FALSE)

      pair_matches_col <- function(col, cf_level, fac_name, lev_a, lev_b) {
        q <- parse_int_col(col)
        if (!identical(q$kind, "2way") || !identical(q$factor, fac_name) ||
            is.null(q$level) || is.na(q$level) || !identical(q$level, cf_level)) return(FALSE)
        qa <- if(!is.null(q$tested_level)) q$tested_level else NA_character_
        qb <- if(!is.null(q$tested_reference)) q$tested_reference else NA_character_
        if (!is.na(qa) && !is.na(qb))
          return(setequal(c(qa, qb), c(lev_a, lev_b)))
        if (!is.na(qa)) {
          ref_key <- paste0("int_ref_", fac_name)
          ref_val <- if (ref_key %in% names(pipeline_info)) pipeline_info[[ref_key]] else NA_character_
          return(!is.na(ref_val) && setequal(c(qa, ref_val), c(lev_a, lev_b)))
        }
        FALSE
      }

      col_is_sig <- function(col) {
        v <- int_sig_vecs[[col]][gid]
        !is.null(v) && !is.na(v) && v == "*"
      }

      if (!interaction_3way) {
        cond_parts <- lapply(treatments, split_condition_components)
        cond_cf <- vapply(cond_parts, function(z) z[1], character(1))
        cf_levels <- unique(cond_cf)
        if (!length(cf_levels)) return(FALSE)
        fac_name <- if (length(int_factors_from_manifest)) int_factors_from_manifest[1] else NA_character_
        fac_pos <- match(fac_name, condition_component_order)
        if (is.na(fac_pos)) fac_pos <- 2L

        block_ok <- vapply(cf_levels, function(cf_lvl) {
          idx <- which(cond_cf == cf_lvl)
          block_sig <- sig_vals[idx]
          block_nes <- nes_vals[idx]
          if (!length(idx) || !any(block_sig == "*")) return(FALSE)

          all_sig_same_sign <- all(block_sig == "*") &&
            length(unique(sign(block_nes))) == 1L
          if (all_sig_same_sign) return(TRUE)

          discordant_pairs <- list()
          if (length(idx) >= 2L) {
            for (ii in seq_along(idx)) for (jj in seq_along(idx)) {
              if (ii >= jj) next
              si <- block_sig[ii]; sj <- block_sig[jj]
              ni <- block_nes[ii]; nj <- block_nes[jj]
              discordant <- ((si == "*") != (sj == "*")) ||
                (si == "*" && sj == "*" && sign(ni) != sign(nj))
              if (discordant) discordant_pairs[[length(discordant_pairs)+1L]] <- c(idx[ii], idx[jj])
            }
          }
          if (!length(discordant_pairs) || is.na(fac_name)) return(FALSE)

          all(vapply(discordant_pairs, function(ij) {
            pa <- cond_parts[[ij[1]]]; pb <- cond_parts[[ij[2]]]
            if (length(pa) < fac_pos || length(pb) < fac_pos) return(FALSE)
            lev_a <- pa[fac_pos]; lev_b <- pb[fac_pos]
            matching_cols <- int_x_cols[vapply(int_x_cols, function(col)
              pair_matches_col(col, cf_lvl, fac_name, lev_a, lev_b), logical(1))]
            length(matching_cols) > 0L && any(vapply(matching_cols, col_is_sig, logical(1)))
          }, logical(1)))
        }, logical(1))

        return(all(block_ok))
      }

      n3way_sig <- if (has_3way_col) {
        any(vapply(int_3way_cols, col_is_sig, logical(1)))
      } else FALSE

      if (n3way_sig) return(TRUE)

      cond_parts <- lapply(treatments, split_condition_components)
      n_parts <- length(cond_parts[[1]])
      if (n_parts < 3L) return(TRUE)

      all_covered <- TRUE
      any_discordance <- FALSE
      for (i in seq_along(treatments)) {
        for (j in seq_along(treatments)) {
          if (i >= j) next
          pi <- cond_parts[[i]]; pj <- cond_parts[[j]]
          diff_pos <- which(pi != pj)
          diff_int <- diff_pos[diff_pos >= 2L]
          if (length(diff_int) != 1L || any(diff_pos == 1L)) next

          sig_i <- sig_vals[i]; nes_i <- nes_vals[i]
          sig_j <- sig_vals[j]; nes_j <- nes_vals[j]
          discordant <- ((sig_i=="*") != (sig_j=="*")) ||
            (sig_i=="*" && sig_j=="*" && sign(nes_i) != sign(nes_j))
          if (!discordant) next
          any_discordance <- TRUE

          fac_pos <- diff_int[1]
          fac_name <- if (fac_pos <= length(condition_component_order))
            condition_component_order[fac_pos] else NA_character_
          cf_level_here <- pi[1]
          if (is.na(fac_name)) { all_covered <- FALSE; next }
          candidate_cols <- int_x_cols[vapply(int_x_cols, function(col)
            pair_matches_col(col, cf_level_here, fac_name, pi[fac_pos], pj[fac_pos]), logical(1))]
          if (!length(candidate_cols) || !any(vapply(candidate_cols, col_is_sig, logical(1))))
            all_covered <- FALSE
        }
      }

      if (!any_discordance) {
        sig_idx <- which(sig_vals == "*")
        return(length(sig_idx) > 0L && length(unique(sign(nes_vals[sig_idx]))) == 1L)
      }
      return(all_covered)
    }

    if (n_treat == 1L) return(sig_vals[1] == "*")

    if (n_treat == 2L) {
      sig1 <- sig_vals[1] == "*"
      sig2 <- sig_vals[2] == "*"

      if (!sig1 && !sig2) return(FALSE)

      same_sign <- sign(nes_vals[1]) == sign(nes_vals[2])

      if (sig1 && sig2 && same_sign) return(TRUE)

      pairwise_supported <- gid %in% names(diff_col_vec) &&
        !is.na(diff_col_vec[gid]) && trimws(diff_col_vec[gid]) == "*"
      return(pairwise_supported)
    }

    if (is.null(cld_mat)) return(FALSE)
    cld_vals  <- as.character(cld_mat[gid, ])
    valid_cld <- !is.na(cld_vals) & nzchar(trimws(cld_vals))
    if (!all(valid_cld)) return(FALSE)
    if (!any(sig_vals == "*")) return(FALSE)

    discordant_pairs <- list()
    for (i in seq_len(n_treat - 1L)) {
      for (j in seq.int(i + 1L, n_treat)) {
        discordant <- ((sig_vals[i] == "*") != (sig_vals[j] == "*")) ||
          (sig_vals[i] == "*" && sig_vals[j] == "*" && sign(nes_vals[i]) != sign(nes_vals[j]))
        if (discordant) discordant_pairs[[length(discordant_pairs) + 1L]] <- c(i, j)
      }
    }
    if (!length(discordant_pairs)) return(TRUE)
    all(vapply(discordant_pairs, function(ij) {
      cld_are_different(cld_vals[ij[1]], cld_vals[ij[2]])
    }, logical(1)))
  }
  
  compute_score <- function(gid) {
    switch(winner_score,
           "score_weighted" = sum(vapply(treatments, function(t) {
             nes <- nes_mat[gid,t]; pv <- get_padj_orig(gid,t)
             if(!is.na(nes)&&!is.na(pv)&&pv>=0) {
               pv_eff <- max(as.numeric(pv), 1e-300)
               abs(nes)*(-log10(pv_eff))
             } else 0
           }, numeric(1))),
           "max_abs_nes"  = max(abs(nes_mat[gid,]),na.rm=TRUE),
           "mean_abs_nes" = mean(abs(nes_mat[gid,]),na.rm=TRUE),
           "n_significant"= sum(!is.na(sig_mat[gid,])&sig_mat[gid,]=="*"),
           0
    )
  }
  
  selected_ids <- rownames(nes_mat)
  n_significant_anywhere <- sum(vapply(selected_ids, function(gid) {
    !is_not_significant_anywhere(gid)
  }, logical(1)))
  cat(sprintf("  Significant in >=1 condition: %d/%d\n",
              n_significant_anywhere, length(selected_ids)))

  # Apply interpretability before redundancy reduction.
  interpretation_filtered_ids <- character(0)
  if (filter_mode=="interpretation" && length(selected_ids)>0) {
    ids_before_interpretation <- selected_ids
    keep_ids <- selected_ids[vapply(selected_ids,is_interpretable_term,logical(1))]
    interpretation_filtered_ids <- setdiff(ids_before_interpretation, keep_ids)
    n_not_significant <- sum(vapply(interpretation_filtered_ids,
                                    is_not_significant_anywhere, logical(1)))
    n_unsupported <- length(interpretation_filtered_ids) - n_not_significant
    cat(sprintf("  Interpretation: %d retained | %d non-significant | %d unsupported\n",
                length(keep_ids), n_not_significant, n_unsupported))
    selected_ids <- keep_ids
    nes_mat <- nes_mat[selected_ids,,drop=FALSE]
    sig_mat <- sig_mat[selected_ids,,drop=FALSE]
    if (!is.null(cld_mat)) cld_mat <- cld_mat[selected_ids,,drop=FALSE]
  }

  elim_map <- list(); similarity_matrix_full <- NULL
  if (do_similarity_reduction && type_suffix=="GO_BP" && length(selected_ids)>1) {
    n <- length(selected_ids)
    similarity_matrix_full <- matrix(0,n,n,dimnames=list(selected_ids,selected_ids))
    diag(similarity_matrix_full) <- 1.0
    conditionwise_le_similarity <- function(term_a, term_b) {
      sims <- vapply(treatments, function(t) {
        le_a <- get_le(term_a, t)
        le_b <- get_le(term_b, t)
        if (!length(le_a) || !length(le_b)) return(NA_real_)
        compute_set_similarity(le_a, le_b, similarity_metric)
      }, numeric(1))
      shared <- is.finite(sims)
      if (!any(shared)) return(0)
      mean(sims[shared])
    }

    for (i in seq_len(n-1)) for (j in seq(i+1,n)) {
      a <- selected_ids[i]
      b <- selected_ids[j]
      v <- conditionwise_le_similarity(a, b)
      similarity_matrix_full[i,j] <- similarity_matrix_full[j,i] <- v
    }
    clusters <- cutree(
      hclust(as.dist(1-similarity_matrix_full), method="complete"),
      h = 1-similarity_threshold
    )
    winners <- unname(vapply(sort(unique(clusters)),function(cid) {
      candidates <- selected_ids[clusters==cid]
      cscores <- vapply(candidates,compute_score,numeric(1))
      candidates[which.max(cscores)]
    },character(1)))
    for (w in winners) {
      m <- selected_ids[clusters==clusters[w]]
      if(length(m)>1) {
        losers <- m[m!=w]
        elim_map[[w]] <- losers
      }
    }
    cat(sprintf("  Redundancy (%s, threshold %.3f): %d -> %d representatives (%d grouped)\n",
                metric_label(similarity_metric), similarity_threshold, length(selected_ids),
                length(winners), length(selected_ids) - length(winners)))
    selected_ids <- winners
    nes_mat <- nes_mat[selected_ids,,drop=FALSE]
    sig_mat <- sig_mat[selected_ids,,drop=FALSE]
    if (!is.null(cld_mat)) cld_mat <- cld_mat[selected_ids,,drop=FALSE]
  }
  
  n_terms <- length(selected_ids)
  if (n_terms==0) {
    cat("  No terms remaining\n")
    return(invisible(NULL))
  }
  eliminated_terms_df <- build_elimination_report()
  
  display_cld_mat_all <- cld_mat_all
  if (!is.null(cld_mat_all) && nrow(cld_mat_all) > 0) {
    chars_in_order <- function(x) {
      x <- x[nzchar(x)]
      if (!length(x)) return(character(0))
      z <- unlist(strsplit(x, ""), use.names = FALSE)
      z <- z[z %in% letters]
      z[!duplicated(z)]
    }
    for (gid in rownames(cld_mat_all)) {
      raw_row <- trimws(as.character(cld_mat_all[gid, ]))
      names(raw_row) <- treatments
      raw_row[is.na(raw_row)] <- ""
      sig_t <- treatments[!is.na(sig_mat_all[gid, treatments]) & sig_mat_all[gid, treatments] == "*"]
      sig_chars <- if (length(sig_t)) chars_in_order(raw_row[sig_t]) else character(0)
      all_chars <- chars_in_order(raw_row)
      ordered_chars <- c(sig_chars, setdiff(all_chars, sig_chars))
      if (!length(ordered_chars)) next
      mapping <- setNames(letters[seq_along(ordered_chars)], ordered_chars)
      for (t in treatments) {
        orig <- raw_row[t]
        if (!nzchar(orig)) next
        display_cld_mat_all[gid, t] <- paste(vapply(
          strsplit(orig, "")[[1]],
          function(ch) if (ch %in% names(mapping)) mapping[[ch]] else ch,
          character(1)
        ), collapse = "")
      }
    }
  }
  display_cld_mat <- if (!is.null(display_cld_mat_all))
    display_cld_mat_all[selected_ids, , drop = FALSE] else NULL
  
  hc_rows <- NULL
  if (n_terms>=2) {
    if (type_suffix=="GO_BP") {
      clust_result <- tryCatch({
        suppressPackageStartupMessages({library(GOSemSim);library(GO.db)})
        if (is.null(.anc_cache_global)) .anc_cache_global <<- build_anc_cache()

        sem_data <- NULL
        invisible(capture.output(
          sem_data <- GOSemSim::godata(
            annoDb    = orgdb,
            ont       = "BP",
            keytype   = key_type,
            computeIC = FALSE
          ),
          type = "output"
        ))
        sim <- GOSemSim::mgoSim(
          selected_ids,
          selected_ids,
          semData = sem_data,
          measure = "Wang",
          combine = NULL
        )
        sim <- as.matrix(sim)
        cl <- cluster_from_similarity(sim, selected_ids)
        if (length(cl$missing_from_similarity)>0) {
          cat(sprintf(
            "  Semantic similarity unavailable for %d selected GO term(s); retained with zero similarity for row ordering.\n",
            length(cl$missing_from_similarity)
          ))
        }
        cl
      },error=function(e){cat(sprintf("  Sim failed (%s) — using correlation\n",conditionMessage(e)));NULL})
      if (!is.null(clust_result)) {
        hc_rows <- clust_result$hc; sel <- clust_result$sel
        selected_ids <- sel
        nes_mat <- nes_mat[sel,,drop=FALSE]
        sig_mat <- sig_mat[sel,,drop=FALSE]
        if (!is.null(cld_mat)) cld_mat <- cld_mat[sel,,drop=FALSE]
        n_terms <- length(sel)
      }
    }
    if (is.null(hc_rows) && type_suffix == "KEGG") {
      clust_result <- tryCatch({
        kegg_t2g_candidates <- c(
          file.path("kegg_cache", paste0(kegg_org, "_term2gene.tsv"))
        )
        kegg_t2g_file <- kegg_t2g_candidates[file.exists(kegg_t2g_candidates)][1]

        pathway_genes <- setNames(vector("list", length(selected_ids)), selected_ids)

        if (!is.na(kegg_t2g_file) && length(kegg_t2g_file) == 1L) {
          kt <- read.table(kegg_t2g_file, header=TRUE, sep="\t",
                           stringsAsFactors=FALSE, check.names=FALSE,
                           quote="", comment.char="")
          term_col <- intersect(c("TERM", "term", "pathway", "ID"), colnames(kt))[1]
          gene_col <- intersect(c("GENE", "gene", "ORF"), colnames(kt))[1]
          if (is.na(term_col) || is.na(gene_col)) stop("KEGG cache columns TERM/GENE not found")
          for (gid in selected_ids) {
            pathway_genes[[gid]] <- unique(as.character(kt[kt[[term_col]] == gid, gene_col]))
            pathway_genes[[gid]] <- pathway_genes[[gid]][nzchar(pathway_genes[[gid]])]
          }
        } else {
          for (gid in selected_ids) {
            pathway_genes[[gid]] <- unique(unlist(lapply(treatments, function(t) get_le(gid, t))))
          }
        }

        n <- length(selected_ids)
        sim <- matrix(0, n, n, dimnames=list(selected_ids, selected_ids))
        diag(sim) <- 1
        if (n > 1L) {
          for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
            v <- compute_set_similarity(pathway_genes[[selected_ids[i]]],
                                        pathway_genes[[selected_ids[j]]],
                                        similarity_metric)
            sim[i, j] <- sim[j, i] <- if (is.na(v)) 0 else v
          }
        }
        cluster_from_similarity(sim, selected_ids)
      }, error=function(e) {
        cat(sprintf("  KEGG similarity clustering failed (%s) — using NES correlation\n",
                    conditionMessage(e)))
        NULL
      })
      if (!is.null(clust_result)) {
        hc_rows <- clust_result$hc
      }
    }

    if (is.null(hc_rows)) {
      nc <- nes_mat; nc[is.na(nc)] <- 0
      if (ncol(nc)>=2) {
        cr <- tryCatch(cor(t(nc),method="pearson"),error=function(e) NULL)
        if (!is.null(cr)) { cr[is.nan(cr)] <- 0; hc_rows <- hclust(as.dist(1-cr),method="average") }
      }
    }
  }
  if (!is.null(hc_rows)&&length(hc_rows$labels)!=nrow(nes_mat)) { cat("  WARNING: hc mismatch\n"); hc_rows <- NULL }
  
  if (!is.null(hc_rows)) {
    go_order <- hc_rows$labels[hc_rows$order]
    go_order <- c(go_order[go_order %in% selected_ids],setdiff(selected_ids,go_order))
  } else go_order <- selected_ids
  
  de_data  <- list()
  de_files <- list.files(
    results_dir,
    pattern = "\\.(txt|tsv|csv)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  de_files <- de_files[file.info(de_files)$isdir %in% FALSE]
  
  read_de_file <- function(path) {
    ext <- tolower(tools::file_ext(path))
    readers <- if (ext == "csv") {
      list(
        function() read.csv(path, header = TRUE, stringsAsFactors = FALSE,
                            check.names = FALSE),
        function() read.table(path, header = TRUE, sep = "\t", quote = "",
                              comment.char = "", stringsAsFactors = FALSE,
                              check.names = FALSE)
      )
    } else {
      list(
        function() read.table(path, header = TRUE, sep = "\t", quote = "",
                              comment.char = "", stringsAsFactors = FALSE,
                              check.names = FALSE),
        function() read.table(path, header = TRUE, sep = "", quote = "",
                              comment.char = "", stringsAsFactors = FALSE,
                              check.names = FALSE),
        function() read.csv(path, header = TRUE, stringsAsFactors = FALSE,
                            check.names = FALSE)
      )
    }
    for (fn in readers) {
      x <- tryCatch(fn(), error = function(e) NULL)
      if (!is.null(x) && nrow(x) > 0L && ncol(x) >= 2L) return(x)
    }
    NULL
  }
  
  norm_contrast_name <- function(x) {
    x <- tools::file_path_sans_ext(basename(x))
    x <- tolower(trimws(x))
    x <- gsub("[[:space:].]+", "_", x)
    x <- gsub("-+", "-", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  split_contrast_name <- function(path) {
    nm <- norm_contrast_name(path)
    pos <- regexpr("_vs_", nm, fixed = TRUE)
    if (pos[1] < 1L) return(NULL)
    lhs <- substr(nm, 1L, pos[1] - 1L)
    rhs <- substr(nm, pos[1] + attr(pos, "match.length"), nchar(nm))
    if (!nzchar(lhs) || !nzchar(rhs)) return(NULL)
    list(lhs = lhs, rhs = rhs)
  }
  
  norm_level <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[[:space:].]+", "_", x)
    x <- gsub("_+", "_", x)
    x
  }
  
  de_ref_level <- if (!is.null(reference_level) && nzchar(reference_level)) {
    reference_level
  } else {
    ref_level
  }
  
  contrast_info <- lapply(de_files, split_contrast_name)
  valid_contrast <- !vapply(contrast_info, is.null, logical(1))
  de_files_valid <- de_files[valid_contrast]
  contrast_info  <- contrast_info[valid_contrast]
  
  de_ref_display <- if (interaction_mode)
    split_condition_components(de_ref_level)[1]
  else de_ref_level
  
  
  for (t in treatments) {
    target_t <- norm_level(t)
    
    if (interaction_mode) {
      cf_ref <- split_condition_components(de_ref_level)[1]
      t_parts   <- split_condition_components(t)
      ref_parts <- t_parts
      ref_parts[1] <- cf_ref
      target_ref <- norm_level(paste(ref_parts, collapse="×"))
    } else {
      target_ref <- norm_level(de_ref_level)
    }
    
    matches <- which(vapply(contrast_info, function(z) {
      identical(norm_level(z$lhs), target_t) &&
        identical(norm_level(z$rhs), target_ref)
    }, logical(1)))
    
    if (!length(matches)) {
      cat(sprintf(
        "  DE data MISSING: no %s_vs_%s table found in %s\n",
        t, if(interaction_mode) paste(split_condition_components(t) |> (\(p){p[1]<-cf_ref;p})(), collapse="×") else de_ref_level,
        results_dir
      ))
      next
    }
    
    candidate_files <- de_files_valid[matches]
    if (length(candidate_files) > 1L) {
      mt <- file.info(candidate_files)$mtime
      candidate_files <- candidate_files[order(mt, decreasing = TRUE)]
      cat(sprintf(
        "  DE data: multiple matches for %s; using newest: %s\n",
        t, basename(candidate_files[1])
      ))
    }
    
    chosen <- candidate_files[1]
    df <- read_de_file(chosen)
    if (is.null(df)) {
      cat(sprintf("  DE data FAILED: %s -> %s\n", t, basename(chosen)))
      next
    }
    
    de_data[[t]] <- df
  }
  
  all_orfs <- unique(unlist(lapply(selected_ids,function(g)
    unlist(lapply(treatments,function(t) get_le(g,t))))))
  all_orfs <- trimws(as.character(all_orfs))
  all_orfs <- unique(all_orfs[!is.na(all_orfs) & nzchar(all_orfs)])

  is_valid_symbol <- function(nm) {
    nm <- trimws(as.character(nm))
    if (length(nm) != 1L || is.na(nm) || !nzchar(nm)) return(FALSE)

    if (is_sgd_yeast) {
      return(nchar(nm) <= 24L && grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", nm))
    }

    nchar(nm) <= 40L && !grepl("[,;()]", nm)
  }

  normalize_annotation_value <- function(x) {
    x <- enc2utf8(as.character(x)); x[is.na(x)] <- ""
    x <- tolower(trimws(x))
    x <- gsub("α", "alpha", x, fixed=TRUE)
    x <- gsub("β", "beta",  x, fixed=TRUE)
    x <- gsub("γ", "gamma", x, fixed=TRUE)
    x <- gsub("[^[:alnum:]]+", " ", x)
    trimws(gsub("[[:space:]]+", " ", x))
  }

  target_keys <- tryCatch(
    AnnotationDbi::keys(orgdb_obj, keytype=key_type),
    error=function(e) character(0)
  )
  # SGD ORF loci may have no ENSEMBL entry. Keep them as report-only
  # identifiers rather than discarding valid leading-edge genes.
  sgd_orf_keys <- if (is_sgd_yeast && "ORF" %in% available_keytypes) {
    tryCatch(AnnotationDbi::keys(orgdb_obj, keytype="ORF"),
             error=function(e) character(0))
  } else character(0)
  report_keys <- unique(c(target_keys, sgd_orf_keys))

  .annotation_reverse_lookup <- NULL
  get_annotation_reverse_lookup <- function() {
    if (!is.null(.annotation_reverse_lookup)) return(.annotation_reverse_lookup)

    base_keytype <- if ("ENSEMBL" %in% available_keytypes) "ENSEMBL" else key_type
    base_keys <- tryCatch(
      AnnotationDbi::keys(orgdb_obj, keytype=base_keytype),
      error=function(e) character(0)
    )
    if (!length(base_keys)) {
      .annotation_reverse_lookup <<- setNames(character(0), character(0))
      return(.annotation_reverse_lookup)
    }

    source_columns <- intersect(
      c("GENENAME", "COMMON", "SYMBOL", "ORF", "ENSEMBL", "SGD",
        "ALIAS", "DESCRIPTION", "REFSEQ", "ENTREZID"),
      available_columns
    )
    requested_columns <- unique(c(key_type, source_columns))
    requested_columns <- setdiff(requested_columns, base_keytype)

    ann <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        orgdb_obj,
        keys=base_keys,
        columns=requested_columns,
        keytype=base_keytype
      )),
      error=function(e) NULL
    )
    if (is.null(ann) || !nrow(ann)) {
      .annotation_reverse_lookup <<- setNames(character(0), character(0))
      return(.annotation_reverse_lookup)
    }

    target_col <- if (key_type %in% colnames(ann)) key_type else base_keytype
    source_columns <- intersect(unique(c(base_keytype, source_columns)), colnames(ann))

    pairs <- lapply(source_columns, function(col) {
      data.frame(
        annotation=normalize_annotation_value(ann[[col]]),
        target=as.character(ann[[target_col]]),
        stringsAsFactors=FALSE
      )
    })
    pairs <- do.call(rbind, pairs)
    pairs <- pairs[nzchar(pairs$annotation) & !is.na(pairs$target) & nzchar(pairs$target), , drop=FALSE]
    pairs <- unique(pairs)

    by_annotation <- split(pairs$target, pairs$annotation)
    unambiguous <- vapply(by_annotation, function(z) {
      z <- unique(z[!is.na(z) & nzchar(z)])
      if (length(z) == 1L) z else NA_character_
    }, character(1))
    unambiguous <- unambiguous[!is.na(unambiguous) & nzchar(unambiguous)]

    .annotation_reverse_lookup <<- unambiguous
    .annotation_reverse_lookup
  }

  to_internal_id <- function(ids) {
    ids <- trimws(as.character(ids)); ids[is.na(ids)] <- ""
    out <- ids
    if (toupper(key_type) == "ENSEMBL") out <- sub("\\.[0-9]+$", "", out, perl=TRUE)
    resolved <- nzchar(out) & out %in% target_keys
    unresolved <- which(nzchar(out) & !resolved)
    if (!length(unresolved)) return(out)

    vals <- out[unresolved]
    candidates <- setdiff(unique(c(symbol_column,"ENTREZID","ENSEMBL","REFSEQ",
                                    "SYMBOL","COMMON","ORF","SGD")), key_type)
    candidates <- candidates[candidates %in% available_keytypes]
    mapped_final <- rep(NA_character_, length(vals))
    for (kt in candidates) {
      need <- is.na(mapped_final) | !nzchar(mapped_final)
      if (!any(need)) break
      keys_kt <- tryCatch(AnnotationDbi::keys(orgdb_obj,keytype=kt),
                          error=function(e) character(0))
      hit_vals <- unique(vals[need & vals %in% keys_kt])
      if (!length(hit_vals)) next
      mp <- tryCatch(suppressMessages(AnnotationDbi::mapIds(
        orgdb_obj, keys=hit_vals, column=key_type, keytype=kt,
        multiVals="CharacterList"
      )), error=function(e) NULL)
      if (is.null(mp)) next
      mp_unique <- vapply(hit_vals, function(vv) {
        z <- unique(as.character(mp[[vv]]))
        z <- z[!is.na(z) & nzchar(z)]
        if (length(z) == 1L) z else NA_character_
      }, character(1))
      names(mp_unique) <- hit_vals
      idx <- which(need & vals %in% names(mp_unique))
      mapped_final[idx] <- unname(mp_unique[vals[idx]])
    }

    need <- is.na(mapped_final) | !nzchar(mapped_final)
    if (any(need)) {
      rev_lu <- get_annotation_reverse_lookup()
      if (length(rev_lu)) {
        normalized_vals <- normalize_annotation_value(vals[need])
        reverse_mapped <- unname(rev_lu[normalized_vals])
        idx_need <- which(need)
        good_reverse <- !is.na(reverse_mapped) & nzchar(reverse_mapped)
        if (any(good_reverse)) {
          mapped_final[idx_need[good_reverse]] <- reverse_mapped[good_reverse]
        }
      }
    }

    good <- !is.na(mapped_final) & nzchar(mapped_final)
    out[unresolved[good]] <- mapped_final[good]
    # Preserve an ORF if it is valid in SGD but has no ENSEMBL mapping.
    # This fallback does not assign differential-expression statistics.
    out
  }

  all_orfs <- unique(to_internal_id(all_orfs))
  invalid_all_orfs <- all_orfs[nzchar(all_orfs) & !all_orfs %in% report_keys]
  if (length(invalid_all_orfs)) {
    cat(sprintf(
      "  WARNING: %d unresolved leading-edge labels were excluded (examples: %s)\n",
      length(unique(invalid_all_orfs)),
      paste(head(unique(invalid_all_orfs), 5L), collapse=", ")
    ))
  }
  all_orfs <- all_orfs[nzchar(all_orfs) & all_orfs %in% report_keys]

  display_name_column <- if (is_sgd_yeast && "GENENAME" %in% available_columns) {
    "GENENAME"
  } else {
    symbol_column
  }

  .ORF2NAME <- setNames(all_orfs, all_orfs)
  if (isTRUE(map_gene_labels)) {
    # Lookup the ordinary analysis IDs and SGD-only ORFs separately.
    for (kt in unique(c(key_type, if (length(sgd_orf_keys)) "ORF"))) {
      ids <- if (identical(kt, key_type)) {
        intersect(all_orfs, target_keys)
      } else {
        setdiff(intersect(all_orfs, sgd_orf_keys), target_keys)
      }
      if (!length(ids)) next
      nms <- tryCatch(suppressMessages(AnnotationDbi::mapIds(
        orgdb_obj, keys=ids, column=display_name_column,
        keytype=kt, multiVals="CharacterList"
      )), error=function(e) NULL)
      if (is.null(nms)) next
      for (id in ids) {
        vals <- unique(as.character(nms[[id]]))
        vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
        if (length(vals) == 1L && is_valid_symbol(vals[[1]])) .ORF2NAME[[id]] <- vals[[1]]
      }
    }
  }

  disp <- function(id) {
    nm <- .ORF2NAME[id]
    if(length(nm)>0&&!is.na(nm)&&nchar(nm)>0) unname(nm) else id
  }

  to_orf <- to_internal_id

  get_le_raw <- get_le
  get_le <- function(gid, t) {
    ids <- get_le_raw(gid, t)
    if (!length(ids)) return(character(0))
    ids <- to_internal_id(ids)
    ids <- ids[!is.na(ids) & nzchar(ids) & ids %in% report_keys]
    unique(ids)
  }

  detect_de_columns <- function(df) {
    nms <- colnames(df)
    gene_cn <- intersect(c("gene","Gene","GENE","ORF","orf","systematic_name",
                           "Systematic.name","ensembl","ENSEMBL","sgdid","SGDID"), nms)[1]
    lfc_cn  <- intersect(c("logFC","log2FoldChange","log2FC","LFC","log2.fold.change"), nms)[1]
    padj_cn <- intersect(c("adj.P.Val","padj","FDR","qvalue","q.value","adj_pvalue",
                           "adjusted_p_value"), nms)[1]
    if (is.na(gene_cn)) {
      cand <- nms[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
      if (length(cand)) gene_cn <- cand[1]
    }
    list(gene=gene_cn, lfc=lfc_cn, padj=padj_cn)
  }
  
  get_de_vectors <- function(df) {
    cc <- detect_de_columns(df)
    if (is.na(cc$lfc) || is.na(cc$padj)) return(NULL)
    input_ids <- if (!is.na(cc$gene)) as.character(df[[cc$gene]]) else rownames(df)
    if (is.null(input_ids)) return(NULL)
    annotation_ids <- if ("original_gene" %in% colnames(df)) as.character(df$original_gene) else input_ids
    orfs <- to_orf(annotation_ids)
    lfc  <- suppressWarnings(as.numeric(df[[cc$lfc]]))
    padj <- suppressWarnings(as.numeric(df[[cc$padj]]))
    ok <- nzchar(orfs) & is.finite(lfc) & is.finite(padj)
    if (!any(ok)) return(NULL)
    out <- data.frame(
      orf=orfs[ok], input_gene=input_ids[ok], original_gene=annotation_ids[ok],
      lfc=lfc[ok], padj=padj[ok], stringsAsFactors=FALSE
    )
    out <- out[order(out$padj, -abs(out$lfc)), , drop=FALSE]
    rownames(out) <- NULL
    out
  }
  
  de_orfs <- unique(unlist(lapply(de_data, function(df) {
    dv <- get_de_vectors(df)
    if (is.null(dv)) character(0) else dv$orf
  })))
  missing_name_orfs <- setdiff(de_orfs, names(.ORF2NAME))
  if (length(missing_name_orfs) && isTRUE(map_gene_labels)) {
    extra_names_raw <- tryCatch({
      suppressMessages(AnnotationDbi::mapIds(
        orgdb_obj, keys=missing_name_orfs, column=display_name_column,
        keytype=key_type, multiVals="CharacterList"
      ))
    }, error=function(e) NULL)
    extra_names <- setNames(vapply(missing_name_orfs, function(o) {
      if (is.null(extra_names_raw)) return(o)
      vals <- unique(as.character(extra_names_raw[[o]]))
      vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
      if (length(vals) == 1L && is_valid_symbol(vals[[1]])) vals[[1]] else o
    }, character(1)), missing_name_orfs)
    .ORF2NAME <- c(.ORF2NAME, extra_names)
  }
  
  if (enable_sgd_pathway_lookup && isTRUE(.SGD_CACHE_AVAILABLE)) {
    # Use all valid SGD ORF mappings while counting the report gene once.
    sgd_ids <- if (identical(key_type, "ENSEMBL")) {
      setNames(lapply(all_orfs, function(x) x), all_orfs)
    } else {
      mp <- tryCatch(suppressMessages(AnnotationDbi::mapIds(
        orgdb_obj, keys=all_orfs, column="ENSEMBL",
        keytype=key_type, multiVals="CharacterList"
      )), error=function(e) NULL)
      setNames(lapply(all_orfs, function(id) {
        if (is.null(mp)) return(character(0))
        vals <- unique(as.character(mp[[id]]))
        vals[!is.na(vals) & nzchar(vals)]
      }), all_orfs)
    }

    lookup_ids <- unique(unlist(sgd_ids, use.names=FALSE))
    lookup_ids <- lookup_ids[!is.na(lookup_ids) & nzchar(lookup_ids)]
    n_in_cache <- sum(vapply(lookup_ids, function(o) !is.null(.SGD_CACHE[[o]]), logical(1)))
    if (n_in_cache < length(lookup_ids))
      cat(sprintf("  SGD pathway cache coverage: %d/%d mapped systematic IDs\n", n_in_cache, length(lookup_ids)))

    cached_by_orf <- setNames(lapply(lookup_ids, fetch_sgd_pathways), lookup_ids)
    orf_pathways <- setNames(lapply(all_orfs, function(id) {
      sids <- sgd_ids[[id]]
      if (is.null(sids) || !length(sids)) return(list())
      merged <- list()
      for (sid in sids) {
        pwlist <- cached_by_orf[[sid]]
        if (is.null(pwlist) || !length(pwlist)) next
        for (pid in names(pwlist)) if (is.null(merged[[pid]])) merged[[pid]] <- pwlist[[pid]]
      }
      merged
    }), all_orfs)

    status_tab <- data.frame(
      SGD_ID=lookup_ids,
      status=vapply(lookup_ids, function(o) {
        st <- .SGD_STATUS[[o]]
        if (is.null(st)) "not_in_cache" else as.character(st)
      }, character(1)),
      message="",
      stringsAsFactors=FALSE
    )
    write.table(status_tab, file.path(output_dir, "SGD_pathway_lookup_status.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  } else {
    orf_pathways <- setNames(rep(list(list()), length(all_orfs)), all_orfs)
  }

  make_heatmap_table <- function(ids, include_grouping = FALSE) {
    ids <- unique(ids[ids %in% rownames(nes_mat_all)])
    out <- data.frame(
      ID = ids,
      Description = unname(desc_lu[ids]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (t in treatments) {
      out[[t]] <- vapply(ids, function(g) {
        cld <- if (!is.null(cld_mat_all) && g %in% rownames(cld_mat_all)) display_cld_mat_all[g, t] else ""
        format_heatmap_value(nes_mat_all[g, t], sig_mat_all[g, t], cld)
      }, character(1))
    }
    if (interaction_mode && length(all_int_cols)) {
      for (col in all_int_cols) {
        excel_label <- if (exists("int_col_labels_full", inherits=FALSE) &&
                           col %in% names(int_col_labels_full))
          unname(int_col_labels_full[col]) else unname(int_col_labels[col])
        out[[excel_label]] <- vapply(ids, function(g) {
          v <- int_sig_vecs[[col]][g]
          if (is.null(v) || is.na(v)) "" else as.character(v)
        }, character(1))
      }
    } else {
      out[[difference_column_label]] <- unname(diff_col_vec[ids])
    }
    if (include_grouping) {
      representative <- setNames(ids, ids)
      status <- setNames(rep("Retained", length(ids)), ids)
      for (w in names(elim_map)) {
        losers <- elim_map[[w]]
        losers <- losers[losers %in% ids]
        if (length(losers)) {
          representative[losers] <- w
          status[losers] <- paste0("Grouped as redundant by ", metric_label(similarity_metric))
        }
      }
      out[["Grouping status"]] <- unname(status[ids])
      out[["Representative ID"]] <- unname(representative[ids])
    }
    out
  }
  
  de_vec_cache <- lapply(de_data, get_de_vectors)
  
  de_lookup <- lapply(de_vec_cache, function(dv) {
    if (is.null(dv)) return(NULL)
    by_orf <- split(seq_len(nrow(dv)), dv$orf)
    one <- vapply(by_orf, length, integer(1)) == 1L
    ids <- names(by_orf)
    lfc <- setNames(rep(NA_real_, length(ids)), ids)
    padj <- setNames(rep(NA_real_, length(ids)), ids)
    for (id in ids[one]) {
      ii <- by_orf[[id]][1]
      lfc[id] <- dv$lfc[ii]
      padj[id] <- dv$padj[ii]
    }
    list(lfc=lfc, padj=padj, n=setNames(lengths(by_orf), ids))
  })

  for (t in treatments) {
    sig_terms_t <- selected_ids[!is.na(sig_mat[selected_ids, t]) & sig_mat[selected_ids, t] == "*"]
    le_t <- unique(unlist(lapply(sig_terms_t, function(g) get_le(g, t)), use.names = FALSE))
    lk <- de_lookup[[t]]
    if (is.null(lk)) {
      cat(sprintf("  Volcano check: %s — DE table not loaded\n", t))
    } else if (!length(le_t)) {
      cat(sprintf("  Volcano check: %s — no leading-edge genes in retained significant terms\n", t))
    } else {
      matched <- intersect(le_t, names(lk$lfc))
      if (length(matched) < length(le_t)) {
        cat(sprintf("  WARNING: volcano mapping for %s matched %d/%d leading-edge genes (%.1f%%)\n",
                    t, length(matched), length(le_t), 100 * length(matched) / length(le_t)))
        if (!length(matched)) {
          cat(sprintf("    Leading-edge examples: %s\n", paste(head(le_t, 5), collapse = ", ")))
          cat(sprintf("    DE-ID examples       : %s\n", paste(head(names(lk$lfc), 5), collapse = ", ")))
        }
      }
    }
  }
  
  leading_rows <- list(); lr_i <- 0L
  for (g in selected_ids) for (t in treatments) {
    term_sig <- !is.na(sig_mat_all[g, t]) && sig_mat_all[g, t] == "*"
    cld <- if (!is.null(cld_mat_all)) display_cld_mat_all[g, t] else ""
    genes <- if (term_sig) get_le(g, t) else character(0)
    if (!length(genes)) {
      lr_i <- lr_i + 1L
      leading_rows[[lr_i]] <- data.frame(
        Term_ID = g, Description = unname(desc_lu[g]), Treatment = t,
        Enrichment = if (term_sig) "Significant; empty leading edge" else "ns",
        NES = if (term_sig) nes_mat_all[g, t] else NA_real_, CLD = cld,
        Gene_ORF = NA_character_, Gene_name = NA_character_,
        Gene_padj = NA_real_, Gene_log2FC = NA_real_, stringsAsFactors = FALSE
      )
    } else {
      dv <- de_vec_cache[[t]]
      for (orf in genes) {
        idx_all <- if (!is.null(dv)) which(dv$orf == orf) else integer(0)
        if (!length(idx_all)) idx_all <- NA_integer_
        for (idx in idx_all) {
          lr_i <- lr_i + 1L
          leading_rows[[lr_i]] <- data.frame(
            Term_ID = g, Description = unname(desc_lu[g]), Treatment = t,
            Enrichment = "Significant", NES = nes_mat_all[g, t], CLD = cld,
            Gene_ORF = orf,
            Gene_name = if (!is.na(idx) && length(idx) && !is.na(dv$input_gene[idx])) dv$input_gene[idx] else disp(orf),
            Gene_padj = if (!is.na(idx)) dv$padj[idx] else NA_real_,
            Gene_log2FC = if (!is.na(idx)) dv$lfc[idx] else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  leading_df <- if (length(leading_rows)) do.call(rbind, leading_rows) else data.frame()

  make_interpretation_filtered_table <- function(ids) {
    ids <- unique(ids[ids %in% rownames(nes_mat_all)])
    out <- data.frame(
      ID = ids,
      Description = unname(desc_lu[ids]),
      Filtering_reason = vapply(ids, function(gid) {
        if (is_not_significant_anywhere(gid)) {
          "Not significant in any condition"
        } else {
          "Pattern lacks required pairwise/interaction support"
        }
      }, character(1)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (t in treatments) {
      out[[paste0(t, "_NES")]] <- as.numeric(nes_mat_all[ids, t])
      out[[paste0(t, "_significance")]] <- as.character(sig_mat_all[ids, t])
      if (!is.null(cld_mat_all) && t %in% colnames(cld_mat_all)) {
        out[[paste0(t, "_CLD")]] <- as.character(display_cld_mat_all[ids, t])
      }
    }
    if (length(ids)) {
      out[["NES_pattern"]] <- apply(nes_mat_all[ids, treatments, drop=FALSE], 1, function(z) {
        paste(vapply(z, function(v) if (is.finite(v)) sprintf("%.3f", v) else "NA", character(1)), collapse=" | ")
      })
      out[["Significance_pattern"]] <- apply(sig_mat_all[ids, treatments, drop=FALSE], 1, function(z) {
        z[is.na(z) | !nzchar(trimws(z))] <- "NA"
        paste(z, collapse=" | ")
      })
    } else {
      out[["NES_pattern"]] <- character(0)
      out[["Significance_pattern"]] <- character(0)
    }
    out
  }
  interpretation_filtered_df <- make_interpretation_filtered_table(interpretation_filtered_ids)
  
  all_valid_ids <- rownames(nes_mat_all)
  all_terms_df <- make_heatmap_table(all_valid_ids, FALSE)

  representative <- setNames(all_valid_ids, all_valid_ids)
  grouping_status <- setNames(rep("Not grouped (singleton)", length(all_valid_ids)), all_valid_ids)
  for (w in names(elim_map)) {
    losers <- intersect(elim_map[[w]], all_valid_ids)
    if (w %in% all_valid_ids) grouping_status[w] <- if (length(losers)) "Jaccard representative retained" else grouping_status[w]
    if (length(losers)) {
      representative[losers] <- w
      grouping_status[losers] <- paste0("Merged by ", metric_label(similarity_metric), " grouping")
    }
  }
  all_terms_df[["Jaccard grouping status"]] <- unname(grouping_status[all_terms_df$ID])
  all_terms_df[["Jaccard representative ID"]] <- unname(representative[all_terms_df$ID])
  all_terms_df[["Present in final heatmap"]] <- ifelse(all_terms_df$ID %in% selected_ids, "yes", "no")

  xlsx_sheets <- list(
    Heatmap = make_heatmap_table(selected_ids, FALSE),
    Leading_edge = leading_df
  )
  if (type_suffix == "GO_BP") xlsx_sheets[["All_terms"]] <- all_terms_df
  xlsx_sheets[["Interpretation_filtered_out"]] <- interpretation_filtered_df
  xlsx_sheets[["Eliminated_terms"]] <- eliminated_terms_df

  elim_files_to_remove <- c(
    file.path(output_dir, paste0(type_suffix, "_eliminated_terms.tsv")),
    file.path(output_dir, paste0(type_suffix, "_eliminated_terms.xlsx")),
    file.path(output_dir, paste0(type_suffix, "_Eliminated_terms.xlsx"))
  )
  elim_files_to_remove <- elim_files_to_remove[file.exists(elim_files_to_remove)]
  if (length(elim_files_to_remove)) unlink(elim_files_to_remove)

  xlsx_path <- file.path(output_dir, paste0(type_suffix, "_heatmap_report.xlsx"))
  workbook_written <- write_analysis_workbook(xlsx_path, xlsx_sheets)
  
  sig_abs_vals <- abs(nes_mat[sig_mat == "*" & is.finite(nes_mat)])
  robust_abs <- if (length(sig_abs_vals)) as.numeric(stats::quantile(sig_abs_vals, 0.95, na.rm=TRUE, names=FALSE)) else NA_real_
  max_abs <- if(!is.null(max_nes_scale)) max_nes_scale else max(1.5, robust_abs, na.rm=TRUE)
  if (!is.finite(max_abs)||max_abs==0) max_abs <- 3
  col_fn <- colorRamp2(c(-max_abs, -0.55*max_abs, 0, 0.55*max_abs, max_abs),
                       c("#3F6F9F", "#86ACCF", "#F7F7F7", "#E49482", "#B34F55"))
  
  nes_display_mat <- nes_mat
  for (g in rownames(nes_display_mat)) for (t in treatments) {
    s <- sig_mat[g,t]; if(!is.na(s)&&s=="ns") nes_display_mat[g,t] <- NA
  }
  
  row_labels_vec <- vapply(selected_ids,function(g) {
    d <- if(g %in% names(desc_lu)) desc_lu[g] else g; d
  },character(1))
  diff_ann <- diff_col_vec[selected_ids]; diff_ann[is.na(diff_ann)] <- ""
  
  lf_img <- if (use_cld&&!is.null(cld_mat)) {
    function(j,i,x,y,w,h,fill) {
      gid <- selected_ids[i]; t <- treatments[j]
      sig <- sig_mat[cbind(gid,t)]; sig[is.na(sig)] <- ""
      ltr <- display_cld_mat[cbind(gid,t)]; ltr[is.na(ltr)] <- ""; ltr <- trimws(ltr)
      ok  <- nchar(ltr) > 0
      if (any(ok)) grid.text(ltr[ok], x[ok], y[ok],
                             gp = gpar(fontsize = 7, col = "#222222"))
    }
  } else NULL
  
  longest_term_label <- row_labels_vec[which.max(nchar(row_labels_vec))]
  measured_label_mm <- tryCatch(
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(longest_term_label, gp=gpar(fontsize=7))),
      "mm", valueOnly=TRUE
    ),
    error=function(e) max(nchar(row_labels_vec), na.rm=TRUE) * 1.72
  )
  term_label_width <- max(30, measured_label_mm + 7)
  
  if (interaction_3way && n_int_cols == 0L && length(int_factors_from_manifest) >= 2L) {
    pfx <- tolower(int_label_prefix)
    background_order <- setdiff(header_factor_order, contrast_factor)
    facs <- background_order[background_order %in% int_factors_from_manifest]
    if (length(facs) < 2L)
      facs <- unique(c(facs, int_factors_from_manifest))[seq_len(2L)]
    all_int_cols <- c(paste0("interaction_x_", facs), "interaction_3way")
    int_col_labels <- setNames(
      c(paste0(tolower(contrast_factor), " × ", facs),
        paste0(tolower(contrast_factor), " × ", paste(facs, collapse=" × "))),
      all_int_cols
    )
    int_col_labels_full <- setNames(
      c(paste0(pfx, " × ", facs), paste0(pfx, " × ", paste(facs, collapse=" × "))),
      all_int_cols
    )
    int_sig_vecs <- setNames(lapply(all_int_cols, function(z)
      setNames(rep(NA_character_, nrow(summ)), summ$ID)), all_int_cols)
    n_int_cols <- length(all_int_cols)
  }
  
  int_col_palette <- c("#4DABF7","#FF8787","#63E6BE","#FFA94D","#9775FA","#F783AC")
  int_col_colors  <- c(
    setNames(rep(int_col_palette, length.out=length(int_x_cols)), int_x_cols),
    if (has_3way_col) c(interaction_3way="#888888") else character(0)
  )
  
  ann_args <- list()
  
  show_static_difference_col <- !show_formal_int_cols && (interaction_mode || n_treat != 2L)
  if (show_static_difference_col) {
    ann_args[["*"]] <- anno_text(
      diff_ann, which="row",
      gp=gpar(fontsize=8, col="#C92A2A", fontface="bold"),
      just="center", width=unit(6,"mm")
    )
  }
  
  if (show_formal_int_cols && n_int_cols > 0) {
    three_way_sig <- if (has_3way_col)
      vapply(selected_ids, function(g) {
        any(vapply(int_3way_cols, function(col) {
          v <- int_sig_vecs[[col]][g]
          !is.null(v) && !is.na(v) && v == "*"
        }, logical(1)))
      }, logical(1)) else rep(FALSE, length(selected_ids))
    
    for (icol in all_int_cols) {
      col_label <- int_col_labels[icol]
      is_3way   <- icol %in% int_3way_cols
      sig_vec   <- vapply(selected_ids, function(g) {
        v <- int_sig_vecs[[icol]][g]
        !is.null(v) && !is.na(v) && v == "*"
      }, logical(1))
      display_vec <- if (!is_3way) ifelse(three_way_sig, "", ifelse(sig_vec, "*", "")) else
        ifelse(sig_vec, "*", "")
      txt_col <- ifelse(display_vec == "*", "#C92A2A", "#DDDDDD")
      ann_args[[col_label]] <- anno_simple(
        rep("white", length(display_vec)), which="row",
        col = c(white = "white"), border = FALSE,
        gp = gpar(col="#D9D9D9", lwd=0.35),
        pch = ifelse(display_vec == "*", "*", ""),
        pt_gp = gpar(col="#C92A2A", fontsize=9, fontface="bold"),
        pt_size = unit(3.2, "mm"),
        width = unit(9, "mm")
      )
    }
  }
  
  ann_args[["Term"]] <- anno_text(
    row_labels_vec, which="row",
    gp=gpar(fontsize=7, col="#111111"),
    just="left", location=0,
    width=unit(term_label_width,"mm")
  )
  n_visual_ann <- length(ann_args)
  ann_args[["gap"]] <- if (n_visual_ann > 1L) {
    unit(c(rep(0, max(0L, n_visual_ann - 2L)), 2.5), "mm")
  } else unit(0.7, "mm")
  ann_args[["annotation_name_gp"]] <- gpar(fontsize=7.5, fontface="bold", col="#222222")
  ann_args[["annotation_name_side"]] <- "top"
  ann_args[["annotation_name_offset"]] <- unit(2.5, "mm")
  ann_args[["annotation_name_rot"]] <- 90
  ann_args[["show_annotation_name"]] <- show_formal_int_cols && n_int_cols > 0
  
  right_ann <- do.call(rowAnnotation, ann_args)
  
  factor_palette_bank <- list(
    c("#0072B2", "#E69F00", "#56B4E9", "#D55E00", "#332288", "#999933"),
    c("#009E73", "#D55E00", "#44AA99", "#CCBB44", "#117733", "#882255"),
    c("#332288", "#DDCC77", "#661100", "#88CCEE", "#999933", "#4477AA"),
    c("#4477AA", "#EE7733", "#228833", "#CCBB44", "#66CCEE", "#7755AA")
  )

  background_factors <- setdiff(header_factor_order, contrast_factor)
  preferred_palette <- c(strain=1L, oxygen=2L)
  assigned_palette <- integer(0)
  used_palette <- integer(0)
  for (fac in background_factors) {
    pref <- unname(preferred_palette[tolower(fac)])
    if (!length(pref) || is.na(pref) || pref %in% used_palette) {
      available <- setdiff(seq_along(factor_palette_bank), used_palette)
      pref <- if (length(available)) available[1] else
        ((length(used_palette) %% length(factor_palette_bank)) + 1L)
    }
    assigned_palette[fac] <- pref
    used_palette <- unique(c(used_palette, pref))
  }

  header_level_colors <- setNames(vector("list", length(background_factors)),
                                  background_factors)
  if (interaction_mode && length(background_factors)) {
    ordered_parts <- unname(header_parts[treatments])
    for (fac in background_factors) {
      pos <- match(fac, header_factor_order)
      levs <- unique(vapply(ordered_parts, function(z)
        if (length(z) >= pos) z[pos] else "", character(1)))
      levs <- levs[nzchar(levs)]
      base_pal <- factor_palette_bank[[assigned_palette[[fac]]]]
      cols <- if (length(levs) <= length(base_pal)) {
        base_pal[seq_along(levs)]
      } else {
        grDevices::colorRampPalette(base_pal)(length(levs))
      }
      header_level_colors[[fac]] <- setNames(cols, levs)
    }
  }

  factor_level_color <- function(factor, level, fallback="#666666") {
    if (is.null(factor) || !length(factor) || is.na(factor) || !nzchar(factor) ||
        !(factor %in% names(header_level_colors))) return(fallback)
    mp <- header_level_colors[[factor]]
    if (is.null(mp) || !length(mp) || is.na(level) || !nzchar(level) ||
        !(level %in% names(mp))) return(fallback)
    unname(mp[[level]])
  }

  header_text_color <- function(fill) {
    if (is.na(fill) || !nzchar(fill)) return("white")
    rgb <- tryCatch(grDevices::col2rgb(fill) / 255, error=function(e) NULL)
    if (is.null(rgb)) return("white")
    luminance <- 0.2126 * rgb[1,1] + 0.7152 * rgb[2,1] + 0.0722 * rgb[3,1]
    if (luminance > 0.62) "#111111" else "white"
  }

  if (interaction_mode && length(background_factors)) {
    inner_factor <- tail(background_factors, 1)
    inner_pos <- match(inner_factor, header_factor_order)
    color_keys <- vapply(unname(header_parts[treatments]), function(z)
      if (length(z) >= inner_pos) z[inner_pos] else z[length(z)], character(1))
    treat_colors <- setNames(vapply(color_keys, function(z)
      factor_level_color(inner_factor, z), character(1)), treatments)
  } else {
    pal_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                  "#F0E442", "#000000")
    level_colors <- setNames(rep(pal_cols, length.out=length(unique(treatments))),
                             unique(treatments))
    treat_colors <- setNames(unname(level_colors[treatments]), treatments)
  }

  cond_parts_png <- unname(header_parts[treatments])

  display_col_labels <- if (interaction_3way) {
    vapply(cond_parts_png, function(p) if (length(p) >= 3) p[length(p)] else p[length(p)], character(1))
  } else if (interaction_mode) {
    vapply(cond_parts_png, function(p) if (length(p) >= 2) p[length(p)] else p[1], character(1))
  } else {
    display_condition(treatments)
  }

  col_split_factor <- NULL
  top_title <- NULL
  top_ann <- NULL
  show_static_column_names <- TRUE

  if (interaction_mode && !interaction_3way) {
    top_levels <- vapply(cond_parts_png, function(p) p[1], character(1))
    second_levels <- vapply(cond_parts_png, function(p) if(length(p)>=2) p[2] else "", character(1))
    top_order <- unique(top_levels)
    top_align <- setNames(lapply(top_order, function(z) which(top_levels == z)), top_order)
    second_align <- setNames(lapply(seq_along(treatments), function(i) i), paste0("condition_", seq_along(treatments)))

    top_ann <- HeatmapAnnotation(
      Level1 = anno_block(
        align_to = top_align,
        panel_fun = function(index, nm) {
          grid.rect(gp=gpar(fill="white", col="#333333", lwd=0.95))
          grid.text(nm, gp=gpar(fontsize=9, fontface="bold"))
        }, height=unit(7.5,"mm")
      ),
      Level2 = anno_block(
        align_to = second_align,
        panel_fun = function(index, nm) {
          i <- index[1]
          second_factor <- if (length(background_factors)) background_factors[1] else ""
          fill_i <- factor_level_color(second_factor, second_levels[i])
          grid.rect(gp=gpar(fill=fill_i, col="white", lwd=0.8))
          grid.text(second_levels[i], gp=gpar(fontsize=6.6, fontface="bold",
                                               col=header_text_color(fill_i)))
        }, height=unit(7,"mm")
      ),
      show_annotation_name=FALSE, gap=unit(0.5,"mm")
    )
    show_static_column_names <- FALSE
  }

  if (interaction_3way) {
    primary_levels <- vapply(cond_parts_png, function(p) p[1], character(1))
    factor2_levels <- vapply(cond_parts_png, function(p) if(length(p)>=2) p[2] else "", character(1))
    factor3_levels <- vapply(cond_parts_png, function(p) if(length(p)>=3) p[3] else "", character(1))
    primary_order <- unique(primary_levels)
    primary_align <- setNames(lapply(primary_order, function(z) which(primary_levels == z)), primary_order)
    factor2_keys <- paste(primary_levels, factor2_levels, sep="
")
    factor2_order <- unique(factor2_keys)
    factor2_align <- setNames(lapply(factor2_order, function(z) which(factor2_keys == z)),
                              vapply(strsplit(factor2_order, "
", fixed=TRUE), `[`, character(1), 2))
    factor3_align  <- setNames(lapply(seq_along(treatments), function(i) i), paste0("factor3_", seq_along(treatments)))

    top_ann <- HeatmapAnnotation(
      Comparison = anno_block(
        align_to = primary_align,
        panel_fun = function(index, nm) {
          grid.rect(gp=gpar(fill="white", col="#333333", lwd=0.95))
          grid.text(nm, gp=gpar(fontsize=9, fontface="bold"))
        },
        height=unit(7.5,"mm")
      ),
      Factor2 = anno_block(
        align_to = factor2_align,
        panel_fun = function(index, nm) {
          factor2_name <- if (length(background_factors) >= 1L) background_factors[1] else ""
          fill_i <- factor_level_color(factor2_name, nm)
          grid.rect(gp=gpar(fill=fill_i, col="white", lwd=0.8))
          grid.text(nm, gp=gpar(fontsize=6.6, fontface="bold",
                                col=header_text_color(fill_i)))
        },
        height=unit(7,"mm")
      ),
      Factor3 = anno_block(
        align_to = factor3_align,
        panel_fun = function(index, nm) {
          factor3_name <- if (length(background_factors) >= 2L) background_factors[2] else
            if (length(background_factors)) background_factors[1] else ""
          level3 <- factor3_levels[index[1]]
          fill_i <- factor_level_color(factor3_name, level3)
          grid.rect(gp=gpar(fill=fill_i, col="white", lwd=0.8))
          grid.text(level3, gp=gpar(fontsize=6.2, fontface="bold",
                                    col=header_text_color(fill_i)))
        },
        height=unit(7,"mm")
      ),
      show_annotation_name=FALSE,
      gap=unit(c(0.5,0.4),"mm")
    )
    show_static_column_names <- FALSE
  }

  longest_treatment <- max(nchar(display_col_labels), na.rm=TRUE)
  min_interaction_col_mm <- if (interaction_mode) 14.2 else 9.0
  column_width_mm <- max(min_interaction_col_mm + 2.4, longest_treatment * 0.82 + 5.0)
  heatmap_body_width <- unit(n_treat * column_width_mm, "mm")

  hc_display <- if (!is.null(hc_rows)) stats::as.dendrogram(hc_rows) else NULL

  ht <- Heatmap(nes_display_mat, name="NES", col=col_fn, na_col="#D0D0D0",
                width=heatmap_body_width,
                column_labels=display_col_labels,
                top_annotation=top_ann,
                column_title=top_title,
                column_title_gp=gpar(fontsize=11, fontface="bold"),
                column_title_side="top",
                column_split=col_split_factor,
                cluster_rows=if(!is.null(hc_display)) hc_display else (n_terms>=2),
                cluster_columns=FALSE, row_dend_side="left",
                row_dend_width=unit(21,"mm"),
                show_row_names=FALSE,
                show_column_names=show_static_column_names, column_names_side="top", column_names_rot=0,
                column_names_gp=gpar(fontsize=6.4, fontface="bold"),
                column_names_centered=TRUE,
                layer_fun=lf_img, rect_gp=gpar(col="#D9D9D9", lwd=0.35),
                right_annotation=right_ann,
                heatmap_legend_param=list(
                  title="NES", title_gp=gpar(fontsize=8, fontface="bold"),
                  labels_gp=gpar(fontsize=7), legend_height=unit(3,"cm"),
                  at=round(seq(-max_abs, max_abs, length.out=5), 1)))
  
  fig_h <- max(10, n_terms*0.34 + 5.0)
  interaction_width_cm <- if (show_formal_int_cols) max(0, n_int_cols * 1.05) else 0.8
  legend_reserve_cm <- 5.8
  fig_w <- max(27, 6.3 + n_treat*column_width_mm/10 + interaction_width_cm +
                 term_label_width/10 + legend_reserve_cm)
  draw_static_heatmap <- function() {
    draw(ht, merge_legend=TRUE, heatmap_legend_side="right",
         annotation_legend_side="right", padding=unit(c(6,18,30,6),"mm"))
    if (show_formal_int_cols && length(all_int_cols) > 0L) {
      for (ii in seq_along(all_int_cols)) {
        nm <- unname(int_col_labels[all_int_cols[ii]])
        try(decorate_annotation(nm, {
          if (ii == 1L) grid.lines(unit(c(0,0),"npc"), unit(c(0,1),"npc"), gp=gpar(col="#888888",lwd=0.55))
          grid.lines(unit(c(1,1),"npc"), unit(c(0,1),"npc"), gp=gpar(col="#888888",lwd=0.55))
          if (n_terms > 0L) {
            for (yy in seq(0, 1, length.out=n_terms + 1L))
              grid.lines(unit(c(0,1),"npc"), unit(c(yy,yy),"npc"),
                         gp=gpar(col="#D9D9D9",lwd=0.35))
          }
        }), silent=TRUE)
      }
    }
  }

  png_path <- file.path(output_dir,paste0(type_suffix,"_heatmap.png"))
  png(png_path,width=fig_w,height=fig_h,units="cm",res=600,type="cairo-png")
  draw_static_heatmap()
  dev.off()
  
  
  svg_inner_labels <- if (interaction_mode) {
    vapply(treatments, function(t) {
      p <- split_condition_components(t)
      if (length(p) >= 2) p[length(p)] else p[1]
    }, character(1))
  } else display_condition(treatments)
  max_header_chars <- max(nchar(svg_inner_labels), na.rm=TRUE)
  cell_w <- as.integer(max(if (interaction_mode) 88 else 66, ceiling(max_header_chars * 5.6 + 22)))
  cell_h <- 25L
  dend_w <- 128L
  if (interaction_3way && n_int_cols == 0L && length(int_factors_from_manifest) >= 2) {
    pfx <- tolower(int_label_prefix)
    background_order <- setdiff(header_factor_order, contrast_factor)
    facs <- background_order[background_order %in% int_factors_from_manifest]
    if (length(facs) < 2L) facs <- unique(c(facs, int_factors_from_manifest))
    f1 <- facs[1]; f2 <- facs[2]
    placeholder_int_cols   <- c(paste0("interaction_x_", f1),
                                paste0("interaction_x_", f2),
                                "interaction_3way")
    placeholder_int_labels <- c(paste0(tolower(contrast_factor), " × ", f1),
                                paste0(tolower(contrast_factor), " × ", f2),
                                paste0(tolower(contrast_factor), " × ", f1, " × ", f2))
    placeholder_int_labels_full <- c(paste0(pfx, " × ", f1),
                                     paste0(pfx, " × ", f2),
                                     paste0(pfx, " × ", f1, " × ", f2))
    all_int_cols   <- placeholder_int_cols
    n_int_cols     <- 3L
    int_col_labels <- setNames(placeholder_int_labels, placeholder_int_cols)
    int_col_labels_full <- setNames(placeholder_int_labels_full, placeholder_int_cols)
    int_sig_vecs   <- setNames(lapply(placeholder_int_cols, function(col)
      setNames(rep(NA_character_, nrow(summ)), summ$ID)), placeholder_int_cols)
    cat("  Interaction columns: unavailable in the treatment summary\n")
  }
  
  n_display_int <- length(all_int_cols)
  star_w  <- if (show_formal_int_cols) as.integer(max(1L, n_display_int) * 80L) else 40L
  header_h <- if (interaction_3way) 126L else if (interaction_mode) 112L else 78L
  dot_r <- 5L; label_w <- 490L; badge_col_w <- 42L; margin_r <- 16L
  cells_start <- dend_w
  star_x      <- dend_w + n_treat * cell_w
  labels_x    <- star_x + star_w + 4L
  svg_w       <- labels_x + label_w + badge_col_w + margin_r
  
  n_go   <- length(go_order)
  svg_h  <- header_h+n_go*cell_h+60L
  leaf_y <- setNames(header_h+(seq_along(go_order)-0.5)*cell_h,go_order)
  
  n_nodes <- 0L; node_y <- numeric(0); node_x <- numeric(0)
  node_members_go <- list(); merge_m <- matrix(0L,0,2); heights <- numeric(0)
  entity_y <- numeric(0); labels_h_dend <- character(0)
  node_lca_desc <- character(0)
  
  if (!is.null(hc_rows)&&n_go>=2) {
    labels_h_dend <- hc_rows$labels; n_leaf_dend <- length(labels_h_dend)
    merge_m <- hc_rows$merge
    n_nodes <- nrow(merge_m)
    heights <- hc_rows$height
    h_max <- max(heights, na.rm=TRUE)
    entity_y <- vapply(seq_len(n_leaf_dend),function(i) {
      lbl <- labels_h_dend[i]; if(lbl %in% names(leaf_y)) leaf_y[lbl] else header_h
    },numeric(1))
    node_y <- numeric(n_nodes); node_x <- numeric(n_nodes)
    node_members_go <- vector("list",n_nodes)
    for (k in seq_len(n_nodes)) {
      li <- merge_m[k,1]; ri <- merge_m[k,2]
      ly <- if(li<0) entity_y[-li] else node_y[li]; ry <- if(ri<0) entity_y[-ri] else node_y[ri]
      node_y[k] <- (ly+ry)/2
      hk <- max(0, heights[k])
      h_draw <- if (is.finite(h_max) && h_max > 0) hk / h_max else 0
      node_x[k] <- dend_w - h_draw * (dend_w - 14)
      lm <- if(li<0) labels_h_dend[-li] else node_members_go[[li]]
      rm <- if(ri<0) labels_h_dend[-ri] else node_members_go[[ri]]
      node_members_go[[k]] <- unique(c(lm,rm))
    }
    node_lca_desc <- character(n_nodes)
    node_category_id <- rep(NA_character_, n_nodes)
    node_category_name <- rep(NA_character_, n_nodes)
    for (k in seq_len(n_nodes)) {
      members <- node_members_go[[k]]
      category_id <- if (!is.null(.anc_cache_global) && type_suffix=="GO_BP")
        find_cluster_category(members, .anc_cache_global) else NA_character_
      category_name <- if (!is.na(category_id) && nzchar(category_id)) {
        tryCatch({
          suppressPackageStartupMessages(library(GO.db))
          go_obj <- GO.db::GOTERM[[category_id]]
          nm <- if (is.null(go_obj)) NA_character_ else as.character(methods::slot(go_obj, "Term"))
          if (is.na(nm) || !nzchar(trimws(nm))) category_id else paste0(trimws(nm), " (", category_id, ")")
        }, error=function(e) category_id)
      } else {
        nm <- unname(desc_lu[members[1]])
        if (is.na(nm) || !nzchar(trimws(nm))) members[1] else trimws(nm)
      }
      node_category_id[k] <- category_id
      node_category_name[k] <- category_name
    }
    keep_node <- rep(TRUE, n_nodes)
    valid_ids <- unique(node_category_id[!is.na(node_category_id) & nzchar(node_category_id)])
    for (cid in valid_ids) {
      idx <- which(node_category_id == cid)
      if (length(idx) > 1L) {
        sizes <- vapply(idx, function(z) length(node_members_go[[z]]), integer(1))
        chosen <- idx[order(sizes, idx)][1]
        keep_node[setdiff(idx, chosen)] <- FALSE
      }
    }
    for (k in seq_len(n_nodes)) {
      if (keep_node[k]) node_lca_desc[k] <- paste0(node_category_name[k], " — ", length(node_members_go[[k]]), " terms")
    }
  }
  
  TAB20 <- c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a',
             '#d62728','#ff9896','#9467bd','#c5b0d5','#8c564b','#c49c94',
             '#e377c2','#f7b6d2','#7f7f7f','#c7c7c7','#bcbd22','#dbdb8d','#17becf','#9edae5')
  
  grid_svg <- paste(c(
    vapply(seq_along(go_order),function(ri){y<-header_h+(ri-1L)*cell_h
    paste0("<line x1='",cells_start,"' y1='",y,"' x2='",star_x,"' y2='",y,"' stroke='#f0f0f0' stroke-width='0.5'/>")}
    ,character(1)),
    vapply(seq_along(treatments),function(ci){x<-cells_start+(ci-1L)*cell_w
    paste0("<line x1='",x,"' y1='",header_h,"' x2='",x,"' y2='",header_h+n_go*cell_h,"' stroke='#f0f0f0' stroke-width='0.5'/>")}
    ,character(1))
  ),collapse="\n")
  
  int_cell_w <- if (n_display_int > 0) as.integer(star_w / max(1L, n_display_int)) else star_w
  
  header_svg <- if (interaction_3way) {
    cond_parts_list <- unname(header_parts[treatments])
    top_levels    <- unique(vapply(cond_parts_list, `[`, character(1), 1))
    mid_levels    <- unique(vapply(cond_parts_list, `[`, character(1), 2))
    inner_labels  <- vapply(cond_parts_list, function(p) paste(p[-(1:2)], collapse=":"), character(1))
    
    h1 <- as.integer(header_h * 0.30)
    h2 <- as.integer(header_h * 0.35)
    h3 <- as.integer(header_h - h1 - h2)
    
    rects <- character(0)
    for (top in top_levels) {
      top_idx <- which(vapply(cond_parts_list, function(p) p[1] == top, logical(1)))
      x1 <- cells_start + (min(top_idx)-1L)*cell_w
      w1 <- length(top_idx)*cell_w
      rects <- c(rects,
                 paste0("<rect x='",x1,"' y='0' width='",w1,"' height='",h1,
                        "' fill='white' stroke='#333333' stroke-width='1'/>") ,
                 paste0("<text x='",x1+w1/2L,"' y='",round(h1*0.67),
                        "' text-anchor='middle' font-size='11' fill='#222222' font-weight='bold'>",
                        htmlesc(top),"</text>"))
    }
    mid_keys <- vapply(cond_parts_list, function(p) paste(p[1], p[2], sep="
"), character(1))
    for (key in unique(mid_keys)) {
      mid_idx <- which(mid_keys == key)
      mid <- strsplit(key, "
", fixed=TRUE)[[1]][2]
      x2 <- cells_start + (min(mid_idx)-1L)*cell_w
      w2 <- length(mid_idx)*cell_w
      factor2_name <- if (length(background_factors) >= 1L) background_factors[1] else ""
      bc <- factor_level_color(factor2_name, mid)
      txt_col <- header_text_color(bc)
      rects <- c(rects,
                 paste0("<rect x='",x2,"' y='",h1,"' width='",w2,"' height='",h2,"' fill='",bc,"' stroke='white' stroke-width='1'/>") ,
                 paste0("<text x='",x2+w2/2L,"' y='",h1+round(h2*0.67),"' text-anchor='middle' font-size='9' fill='",txt_col,"' font-weight='bold'>",htmlesc(mid),"</text>"))
    }
    for (ci in seq_along(treatments)) {
      t  <- treatments[ci]; x3 <- cells_start+(ci-1L)*cell_w
      lbl <- inner_labels[ci]
      factor3_name <- if (length(background_factors) >= 2L) background_factors[2] else
        if (length(background_factors)) background_factors[1] else ""
      bc <- factor_level_color(factor3_name, lbl)
      txt_col <- header_text_color(bc)
      rects <- c(rects,
                 paste0("<rect x='",x3,"' y='",h1+h2,"' width='",cell_w,"' height='",h3,"' fill='",bc,"' stroke='white' stroke-width='1'/>") ,
                 paste0("<text x='",x3+cell_w/2L,"' y='",h1+h2+round(h3*0.72),
                        "' text-anchor='middle' font-size='8' fill='",txt_col,"' font-weight='bold'>",
                        htmlesc(lbl),"</text>"))
    }
    int_x0    <- star_x
    col_w_hdr <- as.integer(star_w / max(1L, n_display_int))
    for (k in seq_along(all_int_cols)) {
      xk    <- int_x0 + (k-1L) * col_w_hdr
      lbl_k <- int_col_labels[all_int_cols[k]]
      xc    <- xk + col_w_hdr / 2L
      rects <- c(rects,
                 paste0("<rect x='",xk,"' y='0' width='",col_w_hdr,"' height='",header_h,"' fill='#444'/>"),
                 render_int_header(lbl_k, xc, header_h),
                 if (k < n_display_int)
                   paste0("<line x1='",xk+col_w_hdr,"' y1='4' x2='",xk+col_w_hdr,"' y2='",header_h-4L,
                          "' stroke='#888' stroke-width='1'/>")
                 else "")
    }
    paste(rects, collapse="\n")
  } else if (interaction_mode) {
    parts2 <- unname(header_parts[treatments])
    top2 <- vapply(parts2, function(p) p[1], character(1))
    inner2 <- vapply(parts2, function(p) if(length(p)>=2) p[2] else p[1], character(1))
    hcap <- as.integer(header_h * 0.42)
    hinner <- header_h - hcap
    rects2 <- character(0)
    for (cap_lbl in unique(top2)) {
      idx <- which(top2 == cap_lbl)
      x0 <- cells_start + (min(idx)-1L)*cell_w
      ww <- length(idx)*cell_w
      rects2 <- c(rects2,
        paste0("<rect x='",x0,"' y='0' width='",ww,"' height='",hcap,"' fill='white' stroke='#333' stroke-width='1'/>",
               "<text x='",x0+ww/2,"' y='",round(hcap*0.68),"' text-anchor='middle' font-size='10' fill='#222' font-weight='bold'>",htmlesc(cap_lbl),"</text>"))
    }
    for (ci in seq_along(treatments)) {
      t <- treatments[ci]; x <- cells_start+(ci-1L)*cell_w
      second_factor <- if (length(background_factors)) background_factors[1] else ""
      bc <- factor_level_color(second_factor, inner2[ci])
      txt_col <- header_text_color(bc)
      rects2 <- c(rects2,
        paste0("<rect x='",x,"' y='",hcap,"' width='",cell_w,"' height='",hinner,"' fill='",bc,"' stroke='white' stroke-width='1'/>",
               "<text x='",x+cell_w/2,"' y='",hcap+round(hinner*0.62),"' text-anchor='middle' font-size='7' fill='",txt_col,"' font-weight='bold'>",htmlesc(inner2[ci]),"</text>"))
    }
    if (show_formal_int_cols && length(all_int_cols)) {
      col_w_hdr <- as.integer(star_w / max(1L, n_display_int))
      for (k in seq_along(all_int_cols)) {
        xk <- star_x + (k-1L) * col_w_hdr
        lbl_k <- int_col_labels[all_int_cols[k]]
        rects2 <- c(rects2,
          paste0("<rect x='",xk,"' y='0' width='",col_w_hdr,"' height='",header_h,"' fill='#444'/>") ,
          render_int_header(lbl_k, xk + col_w_hdr/2L, header_h),
          if (k < n_display_int) paste0("<line x1='",xk+col_w_hdr,"' y1='4' x2='",xk+col_w_hdr,"' y2='",header_h-4L,"' stroke='#888' stroke-width='1'/>") else "")
      }
    } else {
      rects2 <- c(rects2, paste0("<rect x='",star_x,"' y='0' width='",star_w,"' height='",header_h,"' fill='white' stroke='none'/>") )
    }
    paste(rects2, collapse="
")
  } else {
    paste(c(
      vapply(seq_along(treatments),function(ci){
        t<-treatments[ci];x<-cells_start+(ci-1L)*cell_w;bc<-treat_colors[t]
        paste0("<rect x='",x,"' y='0' width='",cell_w,"' height='",header_h,"' fill='",bc,"'/>",
               "<text x='",x+cell_w/2L,"' y='",header_h/2L+5L,
               "' text-anchor='middle' font-size='9' fill='white' font-weight='bold'>",
               htmlesc(display_condition(t)),"</text>")
      },character(1)),
      paste0("<rect x='",star_x,"' y='0' width='",star_w,"' height='",header_h,"' fill='#555'/>",
             "<text x='",star_x+star_w/2L,"' y='",header_h/2L-8L,
             "' text-anchor='middle' font-size='7' fill='white' font-weight='bold'>",
             htmlesc(format_factor_label(contrast_factor)),"</text>",
             paste0("<text x='",star_x+star_w/2L,"' y='",header_h/2L+6L,
                    "' text-anchor='middle' font-size='7' fill='white'>differences</text>"))
    ), collapse="\n")
  }
  
  cells_svg <- paste(unlist(lapply(seq_along(go_order),function(ri){
    go_id<-go_order[ri];y<-header_h+(ri-1L)*cell_h
    vapply(seq_along(treatments),function(ci){
      t<-treatments[ci];x<-cells_start+(ci-1L)*cell_w
      nv<-nes_mat[go_id,t];sig<-sig_mat[go_id,t];if(is.na(sig)) sig<-""
      bg<-cfn_svg(nv,sig,max_abs)
      nes_txt<-if(!is.na(nv)&&sig=="*") sprintf("%.2f",nv) else ""
      ns_txt <-if(!is.na(nv)&&sig=="ns") sprintf("%.2f",nv) else ""
      nes_txt_col <- if(!is.na(nv)&&sig=="*"&&abs(nv/max_abs)>=0.42) "#FFFFFF" else "#333333"
      cld_ltr<-if(use_cld&&!is.null(display_cld_mat)){
        l<-display_cld_mat[go_id,t];if(!is.na(l)&&nchar(trimws(l))>0) trimws(l) else ""} else ""
      paste0("<rect x='",x,"' y='",y,"' width='",cell_w,"' height='",cell_h,"' fill='",bg,"'/>",
             if(nchar(nes_txt)>0) paste0("<text x='",x+cell_w/2L,"' y='",y+round(cell_h*0.47),
                                         "' text-anchor='middle' font-size='8' fill='",nes_txt_col,"'>",nes_txt,"</text>") else
                                           if(nchar(ns_txt)>0) paste0("<text x='",x+cell_w/2L,"' y='",y+round(cell_h*0.55),
                                                                      "' text-anchor='middle' font-size='7' fill='#999' font-style='italic'>",ns_txt,"</text>") else "",
             if(nchar(cld_ltr)>0) paste0("<text x='",x+cell_w/2L,"' y='",y+round(cell_h*0.88),
                                         "' text-anchor='middle' font-size='8' fill='#333'>",cld_ltr,"</text>") else "")
    },character(1))
  })),collapse="\n")
  
  stars_svg <- paste(vapply(seq_along(go_order),function(ri){
    go_id<-go_order[ri];y<-header_h+(ri-1L)*cell_h
    parts <- paste0("<rect x='",star_x,"' y='",y,"' width='",star_w,"' height='",cell_h,"' fill='#F8F8F8'/>")
    
    if (!show_formal_int_cols) {
      da<-diff_col_vec[go_id];if(is.na(da)) da<-""
      if(nchar(da)>0)
        parts <- paste0(parts,"<text x='",star_x+star_w/2L,"' y='",y+round(cell_h*0.72),
                        "' text-anchor='middle' font-size='11' fill='#C92A2A' font-weight='bold'>",htmlesc(da),"</text>")
    } else {
      n_int_total_svg <- max(1L, n_display_int)
      col_w_each_svg  <- as.integer(star_w / n_int_total_svg)
      three_sig <- if (has_3way_col) any(vapply(int_3way_cols, function(col) {
        v <- int_sig_vecs[[col]][go_id]
        !is.null(v) && !is.na(v) && v == "*"
      }, logical(1))) else FALSE
      
      for (k in seq_along(all_int_cols)) {
        icol  <- all_int_cols[k]
        xk    <- star_x + (k-1L) * col_w_each_svg
        sv    <- int_sig_vecs[[icol]][go_id]
        is_sig <- !is.null(sv) && !is.na(sv) && sv == "*"
        is_3way_col <- icol %in% int_3way_cols
        
        display_char <- if (!is_3way_col && three_sig) "" else
          if (is_sig) "*" else ""
        txt_col <- if (display_char == "*") "'#C92A2A'" else "'#BBBBBB'"
        fw <- if (display_char == "*") "bold" else "normal"
        
        if (nchar(display_char) > 0)
          parts <- paste0(parts,
                          "<text x='",xk+col_w_each_svg/2L,"' y='",y+round(cell_h*0.72),
                          "' text-anchor='middle' font-size='11' fill=",txt_col,
                          " font-weight='",fw,"'>",display_char,"</text>")
      }
    }
    parts
  },character(1)),collapse="\n")

  star_row_lines_svg <- paste(vapply(seq_len(n_go + 1L), function(ri) {
    y <- header_h + (ri - 1L) * cell_h
    paste0("<line x1='", star_x, "' y1='", y, "' x2='", star_x + star_w,
           "' y2='", y, "' stroke='#D9D9D9' stroke-width='0.5'/>")
  }, character(1)), collapse="\n")
  
  display_labels_html <- vapply(go_order, function(go_id) {
    desc <- if (go_id %in% names(desc_lu)) desc_lu[go_id] else go_id
    lbl  <- paste0(desc, " (", go_id, ")")
    if (nchar(lbl) > 70) lbl <- paste0(substr(lbl, 1, 68), "...")
    lbl
  }, character(1))
  longest_label_chars <- max(nchar(display_labels_html), na.rm = TRUE)
  badge_x_common <- labels_x + min(label_w - 38L,
                                   max(126L, as.integer(longest_label_chars * 5.35 + 10L)))
  
  labels_svg <- paste(vapply(seq_along(go_order),function(ri){
    go_id<-go_order[ri];y<-header_h+(ri-1L)*cell_h
    lbl<-display_labels_html[go_id]
    gs<-htmlesc(go_id); n_elim<-length(elim_map[[go_id]])
    base<-paste0("<rect id='lb-",gs,"' x='",labels_x,"' y='",y,"' width='",label_w,"' height='",cell_h,"' fill='#FAFAFA'/>",
                 "<text x='",labels_x+6L,"' y='",y+round(cell_h*0.72),
                 "' font-size='9' fill='#111' class='lbl-txt' data-goid='",gs,"' style='cursor:pointer'>",htmlesc(lbl),"</text>")
    badge<-if(n_elim>0){
      bx <- badge_x_common
      paste0("<rect x='",bx,"' y='",y+4L,"' width='36' height='",cell_h-8L,
             "' rx='3' fill='#1971C2' class='elim-btn' data-goid='",gs,"' style='cursor:pointer'/>",
             "<text x='",bx+18L,"' y='",y+round(cell_h*0.70),
             "' text-anchor='middle' font-size='9' fill='white' pointer-events='none'>+",n_elim,"</text>")
    } else ""
    paste0(base,badge)
  },character(1)),collapse="\n")
  
  dend_svg <- ""
  if (n_nodes>0L) {
    dend_lines <- unlist(lapply(seq_len(n_nodes),function(k){
      li<-merge_m[k,1];ri<-merge_m[k,2]
      lx<-if(li<0) dend_w else node_x[li];ly2<-if(li<0) entity_y[-li] else node_y[li]
      rx<-if(ri<0) dend_w else node_x[ri];ry2<-if(ri<0) entity_y[-ri] else node_y[ri]
      nx<-node_x[k];ny<-node_y[k]; members<-node_members_go[[k]]
      ttl<-if(k<=length(node_lca_desc)&&nchar(node_lca_desc[k])>0) htmlesc(node_lca_desc[k]) else ""
      dot<-if(type_suffix == "GO_BP" && length(members)>=2 && nchar(ttl)>0) paste0(
        "<circle cx='",nx,"' cy='",ny,"' r='",dot_r,
        "' fill='#999' stroke='white' stroke-width='1.5' class='node-dot'",
        " data-lca='node-",k,"' data-leaves='",htmlesc(paste(members,collapse=",")),"'",
        " data-desc='",ttl,"'><title>",ttl,"</title></circle>") else ""
      c(paste0("<line x1='",nx,"' y1='",ly2,"' x2='",nx,"' y2='",ry2,"' stroke='#555' stroke-width='1.25'/>"),
        paste0("<line x1='",nx,"' y1='",ly2,"' x2='",lx,"' y2='",ly2,"' stroke='#555' stroke-width='1.25'/>"),
        paste0("<line x1='",nx,"' y1='",ry2,"' x2='",rx,"' y2='",ry2,"' stroke='#555' stroke-width='1.25'/>"),
        dot)
    }))
    dend_svg <- paste(dend_lines,collapse="\n")
  }
  
  GRP_COLORS <- c("#4DABF7","#FF8787","#63E6BE","#FFA94D","#9775FA","#F783AC","#74C0FC","#FFC9C9")
  
  mk_th <- function(txt,bcol) paste0(
    "<th style='background:",bcol,";color:white;font-size:8px;padding:3px 5px;text-align:center;white-space:nowrap'>",txt,"</th>")
  
  tbl_hdr <- paste0(
    "<thead><tr>",
    "<th style='background:#444;color:white;font-size:8px;padding:3px 5px;text-align:left'>GO term</th>",
    paste(vapply(treatments,function(t) mk_th(paste0("NES<br><small>",htmlesc(t),"</small>"),treat_colors[t]),character(1)),collapse=""),
    if(use_cld) paste(vapply(treatments,function(t) mk_th(paste0("CLD<br><small>",htmlesc(t),"</small>"),"#888"),character(1)),collapse="") else "",
    if (interaction_mode && length(all_int_cols) > 0) {
      paste(vapply(all_int_cols, function(col) {
        lbl <- int_col_labels[col]
        mk_th(htmlesc(lbl), "#555")
      }, character(1)), collapse="")
    } else mk_th("Diff","#555"),
    "</tr></thead>")
  
  make_row_tds <- function(gid,winner=FALSE) {
    rb<-if(winner)"background:#EFF6FF;" else ""
    nes_tds <- paste(vapply(treatments,function(t){
      val<-if(gid %in% rownames(nes_mat_all)) nes_mat_all[gid,t] else NA_real_
      sig<-if(gid %in% rownames(sig_mat_all)) sig_mat_all[gid,t] else "";if(is.na(sig)) sig<-""
      bg<-cfn_svg(val,sig,max_abs); txt<-if(!is.na(val)) sprintf("%.2f",val) else "—"
      it<-if(sig=="ns")" font-style:italic;color:#999" else ""
      paste0("<td style='",rb,"background:",bg,";text-align:center;font-size:10px;padding:3px;",it,"'>",txt,"</td>")
    },character(1)),collapse="")
    cld_tds<-if(use_cld&&!is.null(cld_mat)) paste(vapply(treatments,function(t){
      l<-if(!is.null(cld_mat_all) && gid %in% rownames(cld_mat_all)) display_cld_mat_all[gid,t] else "";if(is.na(l)) l<-""
      paste0("<td style='",rb,"text-align:center;font-size:10px;padding:3px'>",trimws(l),"</td>")
    },character(1)),collapse="") else ""
    int_tds <- if (interaction_mode && length(all_int_cols) > 0) {
      three_sig_popup <- if (has_3way_col) {
        any(vapply(int_3way_cols, function(col) {
          v <- int_sig_vecs[[col]][gid]
          !is.null(v) && !is.na(v) && v == "*"
        }, logical(1)))
      } else FALSE
      paste(vapply(all_int_cols, function(col) {
        is_3way_col <- col %in% int_3way_cols
        v <- if (gid %in% names(int_sig_vecs[[col]])) int_sig_vecs[[col]][gid] else ""
        if (is.null(v) || is.na(v)) v <- ""
        disp <- if (!is_3way_col && three_sig_popup) "" else v
        paste0("<td style='",rb,"text-align:center;font-size:11px;padding:3px;",
               if(disp=="*") "color:#C92A2A;font-weight:bold" else "color:#aaa",
               "'>",htmlesc(disp),"</td>")
      }, character(1)), collapse="")
    } else {
      dv <- diff_col_vec[gid]; if(is.na(dv)) dv <- ""
      paste0("<td style='",rb,"text-align:center;font-size:11px;padding:3px;color:#C92A2A;font-weight:bold'>",htmlesc(dv),"</td>")
    }
    paste0(nes_tds,cld_tds,int_tds)
  }
  
  panels_html <- paste(vapply(go_order,function(go_id){
    gs<-htmlesc(go_id); desc<-if(go_id %in% names(desc_lu)) desc_lu[go_id] else go_id
    da<-diff_col_vec[go_id];if(is.na(da)) da<-""
    losers<-elim_map[[go_id]]
    gcol<-GRP_COLORS[((match(go_id,selected_ids,nomatch=1L)-1L)%%length(GRP_COLORS))+1L]
    loser_rows<-if(!is.null(losers)&&length(losers)>0) paste(vapply(losers,function(lg){
      sim_val<-if(!is.null(similarity_matrix_full)&&lg %in% rownames(similarity_matrix_full)&&go_id %in% rownames(similarity_matrix_full))
        round(similarity_matrix_full[go_id,lg],3) else NA
      dsc<-if(lg %in% names(desc_lu)) desc_lu[lg] else lg
      paste0("<tr><td style='padding:2px 5px;font-size:9px;border-left:3px solid ",gcol,
             ";background:",gcol,"11'>",htmlesc(dsc),
             "<br><small style='color:#aaa'>",htmlesc(lg),
             if(!is.na(sim_val)) paste0(" similarity=",sim_val) else "","</small></td>",
             make_row_tds(lg,FALSE),"</tr>")
    },character(1)),collapse="") else ""
    le_btn<-paste0("<button onclick=\"showLeadingEdge('",gs,"')\"",
                   " style='float:right;flex-shrink:0;background:#2ca02c;color:white;border:none;",
                   "border-radius:4px;padding:5px 11px;font-size:9px;cursor:pointer'>&#9654; Leading edge</button>")
    paste0("<div id='pd-",gs,"' style='display:none'>",
           "<div style='display:flex;justify-content:space-between;align-items:flex-start;gap:10px;margin-bottom:8px'>",
           "<div><b style='font-size:12px'>",htmlesc(desc),"</b>",
           "<br><small style='color:#888'>",htmlesc(go_id),"</small>",
           if(nchar(da)>0) paste0(" <b style='color:#C92A2A'>",htmlesc(da),"</b>") else "",
           "</div>",le_btn,"</div>",
           "<p style='font-size:8px;color:#888;margin:0 0 6px'>Grey cell = not sig. vs ref</p>",
           "<div style='overflow-x:auto'><table style='border-collapse:collapse;font-size:10px;width:100%'>",
           tbl_hdr,"<tbody>",
           "<tr><td style='padding:2px 5px;font-size:9px;font-weight:bold;border-left:3px solid ",gcol,
           ";background:#EFF6FF'>",htmlesc(desc),"<br><small style='color:#aaa'>",htmlesc(go_id)," [rep]</small></td>",
           make_row_tds(go_id,TRUE),"</tr>",
           if(!is.null(losers)&&length(losers)>0)
             paste0("<tr><td colspan='",n_treat+(if(use_cld) n_treat else 0)+2,
                    "' style='padding:3px 5px;font-size:8px;font-weight:bold;color:#666;border-top:2px solid #ddd'>",
                    "Eliminated (",length(losers),"):</td></tr>",loser_rows) else "",
           "</tbody></table></div></div>")
  },character(1)),collapse="\n")
  
  all_go_js <- unique(c(selected_ids,unlist(elim_map,use.names=FALSE)))
  
  go_data_js <- paste0("{",paste(vapply(all_go_js,function(g){
    d<-jesc(if(g %in% names(desc_lu)) desc_lu[g] else g)
    np<-paste(vapply(treatments,function(t){
      v<-if(g %in% rownames(nes_mat_all)) nes_mat_all[g,t] else NA_real_
      s<-if(g %in% rownames(sig_mat_all)) sig_mat_all[g,t] else "";if(is.na(s)) s<-""
      paste0('"',jesc(t),'":',if(is.na(v))"null" else sprintf("%.4f",v),
             ',"',jesc(t),'_sig":',jstr(s))
    },character(1)),collapse=",")
    dv<-diff_col_vec[g];if(is.na(dv)) dv<-""
    paste0('"',g,'":{"desc":"',d,'",',np,',"diff":',jstr(dv),'}')
  },character(1)),collapse=","),"}")
  
  leading_js <- paste0("{",paste(vapply(go_order,function(g){
    cond_parts<-paste(vapply(treatments,function(t){
      orfs<-get_le(g,t)
      if(!length(orfs)) return(paste0('"',t,'":null'))
      gene_objs<-paste(vapply(orfs,function(orf){
        nm<-disp(orf)
        pids<-names(orf_pathways[[orf]])
        pids_j<-if(length(pids)>0) paste0('["',paste(jesc(pids),collapse='","'),'"]') else "[]"
        lfc_v<-NA_real_;padj_v<-NA_real_
        lk <- de_lookup[[t]]
        if(!is.null(lk)){
          v_lfc  <- lk$lfc[orf];  if(!is.null(v_lfc)  && !is.na(v_lfc))  lfc_v  <- v_lfc
          v_padj <- lk$padj[orf]; if(!is.null(v_padj) && !is.na(v_padj)) padj_v <- v_padj
        }
        gene_url <- if (enable_sgd_gene_links)
          paste0("https://www.yeastgenome.org/locus/", orf) else ""
        sid_status <- if (enable_sgd_pathway_lookup) {
          sids <- sgd_ids[[orf]]
          if (is.null(sids) || !length(sids)) {
            "unknown"
          } else {
            sts <- vapply(sids, function(sid) {
              st <- .SGD_STATUS[[sid]]
              if (is.null(st)) "unknown" else as.character(st)
            }, character(1))
            if (any(startsWith(sts, "request_failed:"))) "request_failed"
            else paste(unique(sts), collapse="|")
          }
        } else "disabled"
        paste0('{"orf":',jstr(orf),',"name":',jstr(nm),
               ',"url":',jstr(gene_url),
               ',"paths":',pids_j,
               ',"sgdStatus":',jstr(sid_status),
               ',"padj":',if(is.finite(padj_v)) sprintf("%.6f",padj_v) else "null",
               ',"logfc":',if(is.finite(lfc_v)) sprintf("%.4f",lfc_v) else "null",'}')
      },character(1)),collapse=",")
      path_map<-list()
      for(orf in orfs) for(pid in names(orf_pathways[[orf]])){
        pw<-orf_pathways[[orf]][[pid]]
        if(is.null(path_map[[pid]])) path_map[[pid]]<-list(id=pw$id,name=pw$name,url=pw$url,orfs=character(0),n=0L)
        if(!orf %in% path_map[[pid]]$orfs){path_map[[pid]]$orfs<-c(path_map[[pid]]$orfs,orf);path_map[[pid]]$n<-path_map[[pid]]$n+1L}
      }
      path_objs<-if(length(path_map)>0) paste(vapply(names(path_map),function(pid){
        pw<-path_map[[pid]]
        gn<-paste0('["',paste(jesc(vapply(pw$orfs,disp,character(1))),collapse='","'),'"]')
        paste0(jstr(pid),':{\"id\":',jstr(pw$id),',"name":',jstr(pw$name),',"url":',jstr(pw$url),',"n":',pw$n,',"geneNames\":',gn,'}')
      },character(1)),collapse=",") else ""
      paste0('"',jesc(t),'":{"genes":[',gene_objs,'],"pathways":{',path_objs,'}}')
    },character(1)),collapse=",")
    paste0('"',g,'":{'  ,cond_parts,'}')
  },character(1)),collapse=","),"}")
  
  volcano_js <- paste0("{",paste(vapply(treatments,function(t){
    df_t <- de_data[[t]]
    if(is.null(df_t)) return(paste0('"',jesc(t),'":null'))
    dv <- get_de_vectors(df_t)
    if(is.null(dv) || !nrow(dv)) return(paste0('"',jesc(t),'":null'))
    dnames <- vapply(dv$orf, disp, character(1))
    dup_orf <- duplicated(dv$orf) | duplicated(dv$orf, fromLast=TRUE)
    dnames[dup_orf] <- dv$input_gene[dup_orf]
    paste0('"',jesc(t),'":{"orfs":[',paste0('"',jesc(dv$orf),'"',collapse=","),'],',
           '"names":[',paste0('"',jesc(dnames),'"',collapse=","),'],',
           '"lfc":[', paste(sprintf("%.4f",dv$lfc),collapse=","),'],',
           '"padj":[',paste(sprintf("%.6f",dv$padj),collapse=","),']}')
  },character(1)),collapse=","),"}")
  
  col_colors_js <- paste0("{",paste(vapply(treatments,function(t)
    paste0('"',jesc(t),'":',jstr(treat_colors[t])),character(1)),collapse=","),"}")
  
  js_code <- paste0(
    "var MA=",sprintf("%.6f",max_abs),",SIG_THR=",format(sig_threshold, scientific=FALSE, trim=TRUE),",MC=[",paste0('"',jesc(treatments),'"',collapse=","),
    "],CC=",col_colors_js,";\n",
    "var GD=",go_data_js,";\nvar LD=",leading_js,";\nvar VD=",volcano_js,";\nvar IS_KEGG=",if(type_suffix=="KEGG") "true" else "false",";\n",
    "var TAB20=[",paste0("'",TAB20,"'",collapse=","),"];\n\n",
    
    "document.addEventListener('DOMContentLoaded',function(){\n",
    "  var tip=document.getElementById('tt'),cur=null;\n",
    "  var lockedPathId=null; window.lockedOverlapGenes=null;\n",
    "  var nextIdx=0,lcaCol={},activeNodes={};\n\n",
    
    "  function colorNes(v,sig){\n",
    "    if(v===null||v===undefined)return'#E0E0E0';\n",
    "    if(sig==='ns')return'#D0D0D0';\n",
    "    var t=Math.max(-1,Math.min(1,v/MA)),A=[-1,-0.55,0,0.55,1],C=[[5,48,97],[33,102,172],[247,247,247],[214,96,77],[103,0,31]],k=t<=A[1]?0:(t<=A[2]?1:(t<=A[3]?2:3)),u=(t-A[k])/(A[k+1]-A[k]);\n",
    "    return'rgb('+Math.round(C[k][0]+u*(C[k+1][0]-C[k][0]))+','+Math.round(C[k][1]+u*(C[k+1][1]-C[k][1]))+','+Math.round(C[k][2]+u*(C[k+1][2]-C[k][2]))+')';}\n\n",
    
    "  function showTip(html,x,y){\n",
    "    tip.innerHTML=html;\n",
    "    var btn=document.createElement('span');btn.textContent='×';\n",
    "    btn.style.cssText='float:right;cursor:pointer;font-size:16px;line-height:1';\n",
    "    btn.onclick=function(){tip.style.display='none';};\n",
    "    tip.insertBefore(btn,tip.firstChild);tip.style.display='block';\n",
    "    var lx=x+12;if(lx+tip.offsetWidth+8>window.innerWidth)lx=Math.max(4,x-tip.offsetWidth-8);\n",
    "    var vy=y;if(vy+tip.offsetHeight+8>window.innerHeight)vy=Math.max(4,window.innerHeight-tip.offsetHeight-8);\n",
    "    tip.style.left=lx+'px';tip.style.top=vy+'px';\n  }\n\n",
    
    "  window.highlightPathway=function(pathId){\n",
    "    document.querySelectorAll('#path-panel .path-row').forEach(function(r){r.style.background=r.getAttribute('data-pid')===pathId?'#FFF3CD':'';});\n",
    "    document.querySelectorAll('#gene-panel .gene-tag').forEach(function(el){\n",
    "      var pids=(el.getAttribute('data-paths')||'').split(',').filter(Boolean);\n",
    "      if(pids.indexOf(pathId)>=0){el.style.opacity='1';el.style.background='#FFF3CD';el.style.boxShadow='0 0 0 2px #FFD43B';}\n",
    "      else{el.style.opacity='0.25';el.style.background=el.getAttribute('data-bg')||'';el.style.boxShadow='none';}\n",
    "    });if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);\n  };\n",
    "  window.previewPathway=function(pathId){window.highlightPathway(pathId);};\n",
    "  window.endPathwayPreview=function(){if(lockedPathId)window.highlightPathway(lockedPathId);else window.resetHighlight();if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);};\n",
    "  window.togglePathway=function(pathId){lockedPathId=(lockedPathId===pathId)?null:pathId;if(lockedPathId)window.highlightPathway(lockedPathId);else window.resetHighlight();if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);};\n",
    "  window.resetHighlight=function(){\n",
    "    document.querySelectorAll('#path-panel .path-row').forEach(function(r){r.style.background='';});\n",
    "    document.querySelectorAll('#gene-panel .gene-tag').forEach(function(el){el.style.opacity='';el.style.background=el.getAttribute('data-bg')||'';el.style.boxShadow='none';el.style.fontWeight='400';});if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);\n",
    "  };\n\n",
    
    "  function selectTopLeadingGenes(vd,leSet,n){\n",
    "    n=Math.max(0,parseInt(n||0));if(!n)return{};var arr=[];\n",
    "    for(var i=0;i<vd.orfs.length;i++){var o=vd.orfs[i],p=vd.padj[i],l=vd.lfc[i];\n",
    "      if(!leSet[o]||p===null||l===null||!isFinite(p)||!isFinite(l))continue;\n",
    "      arr.push({i:i,p:p,a:Math.abs(l)});}\n",
    "    arr.sort(function(x,y){if(x.p!==y.p)return x.p-y.p;return y.a-x.a;});\n",
    "    var out={};arr.slice(0,n).forEach(function(x){out[x.i]=1;});return out;\n  }\n\n",
    
    "  function getCommonVolcanoScale(){\n",
    "    var rawL=0,rawN=0;if(!window.currentLeadingData)return{maxL:2,maxN:4};\n",
    "    MC.forEach(function(col){var cd=window.currentLeadingData[col],sg=window.currentGoData?window.currentGoData[col+'_sig']:'';if(!cd||sg!=='*'||!VD[col])return;\n",
    "      VD[col].lfc.forEach(function(v){if(v!==null&&isFinite(v))rawL=Math.max(rawL,Math.abs(v));});\n",
    "      VD[col].padj.forEach(function(p){if(p!==null&&isFinite(p)&&p>0)rawN=Math.max(rawN,-Math.log10(p));});\n",
    "    });\n",
    "    return{maxL:Math.max(2,Math.ceil(rawL/2)*2),maxN:Math.max(4,Math.ceil(rawN))};\n  }\n\n",
    
    "  function buildVolcanoSVG(colName,leGenes,large,topN,scale){\n",
    "    var vd=VD[colName];if(!vd||!vd.orfs.length)return'';scale=scale||getCommonVolcanoScale();\n",
    "    var W=large?430:250,H=large?330:175,pl=large?54:36,pr=large?30:10,pt=large?32:20,pb=large?46:30,pw2=W-pl-pr,ph=H-pt-pb;\n",
    "    var maxL=scale.maxL,maxN=scale.maxN,thr=-Math.log10(SIG_THR);\n",
    "    var xs=function(v){return pl+(v+maxL)/(2*maxL)*pw2;},ys=function(v){return pt+ph*(1-v/maxN);};\n",
    "    var leSet={};(leGenes||[]).forEach(function(g){leSet[g.orf]=g;});\n",
    "    var labels=large?selectTopLeadingGenes(vd,leSet,topN):{};var bgPts=[],lePts=[],labelPts=[];\n",
    "    for(var i=0;i<vd.orfs.length;i++){var o=vd.orfs[i],l=vd.lfc[i],p=vd.padj[i];if(l===null||p===null||!isFinite(l)||!isFinite(p))continue;\n",
    "      var nlp=(p>0)?-Math.log10(p):maxN,isLe=!!leSet[o],fill=isLe?'#D73027':(nlp>=thr?'#555':'#C7C7C7'),r=isLe?(large?3.0:2.35):(large?2.0:1.55);\n",
    "      var handlers=isLe?' onmouseenter=\"highlightGeneFromVolcano(\\\''+o+'\\\',\\\''+colName+'\\\',true)\" onmouseleave=\"highlightGeneFromVolcano(\\\''+o+'\\\',\\\''+colName+'\\\',false)\"':'';\n",
    "      var cx=xs(l),cy=ys(Math.min(nlp,maxN));\n",
    "      var dot='<circle class=\"vlc-dot\" data-orf=\"'+o+'\" data-treatment=\"'+colName+'\" data-le=\"'+(isLe?'1':'0')+'\" data-base-r=\"'+r+'\" data-base-fill=\"'+fill+'\" data-base-stroke=\"'+(isLe?'#D73027':'none')+'\" cx=\"'+cx.toFixed(1)+'\" cy=\"'+cy.toFixed(1)+'\" r=\"'+r+'\" fill=\"'+fill+'\" opacity=\"'+(isLe?'1':'0.55')+'\" stroke=\"'+(isLe?'#D73027':'none')+'\" stroke-width=\"1.4\"'+handlers+'><title>'+((vd.names&&vd.names[i])||o)+' | log2FC='+l.toFixed(2)+' | padj='+p.toExponential(2)+'<\\/title><\\/circle>';if(isLe)lePts.push(dot);else bgPts.push(dot);\n",
    "      if(large&&labels[i]){var nm=(vd.names&&vd.names[i])||o,tx=cx+8,ty=cy-6,anchor='start';if(tx+75>W-pr){tx=cx-8;anchor='end';}\n",
    "        labelPts.push('<line class=\"label-link\" data-orf=\"'+o+'\" x1=\"'+cx.toFixed(1)+'\" y1=\"'+cy.toFixed(1)+'\" x2=\"'+tx.toFixed(1)+'\" y2=\"'+ty.toFixed(1)+'\" stroke=\"#AAAAAA\" stroke-width=\"0.8\" visibility=\"hidden\"/>');\n",
    "        labelPts.push('<text class=\"drag-label\" data-orf=\"'+o+'\" data-point-x=\"'+cx.toFixed(1)+'\" data-point-y=\"'+cy.toFixed(1)+'\" x=\"'+tx.toFixed(1)+'\" y=\"'+ty.toFixed(1)+'\" text-anchor=\"'+anchor+'\" font-size=\"9\" font-family=\"Helvetica,Arial,sans-serif\" font-style=\"italic\" font-weight=\"500\" fill=\"#111111\" style=\"cursor:move;user-select:none\">'+nm+'<\\/text>');}\n",
    "    }\n",
    "    var tyThr=ys(thr).toFixed(1),zx=xs(0).toFixed(1),tickFS=large?8:6.5;var ax='<line x1=\"'+pl+'\" y1=\"'+tyThr+'\" x2=\"'+(W-pr)+'\" y2=\"'+tyThr+'\" stroke=\"#999\" stroke-dasharray=\"3,2\"/>'+\n",
    "      '<line x1=\"'+zx+'\" y1=\"'+pt+'\" x2=\"'+zx+'\" y2=\"'+(pt+ph)+'\" stroke=\"#ddd\"/>'+\n",
    "      '<line x1=\"'+pl+'\" y1=\"'+(pt+ph)+'\" x2=\"'+(W-pr)+'\" y2=\"'+(pt+ph)+'\" stroke=\"#666\"/>'+\n",
    "      '<line x1=\"'+pl+'\" y1=\"'+pt+'\" x2=\"'+pl+'\" y2=\"'+(pt+ph)+'\" stroke=\"#666\"/>';\n",
    "    for(var xv=-maxL;xv<=maxL+0.001;xv+=2){var x=xs(xv);ax+='<line x1=\"'+x.toFixed(1)+'\" y1=\"'+(pt+ph)+'\" x2=\"'+x.toFixed(1)+'\" y2=\"'+(pt+ph+3)+'\" stroke=\"#666\"/><text x=\"'+x.toFixed(1)+'\" y=\"'+(pt+ph+11)+'\" text-anchor=\"middle\" font-size=\"'+tickFS+'\" fill=\"#555\">'+xv.toFixed(0)+'</text>';}\n",
    "    var yt=[0,thr];for(var yv=2;yv<=maxN;yv+=2)yt.push(yv);if(yt.indexOf(maxN)<0)yt.push(maxN);yt.sort(function(a,b){return a-b;}).filter(function(v,i,a){return i===0||Math.abs(v-a[i-1])>0.15;}).forEach(function(v){var y=ys(v);ax+='<line x1=\"'+(pl-3)+'\" y1=\"'+y.toFixed(1)+'\" x2=\"'+pl+'\" y2=\"'+y.toFixed(1)+'\" stroke=\"#666\"/><text x=\"'+(pl-5)+'\" y=\"'+(y+2).toFixed(1)+'\" text-anchor=\"end\" font-size=\"'+tickFS+'\" fill=\"#555\">'+(Math.abs(v-thr)<0.01?(-Math.log10(SIG_THR)).toFixed(2):v.toFixed(0))+'</text>';});\n",
    "    // Clean Y-axis: integer ticks + configured FDR label only\n",
    "    var yStep=maxN<=5?1:2,yTkSvg='',yv=0;\n",
    "    while(yv<=maxN){yTkSvg+='<line x1=\"'+(pl-3)+'\" y1=\"'+ys(yv).toFixed(1)+'\" x2=\"'+pl+'\" y2=\"'+ys(yv).toFixed(1)+'\" stroke=\"#888\" stroke-width=\"0.5\"/>'+\n",
    "      '<text x=\"'+(pl-5)+'\" y=\"'+(parseFloat(ys(yv))+3).toFixed(1)+'\" text-anchor=\"end\" font-size=\"7\" fill=\"#666\">'+yv+'<\\/text>';yv+=yStep;}\n",
    "    ax=ax.replace(/<line[^/]+stroke-dasharray[^/]+\\/>/,'');// remove old threshold line if present\n",
    "    ax+='<line x1=\"'+pl+'\" y1=\"'+ys(thr).toFixed(1)+'\" x2=\"'+(W-pr)+'\" y2=\"'+ys(thr).toFixed(1)+'\" stroke=\"#bbb\" stroke-dasharray=\"3,2\" stroke-width=\"0.8\"/>'+\n",
    "      '<text x=\"'+(W-pr-2)+'\" y=\"'+(parseFloat(ys(thr))-2).toFixed(1)+'\" text-anchor=\"end\" font-size=\"6\" fill=\"#aaa\">FDR='+SIG_THR+'<\\/text>'+\n",
    "      yTkSvg+\n",
    "      '<text x=\"'+(pl+pw2/2)+'\" y=\"'+(H-5)+'\" text-anchor=\"middle\" font-size=\"8\">log\u2082FC<\\/text>'+\n",
    "      '<text x=\"11\" y=\"'+(pt+ph/2)+'\" text-anchor=\"middle\" font-size=\"8\" transform=\"rotate(-90,11,'+(pt+ph/2)+')\">-log\u2081\u2080(adj.P)<\\/text>'+\n",
    "      '<text x=\"'+(W/2)+'\" y=\"17\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"bold\">'+String(colName).replace(/X/g,' × ')+'<\\/text>';\n",
    "    return'<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"'+W+'\" height=\"'+H+'\" onmouseleave=\"resetVolcanoGeneHover()\" style=\"background:white;border:1px solid #ccc;border-radius:4px\">'+bgPts.join('')+lePts.join('')+labelPts.join('')+ax+'<\\/svg>';\n  }\n\n",
    
    "  window.resetVolcanoGeneHover=function(){\n",
    "    document.querySelectorAll('.vlc-dot').forEach(function(d){\n",
    "      var isLe=d.getAttribute('data-le')==='1';\n",
    "      d.setAttribute('r',d.getAttribute('data-base-r')||(isLe?'2.35':'1.55'));\n",
    "      d.setAttribute('stroke',d.getAttribute('data-base-stroke')||(isLe?'#D73027':'none'));\n",
    "      d.setAttribute('stroke-width',isLe?'1.4':'0');d.style.filter='';\n",
    "    });\n",
    "    document.querySelectorAll('#gene-panel .gene-tag').forEach(function(gt){\n",
    "      gt.style.boxShadow='';gt.style.background=gt.getAttribute('data-bg')||'#f0f7ff';\n",
    "    });\n",
    "  };\n\n",
    "  window.highlightGeneFromVolcano=function(orf,col,on){\n",
    "    resetVolcanoGeneHover();if(!on)return;\n",
    "    document.querySelectorAll('.vlc-dot[data-orf=\"'+orf+'\"][data-treatment=\"'+col+'\"][data-le=\"1\"]').forEach(function(d){\n",
    "      var br=parseFloat(d.getAttribute('data-base-r')||'2.35');\n",
    "      d.setAttribute('r',(br*1.7).toFixed(1));d.setAttribute('stroke','#FFE066');d.setAttribute('stroke-width','2.4');\n",
    "      d.style.filter='drop-shadow(0 0 3px #FFF3BF)';\n",
    "    });\n",
    "    document.querySelectorAll('#gene-panel .gene-tag[data-orf=\"'+orf+'\"][data-treatment=\"'+col+'\"]').forEach(function(gt){\n",
    "      gt.style.boxShadow='0 0 0 2px #74C0FC';gt.style.background='#E7F5FF';\n",
    "    });\n",
    "  };\n\n",
    "  function enableLabelDragging(){\n",
    "    document.querySelectorAll('#expanded-plots .drag-label').forEach(function(lbl){var svg=lbl.ownerSVGElement,line=svg.querySelector('.label-link[data-orf=\"'+lbl.getAttribute('data-orf')+'\"]'),drag=false,dx=0,dy=0;\n",
    "      function pt(ev){var p=svg.createSVGPoint();p.x=ev.clientX;p.y=ev.clientY;return p.matrixTransform(svg.getScreenCTM().inverse());}\n",
    "      lbl.addEventListener('pointerdown',function(ev){drag=true;lbl.setPointerCapture(ev.pointerId);var p=pt(ev);dx=parseFloat(lbl.getAttribute('x'))-p.x;dy=parseFloat(lbl.getAttribute('y'))-p.y;ev.preventDefault();});\n",
    "      lbl.addEventListener('pointermove',function(ev){if(!drag)return;var p=pt(ev),x=p.x+dx,y=p.y+dy;lbl.setAttribute('x',x);lbl.setAttribute('y',y);if(line){line.setAttribute('x2',x);line.setAttribute('y2',y);line.setAttribute('visibility','visible');}ev.preventDefault();});\n",
    "      lbl.addEventListener('pointerup',function(ev){drag=false;try{lbl.releasePointerCapture(ev.pointerId);}catch(e){}});lbl.addEventListener('pointercancel',function(){drag=false;});\n",
    "    });\n  }\n\n",
    
    "  window.openExpandedVolcano=function(){\n",
    "    var ov=document.getElementById('volcano-expanded');ov.style.display='flex';document.getElementById('side-panels').style.visibility='hidden';document.getElementById('tt').style.display='none';renderExpandedVolcano();\n  };\n",
    "  window.closeExpandedVolcano=function(){document.getElementById('volcano-expanded').style.display='none';document.getElementById('side-panels').style.visibility='visible';};\n",
    "  window.renderExpandedVolcano=function(){var n=parseInt(document.getElementById('topn-input').value||'10');if(n<0)n=0;var box=document.getElementById('expanded-plots');box.innerHTML='';\n",
    "    if(!window.currentLeadingData)return;var scale=getCommonVolcanoScale();MC.forEach(function(col){var cd=window.currentLeadingData[col],sg=window.currentGoData?window.currentGoData[col+'_sig']:'';if(!cd||sg!=='*')return;box.innerHTML+=buildVolcanoSVG(col,cd.genes,true,n,scale);});enableLabelDragging();};\n\n",
    
    "  window.downloadExpandedJPG=async function(){\n",
    "    var svgs=Array.from(document.querySelectorAll('#expanded-plots svg'));if(!svgs.length)return;\n",
    "    var gap=20,scale=600/96;\n",
    "    var widths=svgs.map(function(s){return parseInt(s.getAttribute('width'));}),heights=svgs.map(function(s){return parseInt(s.getAttribute('height'));});\n",
    "    var cssW=widths.reduce(function(a,b){return a+b;},0)+gap*(svgs.length-1),cssH=Math.max.apply(null,heights);\n",
    "    var canvas=document.createElement('canvas');canvas.width=Math.round(cssW*scale);canvas.height=Math.round(cssH*scale);\n",
    "    var ctx=canvas.getContext('2d');ctx.setTransform(scale,0,0,scale,0,0);ctx.fillStyle='white';ctx.fillRect(0,0,cssW,cssH);var x=0;\n",
    "    for(var i=0;i<svgs.length;i++){var clone=svgs[i].cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');var xml=new XMLSerializer().serializeToString(clone);var blob=new Blob([xml],{type:'image/svg+xml;charset=utf-8'});var url=URL.createObjectURL(blob);\n",
    "      await new Promise(function(resolve,reject){var img=new Image();img.onload=function(){ctx.drawImage(img,x,0,widths[i],heights[i]);x+=widths[i]+gap;URL.revokeObjectURL(url);resolve();};img.onerror=reject;img.src=url;});}\n",
    "    canvas.toBlob(async function(blob){if(!blob)return;var buf=new Uint8Array(await blob.arrayBuffer());\n",
    "      if(buf.length>18&&buf[0]===0xFF&&buf[1]===0xD8&&buf[2]===0xFF&&buf[3]===0xE0&&String.fromCharCode(buf[6],buf[7],buf[8],buf[9])==='JFIF'){buf[13]=1;buf[14]=(600>>8)&255;buf[15]=600&255;buf[16]=(600>>8)&255;buf[17]=600&255;}\n",
    "      var out=new Blob([buf],{type:'image/jpeg'}),a=document.createElement('a');a.href=URL.createObjectURL(out);a.download='leading_edge_volcano_plots_600dpi.jpg';a.click();setTimeout(function(){URL.revokeObjectURL(a.href);},1000);},'image/jpeg',0.98);\n",
    "  };\n\n",
    
    "function displayConditionName(x){return String(x==null?'':x).replace(/X/g,'\u00d7');}\n",
    "function splitConditionName(x){var z=String(x||'').split('X');return {primary:z.length>1?z[0]:'',secondary:z.length>1?z.slice(1).join('\u00d7'):String(x||'')};}\n\n",
    "  window.vennGeneMode='all';\n",
    "  function getSignificantLeadingSets(mode){\n",
    "    mode=mode||window.vennGeneMode||'all';var sets=[];if(!window.currentLeadingData||!window.currentGoData)return sets;\n",
    "    MC.forEach(function(col){var cd=window.currentLeadingData[col],sg=window.currentGoData[col+'_sig'];if(!cd||sg!=='*')return;var genes=(cd.genes||[]).filter(function(g){return mode!=='significant'||(typeof g.padj==='number'&&isFinite(g.padj)&&g.padj<SIG_THR);}).map(function(g){return g.orf;});sets.push({name:col,color:CC[col]||'#777',genes:Array.from(new Set(genes))});});\n",
    "    return sets;\n  }\n\n",
    "  window.setVennGeneMode=function(mode){window.vennGeneMode=mode==='significant'?'significant':'all';renderVennDiagrams();};\n",
    "  function updateVennModeButtons(){var allBtn=document.getElementById('venn-mode-all'),sigBtn=document.getElementById('venn-mode-significant'),isSig=window.vennGeneMode==='significant';if(allBtn){allBtn.style.background=isSig?'white':'#444';allBtn.style.color=isSig?'#444':'white';allBtn.style.borderColor='#444';}if(sigBtn){sigBtn.style.background=isSig?'#FFF5F5':'white';sigBtn.style.color='#C92A2A';sigBtn.style.borderColor='#E03131';sigBtn.style.boxShadow=isSig?'0 0 0 2px rgba(224,49,49,.16)':'none';}}\n\n",
    "  function computeExclusiveIntersections(sets){var universe={};sets.forEach(function(s,i){s.genes.forEach(function(g){if(!universe[g])universe[g]=0;universe[g]|=(1<<i);});});var out={};Object.keys(universe).forEach(function(g){var m=universe[g];if(!out[m])out[m]=[];out[m].push(g);});return out;}\n\n",
    "  window.selectOverlapGenes=function(encoded){var genes=[];try{genes=JSON.parse(decodeURIComponent(encoded));}catch(e){}window.lockedOverlapGenes=genes;document.getElementById('venn-expanded').style.display='none';window.applyOverlapGeneSelection();};\n",
    "  window.applyOverlapGeneSelection=function(preserveOthers){var genes=window.lockedOverlapGenes||[];document.querySelectorAll('#gene-panel .gene-tag').forEach(function(el){var hit=genes.indexOf(el.getAttribute('data-orf'))>=0;if(hit){el.style.opacity='1';el.style.background='#FFF3BF';el.style.boxShadow='0 0 0 2px #FFD43B';el.style.fontWeight='700';}else if(!preserveOthers){el.style.opacity='1';el.style.background=el.getAttribute('data-bg')||'';el.style.boxShadow='';el.style.fontWeight='400';}});};\n",
    "  window.clearOverlapGenes=function(){window.lockedOverlapGenes=null;window.resetHighlight();};\n\n",
    "  function overlapNumber(x,y,genes,label){var enc=encodeURIComponent(JSON.stringify(genes));return '<g class=\"venn-count\" onclick=\"selectOverlapGenes(\\''+enc+'\\')\" style=\"cursor:pointer\"><circle cx=\"'+x+'\" cy=\"'+(y-4)+'\" r=\"13\" fill=\"white\" fill-opacity=\"0.9\" stroke=\"#777\" stroke-width=\"0.6\"/><text x=\"'+x+'\" y=\"'+y+'\" text-anchor=\"middle\" font-size=\"10\" font-weight=\"bold\" fill=\"#111\">'+genes.length+'<\\/text><title>'+label+': '+genes.length+' genes<\\/title><\\/g>';}\n\n",
    "  window.openVennDiagrams=function(){window.vennGeneMode='all';document.getElementById('venn-expanded').style.display='flex';renderVennDiagrams();};\n",
    "  function splitVennConditionName(x){return splitConditionName(x);}\n",
    "  window.renderVennDiagrams=function(){var sets=getSignificantLeadingSets(window.vennGeneMode),ov=document.getElementById('venn-expanded'),box=document.getElementById('venn-content');box.innerHTML='';updateVennModeButtons();if(!sets.length){box.innerHTML='<p style=\"color:#888\">No significantly enriched conditions with leading-edge genes.<\\/p>';ov.style.display='flex';return;}var inter=computeExclusiveIntersections(sets),parts=sets.map(function(s){return splitVennConditionName(s.name);}),primaryNames=parts.map(function(p){return p.primary;}).filter(function(x){return x;}),sharedPrimary=primaryNames.length===sets.length&&primaryNames.every(function(x){return x===primaryNames[0];})?primaryNames[0]:'';\n",
    "    if(sets.length<=3){var W=sets.length===1?360:sets.length===2?520:600,H=sets.length===3?430:390,svg='<svg id=\"venn-export-svg\" width=\"'+W+'\" height=\"'+H+'\" viewBox=\"0 0 '+W+' '+H+'\" xmlns=\"http://www.w3.org/2000/svg\">';\n",
    "      if(sets.length===1){var s=sets[0],g=inter[1]||[],p=parts[0],top=sharedPrimary||'',lab=sharedPrimary?p.secondary:displayConditionName(s.name);if(top)svg+='<text x=\"180\" y=\"28\" text-anchor=\"middle\" font-size=\"15\" font-weight=\"bold\">'+top+'<\\/text>';svg+='<circle cx=\"180\" cy=\"205\" r=\"120\" fill=\"'+s.color+'\" fill-opacity=\"0.30\" stroke=\"'+s.color+'\" stroke-width=\"3\"/><text x=\"180\" y=\"62\" text-anchor=\"middle\" font-size=\"12\" font-weight=\"bold\">'+lab+'<\\/text>'+overlapNumber(180,210,g,displayConditionName(s.name));}\n",
    "      if(sets.length===2){var a=sets[0],b=sets[1],g1=inter[1]||[],g2=inter[2]||[],g12=inter[3]||[],la=sharedPrimary?parts[0].secondary:displayConditionName(a.name),lb=sharedPrimary?parts[1].secondary:displayConditionName(b.name);if(sharedPrimary)svg+='<text x=\"260\" y=\"25\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">'+sharedPrimary+'<\\/text>';svg+='<circle cx=\"215\" cy=\"215\" r=\"135\" fill=\"'+a.color+'\" fill-opacity=\"0.30\" stroke=\"'+a.color+'\" stroke-width=\"3\"/><circle cx=\"305\" cy=\"215\" r=\"135\" fill=\"'+b.color+'\" fill-opacity=\"0.30\" stroke=\"'+b.color+'\" stroke-width=\"3\"/><text x=\"145\" y=\"58\" text-anchor=\"middle\" font-size=\"12\" font-weight=\"bold\">'+la+'<\\/text><text x=\"375\" y=\"58\" text-anchor=\"middle\" font-size=\"12\" font-weight=\"bold\">'+lb+'<\\/text>'+overlapNumber(150,220,g1,displayConditionName(a.name)+' only')+overlapNumber(370,220,g2,displayConditionName(b.name)+' only')+overlapNumber(260,220,g12,displayConditionName(a.name)+' ∩ '+displayConditionName(b.name));}\n",
    "      if(sets.length===3){var a=sets[0],b=sets[1],c=sets[2];svg+='<circle cx=\"245\" cy=\"180\" r=\"125\" fill=\"'+a.color+'\" fill-opacity=\"0.27\" stroke=\"'+a.color+'\" stroke-width=\"3\"/><circle cx=\"355\" cy=\"180\" r=\"125\" fill=\"'+b.color+'\" fill-opacity=\"0.27\" stroke=\"'+b.color+'\" stroke-width=\"3\"/><circle cx=\"300\" cy=\"265\" r=\"125\" fill=\"'+c.color+'\" fill-opacity=\"0.27\" stroke=\"'+c.color+'\" stroke-width=\"3\"/><text x=\"170\" y=\"28\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"bold\">'+displayConditionName(a.name)+'<\\/text><text x=\"430\" y=\"28\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"bold\">'+displayConditionName(b.name)+'<\\/text><text x=\"300\" y=\"418\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"bold\">'+displayConditionName(c.name)+'<\\/text>';var pos={1:[190,165],2:[410,165],4:[300,325],3:[300,125],5:[230,255],6:[370,255],7:[300,215]};Object.keys(pos).forEach(function(k){var m=parseInt(k),names=[];sets.forEach(function(s,i){if(m&(1<<i))names.push(displayConditionName(s.name));});svg+=overlapNumber(pos[k][0],pos[k][1],inter[m]||[],names.join(' ∩ '));});}\n",
    "      svg+='</svg><p style=\"font-size:9px;color:#777;margin:6px 0 0\">Click a number to highlight those genes in all leading-edge lists.</p>';box.innerHTML=svg;}else{var rows=[];Object.keys(inter).forEach(function(k){if(inter[k].length)rows.push({mask:parseInt(k),genes:inter[k]});});rows.sort(function(a,b){return b.genes.length-a.genes.length;});var labels='';sets.forEach(function(s){labels+='<span class=\"upset-col-label\" title=\"'+displayConditionName(s.name)+'\">'+displayConditionName(s.name)+'</span>';});var h='<div class=\"upset-table\"><div class=\"upset-head\"><div class=\"upset-labels\">'+labels+'</div><div></div></div>';rows.forEach(function(r){var enc=encodeURIComponent(JSON.stringify(r.genes)),dots='';sets.forEach(function(s,i){var on=!!(r.mask&(1<<i));dots+='<span title=\"'+displayConditionName(s.name)+'\" class=\"upset-dot\" style=\"background:'+(on?s.color:'#ddd')+';border-color:'+(on?s.color:'#bbb')+'\"></span>';});h+='<div class=\"upset-row\"><div class=\"upset-dots\">'+dots+'</div><button onclick=\"selectOverlapGenes(\\''+enc+'\\')\" class=\"upset-count\">'+r.genes.length+' genes</button></div>';});h+='</div>';box.innerHTML=h;}ov.style.display='flex';};\n",
    "  function vennGeneNameMap(){var out={};if(!window.currentLeadingData)return out;Object.keys(window.currentLeadingData).forEach(function(k){var cd=window.currentLeadingData[k];(cd&&cd.genes||[]).forEach(function(g){out[g.orf]=g.name||g.orf;});});return out;}\n",
    "  function vennExportRows(){var sets=getSignificantLeadingSets(window.vennGeneMode),inter=computeExclusiveIntersections(sets),nm=vennGeneNameMap(),rows=[];Object.keys(inter).forEach(function(k){var mask=parseInt(k),conds=[],parts=[];sets.forEach(function(s,i){if(mask&(1<<i)){conds.push(displayConditionName(s.name));parts.push(splitVennConditionName(s.name));}});var prim=parts.map(function(p){return p.primary;}).filter(function(x){return x;}),sharedPrimary=prim.length&&prim.every(function(x){return x===prim[0];})?prim[0]:'',secondary=parts.map(function(p){return p.secondary;});(inter[k]||[]).forEach(function(g){rows.push({Tested_level:sharedPrimary,Conditions:secondary.join('; '),Gene_ORF:g,Gene_name:nm[g]||g,Mode:window.vennGeneMode==='significant'?'Significant genes only':'All leading-edge genes'});});});return rows;}\n",
    "  function escXml(x){return String(x==null?'':x).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;');}\n",
    "  window.exportVennData=async function(){var svg=document.getElementById('venn-export-svg');if(svg){var clone=svg.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');var W=parseFloat(clone.getAttribute('width'))||600,H=parseFloat(clone.getAttribute('height'))||390,scale=600/96,canvas=document.createElement('canvas');canvas.width=Math.round(W*scale);canvas.height=Math.round(H*scale);var ctx=canvas.getContext('2d');ctx.setTransform(scale,0,0,scale,0,0);ctx.fillStyle='white';ctx.fillRect(0,0,W,H);var xml=new XMLSerializer().serializeToString(clone),blob=new Blob([xml],{type:'image/svg+xml;charset=utf-8'}),url=URL.createObjectURL(blob);await new Promise(function(resolve,reject){var img=new Image();img.onload=function(){ctx.drawImage(img,0,0,W,H);URL.revokeObjectURL(url);resolve();};img.onerror=reject;img.src=url;});canvas.toBlob(async function(jpg){if(!jpg)return;var a=document.createElement('a');a.href=URL.createObjectURL(jpg);a.download='leading_edge_venn.jpg';a.click();},'image/jpeg',0.97);}var rows=vennExportRows(),cols=['Tested_level','Conditions','Gene_ORF','Gene_name','Mode'],tsv=[cols.join('\\t')].concat(rows.map(function(r){return cols.map(function(c){return String(r[c]==null?'':r[c]).replace(/\\t/g,' ').replace(/\\r?\\n/g,' ');}).join('\\t');})).join('\\n'),tb=new Blob([tsv],{type:'text/tab-separated-values;charset=utf-8'}),ta=document.createElement('a');ta.href=URL.createObjectURL(tb);ta.download='leading_edge_overlap_genes.tsv';ta.click();};\n",
    "  window.closeVennDiagrams=function(){document.getElementById('venn-expanded').style.display='none';};\n\n",
    
    "  window.showLeadingEdge=function(go_id){\n",
    "    lockedPathId=null;window.lockedOverlapGenes=null;window.resetHighlight();document.getElementById('venn-expanded').style.display='none';document.getElementById('tt').style.display='none';\n",
    "    if(window.activeLeadingGo){document.querySelectorAll('.lbl-txt[data-goid=\"'+window.activeLeadingGo+'\"]').forEach(function(el){el.style.fontWeight='normal';el.style.fontSize='9px';});}\n",
    "    window.activeLeadingGo=go_id;document.querySelectorAll('.lbl-txt[data-goid=\"'+go_id+'\"]').forEach(function(el){el.style.fontWeight='bold';el.style.fontSize='10.5px';});\n",
    "    var ld=LD[go_id],d=GD[go_id];window.currentLeadingData=ld;window.currentGoData=d;\n",
    "    var gp=document.getElementById('gene-panel'),pp=document.getElementById('path-panel');\n",
    "    var closeBtn='';\n",
    "    if(!ld){gp.innerHTML='<div style=\"display:flex;justify-content:space-between\"><b>No data<\\/b>'+closeBtn+'<\\/div>';gp.style.display='block';document.getElementById('side-panels').style.display='block';if(pp)pp.style.display='none';return;}\n",
    "    var hasSigVenn=MC.some(function(col){return d&&d[col+'_sig']==='*'&&ld[col]&&(ld[col].genes||[]).length>0;});\n",
    "    var vennBtn=hasSigVenn?'<button onclick=\"openVennDiagrams()\" style=\"border:0;background:#6F42C1;color:white;border-radius:3px;padding:3px 7px;cursor:pointer;font-size:8px\">Venn &#8599;<\\/button>':'';\n",
    "    var ghtml='<div style=\"display:flex;justify-content:space-between;align-items:center;margin-bottom:6px\">';\n",
    "    ghtml+='<b style=\"font-size:10px\">Leading edge<\\/b><div style=\"display:flex;align-items:center;gap:7px\">'+vennBtn+closeBtn+'<\\/div><\\/div>';\n",
    "    if(d)ghtml+='<div style=\"font-size:9px;color:#555;margin-bottom:6px\">'+d.desc+'<\\/div>';\n",
    "    ghtml+='<div style=\"font-size:7px;color:#888;margin-bottom:5px\"><span style=\"display:inline-block;width:9px;height:9px;border:2px solid #e03131;border-radius:2px;vertical-align:middle;margin-right:3px\"><\\/span>red border = adj.P&lt;'+SIG_THR+' in that condition<\\/div>';\n",
    "    var pathwayMap={};\n",
    "    MC.forEach(function(col){\n",
    "      var cdata=ld[col]||{genes:[]},termSig=d?d[col+'_sig']:'';var bc=CC[col]||'#888';\n",
    "      ghtml+='<div style=\"margin-bottom:10px\">';\n",
    "      ghtml+='<div style=\"background:'+bc+';color:white;padding:3px 7px;font-size:8px;border-radius:3px;margin-bottom:4px;font-weight:bold\">'+displayConditionName(col)+(termSig==='*'?' \\u2014 '+cdata.genes.length+' leading edge genes':'')+'<\\/div>';\n",
    "      if(termSig!=='*'){ghtml+='<div style=\"background:#E6E6E6;color:#666;padding:9px 7px;font-size:9px;border-radius:3px;text-align:center\">Not significantly enriched<\\/div><\\/div>';return;}\n",
    "      if(cdata.genes.length>0){ghtml+='<div style=\"font-size:8px;line-height:2;word-break:break-word\">';\n",
    "        cdata.genes.forEach(function(g){var safeUrl=g.url.replace(/'/g,\"&#39;\"),pathsStr=(g.paths||[]).join(','),gpadj=typeof g.padj==='number'?g.padj:null;\n",
    "          var bgcol='#f0f7ff';var bw='2px',bcstr=gpadj!==null&&gpadj<SIG_THR?'#e03131':'#d0e8ff';\n",
    "          var tip2=g.orf+(gpadj!==null?' padj='+gpadj.toExponential(2):'')+(g.logfc!==null?' lfc='+g.logfc.toFixed(2):'');\n",
    "          ghtml+='<span class=\"gene-tag\" data-orf=\"'+g.orf+'\" data-treatment=\"'+col+'\" data-paths=\"'+pathsStr+'\" data-bg=\"'+bgcol+'\" onclick=\"window.open(\\\''+safeUrl+'\\\')\" title=\"'+tip2+'\" onmouseover=\"hlV(\\\''+g.orf+'\\\')\" onmouseout=\"hlV(null)\" style=\"cursor:pointer;color:#111;margin:1px 6px 1px 0;display:inline-block;box-sizing:border-box;padding:1px 4px;border-radius:3px;line-height:1.7;border:'+bw+' solid '+bcstr+';background:'+bgcol+'\">'+g.name+'<\\/span>';\n",
    "        });ghtml+='<\\/div>';\n",
    "      }else{ghtml+='<span style=\"color:#aaa;font-size:8px;font-style:italic\">empty<\\/span>';}ghtml+='<\\/div>';\n",
    "      if(!IS_KEGG&&cdata.pathways){Object.keys(cdata.pathways).forEach(function(pid){var pw=cdata.pathways[pid];if(!pathwayMap[pid])pathwayMap[pid]={id:pw.id,name:pw.name,url:pw.url,orfs:[],n:0};pw.geneNames.forEach(function(nm){if(pathwayMap[pid].orfs.indexOf(nm)<0){pathwayMap[pid].orfs.push(nm);pathwayMap[pid].n++;}});});}\n",
    "    });\n",
    "    gp.innerHTML=ghtml;gp.style.display='block';document.getElementById('side-panels').style.display='block';var vp=document.getElementById('volcano-panel'),vhtml='',commonScale=getCommonVolcanoScale(),volcanoReason=[];\n",
    "    MC.forEach(function(col){var cd=ld[col]||{genes:[]},sg=d?d[col+'_sig']:'';if(sg!=='*')return;if(!cd.genes.length){volcanoReason.push(displayConditionName(col)+': empty leading edge');return;}if(!VD[col]){volcanoReason.push(displayConditionName(col)+': DE table not loaded');return;}var deSet={};VD[col].orfs.forEach(function(o){deSet[o]=1;});var nMatch=cd.genes.reduce(function(n,g){return n+(deSet[g.orf]?1:0);},0);if(!nMatch){volcanoReason.push(displayConditionName(col)+': 0/'+cd.genes.length+' leading-edge genes matched to DE IDs');return;}vhtml+=buildVolcanoSVG(col,cd.genes,false,0,commonScale);});\n",
    "    if(vhtml){vp.innerHTML='<div style=\"display:flex;justify-content:space-between;align-items:center;margin-bottom:5px\"><b style=\"font-size:9px\">Volcano plots</b><button onclick=\"openExpandedVolcano()\" style=\"border:0;background:#444;color:white;border-radius:3px;padding:3px 7px;cursor:pointer;font-size:8px\">Expand &#8599;</button></div>'+vhtml+(volcanoReason.length?'<div style=\"font-size:8px;color:#999;margin-top:4px\">'+volcanoReason.join('<br>')+'</div>':'');vp.style.display='block';}\n",
    "    else{vp.innerHTML='<p style=\"color:#888;font-size:9px\">'+(volcanoReason.length?volcanoReason.join('<br>'):'No significantly enriched condition with a non-empty leading edge')+'</p>';vp.style.display='block';}\n",
    "    if(IS_KEGG||!pp){if(pp)pp.style.display='none';return;}var pids=Object.keys(pathwayMap);if(!pids.length){pp.style.display='none';return;}\n",
    "    pids.sort(function(a,b){return pathwayMap[b].n-pathwayMap[a].n;});var phtml='<div style=\"display:flex;justify-content:space-between;align-items:center;margin-bottom:8px\"><b style=\"font-size:10px\">SGD Pathways<\\/b><span style=\"font-size:8px;color:#aaa\">('+pids.length+')<\\/span><\\/div><div class=\"path-grid\">';\n",
    "    pids.forEach(function(pid){var pw=pathwayMap[pid],safeUrl=pw.url.replace(/'/g,\"&#39;\"),safePid=pid.replace(/'/g,\"&#39;\");phtml+='<div class=\"path-row\" data-pid=\"'+safePid+'\"><span class=\"path-name\" onclick=\"event.stopPropagation();window.open(\\\''+safeUrl+'\\\')\">'+pw.name+'<\\/span><span class=\"path-count\" title=\"Hover to preview; click to keep genes highlighted\" onmouseover=\"event.stopPropagation();window.previewPathway(\\\''+safePid+'\\\')\" onmouseout=\"event.stopPropagation();window.endPathwayPreview()\" onclick=\"event.stopPropagation();window.togglePathway(\\\''+safePid+'\\\')\">'+pw.n+' genes<\\/span><\\/div>';});\n",
    "    phtml+='<\\/div>';pp.innerHTML=phtml;pp.style.display='block';\n  };\n\n",
    
    "  window.hlV=function(orf){\n",
    "    // Reset ALL dots first (size + border) — robust against sticky state\n",
    "    document.querySelectorAll('.vlc-dot').forEach(function(el){\n",
    "      var isLe=el.getAttribute('data-le')==='1';\n",
    "      el.setAttribute('r',el.getAttribute('data-base-r')||(isLe?'2.35':'1.55'));\n",
    "      el.setAttribute('stroke',isLe?'#D73027':'none');\n",
    "      el.setAttribute('stroke-width',isLe?'1.4':'0');\n",
    "      el.style.filter='';\n",
    "    });\n",
    "    // Apply yellow highlight + size increase to matching LE dots\n",
    "    if(orf){\n",
    "      document.querySelectorAll('.vlc-dot[data-orf=\"'+orf+'\"]').forEach(function(el){\n",
    "        if(el.getAttribute('data-le')==='1'){\n",
    "          var br=parseFloat(el.getAttribute('data-base-r')||'2.35');\n",
    "          el.setAttribute('r',(br*1.7).toFixed(1));\n",
    "          el.setAttribute('stroke','#FFE066');el.setAttribute('stroke-width','2.4');el.style.filter='drop-shadow(0 0 3px #FFF3BF)';\n",
    "        }\n",
    "      });\n",
    "    }\n",
    "    document.querySelectorAll('#gene-panel .gene-tag').forEach(function(gt){\n",
    "      var match=orf&&gt.getAttribute('data-orf')===orf;\n",
    "      gt.style.boxShadow=match?'0 0 0 2px #74C0FC':'';\n",
    "      gt.style.background='#f0f7ff';\n",
    "    });\n",
    "    if(!orf){\n",
    "      document.querySelectorAll('#path-panel .path-name').forEach(function(el){el.style.fontWeight='500';el.style.color='#1971C2';});if(lockedPathId)window.highlightPathway(lockedPathId);if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);\n",
    "      return;\n",
    "    }\n",
    "    var activePids=[];\n",
    "    document.querySelectorAll('#gene-panel .gene-tag[data-orf=\"'+orf+'\"]').forEach(function(gt){\n",
    "      (gt.getAttribute('data-paths')||'').split(',').filter(Boolean).forEach(function(p){if(activePids.indexOf(p)<0)activePids.push(p);});\n",
    "    });\n",
    "    document.querySelectorAll('#path-panel .path-name').forEach(function(el){\n",
    "      var row=el.closest('.path-row'),hit=row&&activePids.indexOf(row.getAttribute('data-pid'))>=0;\n",
    "      el.style.fontWeight=hit?'bold':'500';el.style.color=hit?'#c92a2a':'#1971C2';\n",
    "    });if(lockedPathId)window.highlightPathway(lockedPathId);if(window.lockedOverlapGenes)window.applyOverlapGeneSelection(true);\n  };\n\n",
    
    "  window.closeLeadingPanels=function(){\n",
    "    lockedPathId=null;window.resetHighlight();\n",
    "    document.getElementById('gene-panel').style.display='none';\n",
    "    var pp=document.getElementById('path-panel');if(pp)pp.style.display='none';\n",
    "    document.getElementById('volcano-panel').style.display='none';document.getElementById('side-panels').style.display='none';\n",
    "    if(window.activeLeadingGo){document.querySelectorAll('.lbl-txt[data-goid=\"'+window.activeLeadingGo+'\"]').forEach(function(el){el.style.fontWeight='normal';el.style.fontSize='9px';});window.activeLeadingGo=null;}\n  };\n\n",
    
    "  document.querySelectorAll('.lbl-txt,.elim-btn').forEach(function(el){\n",
    "    el.addEventListener('click',function(e){\n",
    "      var g=el.getAttribute('data-goid');\n",
    "      if(cur===g&&tip.style.display==='block'){tip.style.display='none';cur=null;return;}\n",
    "      cur=g;e.stopPropagation();\n",
    "      var pd=document.getElementById('pd-'+g);if(!pd)return;\n",
    "      var r=el.getBoundingClientRect();showTip(pd.innerHTML,r.right,r.top);\n",
    "    });\n  });\n\n",
    
    "  function applyNodeColors(){\n",
    "    var lc={};\n",
    "    Object.values(activeNodes).forEach(function(nd){nd.leafSet.forEach(function(g){lc[g]=nd.col+'bb';});});\n",
    "    document.querySelectorAll('rect[id^=lb-]').forEach(function(el){el.setAttribute('fill',lc[el.id.slice(3)]||'#FAFAFA');});\n",
    "    document.querySelectorAll('.node-dot').forEach(function(el){\n",
    "      var k=el.getAttribute('data-lca');el.setAttribute('fill',activeNodes[k]?activeNodes[k].col:'#999');});\n  }\n",
    
    "  document.querySelectorAll('.node-dot').forEach(function(el){\n",
    "    el.addEventListener('click',function(e){\n",
    "      var lca=el.getAttribute('data-lca'),desc=el.getAttribute('data-desc');\n",
    "      var leafArr=(el.getAttribute('data-leaves')||'').split(',').filter(Boolean);\n",
    "      if(activeNodes[lca]){delete activeNodes[lca];applyNodeColors();tip.style.display='none';cur=null;e.stopPropagation();return;}\n",
    "      if(!lcaCol[lca]){lcaCol[lca]=TAB20[nextIdx%TAB20.length];nextIdx++;}\n",
    "      var col=lcaCol[lca];\n",
    "      Object.keys(activeNodes).forEach(function(k){\n",
    "        if(activeNodes[k].leafSet.some(function(g){return new Set(leafArr).has(g);}))delete activeNodes[k];});\n",
    "      activeNodes[lca]={col:col,leafSet:leafArr};applyNodeColors();\n",
    "      var r=el.getBoundingClientRect();\n",
    "      showTip('<b style=\"font-size:11px\">'+(desc||lca)+'<\\/b>',r.right,r.top);\n",
    "      cur='node-'+lca;e.stopPropagation();\n",
    "    });\n  });\n\n",
    
    "  document.addEventListener('click',function(e){\n",
    "    var inside=e.target.closest('#tt')||e.target.closest('#gene-panel')||(document.getElementById('path-panel')&&e.target.closest('#path-panel'))||e.target.closest('#volcano-panel');\n",
    "    if(!inside){tip.style.display='none';cur=null;}\n",
    "  });\n",
    "});\n"
  )
  
  html_out <- paste0(
    "<!DOCTYPE html>\n<html lang='en'>\n<head>\n<meta charset='UTF-8'>\n",
    "<title>",htmlesc(type_suffix)," GSEA</title>\n",
    "<style>\n",
    "body{font-family:Arial,sans-serif;font-size:12px;margin:16px;color:#222}\n",
    ".heatmap-wrap{display:flex;align-items:flex-start;gap:12px}\n",
    "svg{border:1px solid #ddd;border-radius:3px;background:white;flex-shrink:0}\n",
    ".lbl-txt{cursor:pointer}.lbl-txt:hover{font-weight:bold}\n",
    ".elim-btn{cursor:pointer}.elim-btn:hover{opacity:.8}\n",
    ".node-dot{cursor:pointer}.node-dot:hover{stroke:#333;stroke-width:2.5}\n",
    ".gene-tag:hover{background:#E7F5FF!important;box-shadow:0 0 0 2px #74C0FC!important;font-weight:normal;opacity:1!important}\n.drag-label{cursor:move;touch-action:none;font-family:Helvetica,Arial,sans-serif;font-style:italic;font-weight:500}\n",
    ".path-row:hover{background:#f0f0f0!important}\n",
    ".path-name:hover{text-decoration:underline}\n",
    "#side-panels{display:none;position:sticky;top:16px;max-width:860px;max-height:88vh;overflow-y:auto;overflow-x:hidden;padding-right:4px}\n",
    "#side-close-wrap{display:flex;justify-content:flex-end;width:min(340px,40%);margin:0 0 4px auto}\n",
    "#side-close{width:24px;height:24px;border:0;border-radius:3px;background:#C92A2A;color:white;font-size:17px;font-weight:bold;line-height:22px;padding:0;cursor:pointer;box-shadow:0 1px 3px rgba(0,0,0,.22)}#side-close:hover{background:#A51111}\n",
    "#leading-main-row{display:grid;grid-template-columns:minmax(300px,1fr) minmax(270px,340px);gap:12px;align-items:start}\n",
    "#gene-panel,#path-panel,#volcano-panel{display:none;background:white;border:1px solid #aaa;border-radius:5px;padding:10px;box-shadow:2px 2px 8px rgba(0,0,0,.15);font-size:10px;box-sizing:border-box}\n",
    ".path-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:7px;align-items:stretch}\n",
    ".path-row{display:flex;justify-content:space-between;align-items:center;padding:7px 8px;border-radius:5px;cursor:pointer;border:1px solid #d9d9d9;background:#fafafa;gap:8px;min-height:34px;box-sizing:border-box}\n",
    ".path-name{color:#1971C2;cursor:pointer;font-weight:500;flex:1;font-size:8px;line-height:1.25}\n",
    ".path-count{background:#1971C2;color:white;border-radius:10px;padding:1px 7px;font-size:8px;flex-shrink:0;cursor:pointer;transition:filter .12s,box-shadow .12s}.path-count:hover{filter:brightness(1.22);box-shadow:0 0 0 2px #FFD43B}\n",
    "#tt{display:none;position:fixed;background:white;border:1px solid #aaa;padding:12px;",
    "min-width:500px;max-width:900px;z-index:9999;",
    "box-shadow:3px 3px 12px rgba(0,0,0,.2);border-radius:5px;font-size:11px;max-height:660px;overflow-y:auto}\n",
    "#tt table{border-collapse:collapse;width:100%}\n",
    "#venn-expanded{display:none;position:fixed;inset:35px;z-index:10020;background:white;border:1px solid #777;border-radius:7px;box-shadow:0 6px 30px rgba(0,0,0,.35);padding:14px;flex-direction:column;overflow:auto}\n",
    ".venn-count:hover circle{stroke:#FFD43B;stroke-width:2}.upset-table{display:inline-grid;gap:5px;min-width:max-content}.upset-head,.upset-row{display:grid;grid-template-columns:max-content 95px;column-gap:12px;align-items:center}.upset-labels,.upset-dots{display:grid;grid-auto-flow:column;grid-auto-columns:22px;justify-content:start;align-items:end}.upset-labels{height:128px;border-bottom:1px solid #bbb;padding-bottom:4px}.upset-col-label{writing-mode:vertical-rl;transform:rotate(180deg);font-size:8px;line-height:1;text-align:left;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;height:118px;padding-left:4px}.upset-row{border:1px solid #ddd;border-radius:5px;padding:7px;background:#fafafa}.upset-dot{display:block;width:11px;height:11px;border-radius:50%;margin:auto;border:1px solid}.upset-count{border:0;background:#444;color:white;border-radius:10px;padding:3px 9px;cursor:pointer;font-size:9px}.upset-count:hover{background:#6F42C1}\n",
    "#volcano-expanded{display:none;position:fixed;inset:20px;z-index:10000;background:white;border:1px solid #777;border-radius:7px;box-shadow:0 6px 30px rgba(0,0,0,.35);padding:14px;flex-direction:column;overflow:auto}\n",
    "#expanded-plots{display:flex;flex-direction:row;gap:12px;align-items:flex-start;flex-wrap:nowrap;overflow-x:auto}\n",
    "</style>\n<script>\n",js_code,"</script>\n</head>\n<body>\n",
    "<div id='tt'></div>\n",
    "<div id='venn-expanded'><div style='display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:10px;gap:12px'><div><b style='font-size:14px'>Leading-edge overlaps</b><div style='font-size:9px;color:#777'>Exclusive gene intersections among significantly enriched conditions</div><div style='display:flex;gap:7px;flex-wrap:wrap;margin-top:9px'><button id='venn-mode-significant' onclick=\"setVennGeneMode('significant')\" style='border:2px solid #E03131;background:white;color:#C92A2A;border-radius:4px;padding:5px 9px;cursor:pointer;font-size:9px;font-weight:bold'>Show only significant genes</button><button id='venn-mode-all' onclick=\"setVennGeneMode('all')\" style='border:2px solid #444;background:#444;color:white;border-radius:4px;padding:5px 9px;cursor:pointer;font-size:9px;font-weight:bold'>Show all leading-edge genes</button></div></div><div style='display:flex;gap:8px;flex-shrink:0'><button onclick='exportVennData()' style='border:0;background:#0072B2;color:white;border-radius:4px;padding:6px 10px;cursor:pointer'>Export data</button><button onclick='clearOverlapGenes();closeVennDiagrams()' style='border:0;background:#888;color:white;border-radius:4px;padding:6px 10px;cursor:pointer'>Clear selection</button><button onclick='closeVennDiagrams()' style='border:0;background:#444;color:white;border-radius:4px;padding:6px 12px;cursor:pointer'>Close</button></div></div><div id='venn-content' style='margin:auto'></div></div>\n",
    "<div id='volcano-expanded'><div style='display:flex;justify-content:space-between;align-items:center;margin-bottom:10px'>",
    "<label style='font-size:12px'><b>Display top n leading-edge gene names</b> ",
    "<input id='topn-input' type='number' min='0' value='5' style='width:65px' oninput='renderExpandedVolcano()'></label>",
    "<div style='display:flex;gap:8px'><button onclick='downloadExpandedJPG()' style='border:0;background:#0072B2;color:white;border-radius:4px;padding:6px 12px;cursor:pointer'>Download JPG</button><button onclick='closeExpandedVolcano()' style='border:0;background:#444;color:white;border-radius:4px;padding:6px 12px;cursor:pointer'>Close</button></div></div>",
    "<div id='expanded-plots'></div></div>\n",
    "<h2>",htmlesc(if(type_suffix=="KEGG") "KEGG pathway enrichment heatmap" else "GO biological-process enrichment heatmap"),"</h2>\n",
    "<p style='font-size:11px;color:#555;margin:0 0 5px'><b>",
    if(interaction_mode) "Reference condition:" else paste0("Reference ",htmlesc(contrast_factor),":"),
    "</b> ",htmlesc(de_ref_display),
    if(interaction_mode && "interaction_factors" %in% names(pipeline_info))
      paste0(" | <b>Interaction:</b> ",
             htmlesc(pipeline_info[["contrast_factor"]]),
             " \u00d7 ",
             htmlesc(gsub(";", " \u00d7 ", pipeline_info[["interaction_factors"]])))
    else "",
    "</p>\n",
    if ("formula_str" %in% names(pipeline_info) || "formula_expanded" %in% names(pipeline_info))
      paste0("<p style='font-size:11px;color:#555;margin:0 0 7px'><b>Tested model:</b> ",
             htmlesc(if("formula_expanded" %in% names(pipeline_info))
               pipeline_info[["formula_expanded"]]
               else pipeline_info[["formula_str"]]),
             "</p>\n")
    else "",
    "<p style='font-size:11px;color:#666;margin:4px 0 10px'>",
    "<b>",n_go,"</b> terms",
    if(do_similarity_reduction && type_suffix=="GO_BP")
      paste0(" | ",htmlesc(metric_label(similarity_metric))," similarity grouping threshold=",similarity_threshold)
    else "",
    if(filter_mode=="interpretation") " | Interpretation filter: active" else " | Interpretation filter: inactive",
    " | <span style='background:#D0D0D0;padding:0 3px;border-radius:2px;font-size:10px'>Grey cells</span> = not significant vs reference",
    " | <span style='background:#1971C2;color:white;border-radius:3px;padding:1px 5px;font-size:10px'>+N</span> = additional terms grouped with the representative</p>\n",
    "<div class='heatmap-wrap'>\n",
    "<div style='overflow-x:auto'>",
    "<svg width='",svg_w,"' height='",svg_h,"' xmlns='http://www.w3.org/2000/svg'>",
    "<rect width='",svg_w,"' height='",svg_h,"' fill='white'/>",
    grid_svg,"\n",header_svg,"\n",
    "<g>",cells_svg,"</g><g>",stars_svg,"</g><g>",star_row_lines_svg,"</g><g>",labels_svg,"</g>",
    dend_svg,"</svg></div>\n",
    "<div id='side-panels'>\n",
    "<div id='side-close-wrap'><button id='side-close' onclick='closeLeadingPanels()' title='Close'>&times;</button></div>\n",
    "<div id='leading-main-row'>\n",
    "<div id='gene-panel'>",
    "<p style='color:#aaa;font-size:10px;margin:0'>Click a term label<br>",
    "then <b style='color:#2ca02c'>&#9654; Leading edge</b></p></div>\n",
    "<div id='volcano-panel'></div>\n",
    "</div>\n",
    if(type_suffix!="KEGG") "<div id='path-panel' style='margin-top:12px;width:100%'></div>\n" else "",
    "</div></div>\n",
    "<div style='display:none'>",panels_html,"</div>\n",
    "</body></html>\n"
  )
  
  html_path <- file.path(output_dir,paste0(type_suffix,"_heatmap.html"))
  writeLines(html_out,html_path,useBytes=FALSE)
  report_names <- c(if (isTRUE(workbook_written)) basename(xlsx_path) else NULL,
                    basename(png_path), basename(html_path))
  cat(sprintf("  Reports: %s\n", paste(report_names, collapse=", ")))
  
  invisible(NULL)
}


process_type("GO_BP")
process_type("KEGG")
build_variance_partition_report()


# Suppress the implicit return value of the last report task in batch runs.
invisible(NULL)
