options(warn = 1)  # Print each warning immediately so it is preserved in the stage log.
# Run GO and KEGG GSEA

results_dir <- "results/"
output_dir  <- "GSEA_results/"
kegg_cache  <- "kegg_cache/"
go_cache    <- "go_cache/"

rank_by <- "t"

organism_config_file <- "organism_config.R"
if (!file.exists(organism_config_file))
  stop(sprintf("Missing organism configuration file: %s", organism_config_file))
source(organism_config_file, local = FALSE)

min_gs      <- 10
max_gs      <- 500
padj_cut    <- 1.0
padj_method <- "BH"
sig_threshold <- 0.05

run_go_bp <- TRUE
run_kegg  <- TRUE


suppressPackageStartupMessages({
  library(clusterProfiler)
  library(tools)
  library(BiocParallel)
})
if (run_go_bp) {
  suppressPackageStartupMessages({
    library(orgdb, character.only = TRUE)
    library(AnnotationDbi)
    library(GO.db)
  })
}
if (run_kegg)  suppressPackageStartupMessages(library(KEGGREST))

if (!requireNamespace(orgdb, quietly = TRUE)) {
  stop(sprintf("OrgDb package '%s' is not installed. Install it with BiocManager::install('%s').",
               orgdb, orgdb))
}
orgdb_obj <- getExportedValue(orgdb, orgdb)
available_keytypes <- AnnotationDbi::keytypes(orgdb_obj)
available_columns  <- AnnotationDbi::columns(orgdb_obj)
if (!key_type %in% available_keytypes) {
  stop(sprintf("key_type '%s' is not available in %s. Available key types: %s",
               key_type, orgdb, paste(available_keytypes, collapse=", ")))
}
if (!symbol_column %in% available_columns) {
  warning(sprintf("symbol_column '%s' is not available in %s; gene IDs will be used as labels.",
                  symbol_column, orgdb))
  symbol_column <- key_type
}

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(kegg_cache, showWarnings = FALSE, recursive = TRUE)
dir.create(go_cache, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(tempdir())) dir.create(tempdir(), recursive = TRUE)


go_term2gene_file <- file.path(go_cache, paste0(orgdb, "_", key_type, "_GOALL_VALID_BP_term2gene.tsv"))
go_term2name_file <- file.path(go_cache, paste0(orgdb, "_GOALL_VALID_BP_term2name.tsv"))

load_or_build_go_bp <- function() {
  if (file.exists(go_term2gene_file) && file.exists(go_term2name_file)) {
    t2g <- tryCatch(read.table(go_term2gene_file, header=TRUE, sep="\t",
                               stringsAsFactors=FALSE, check.names=FALSE,
                               quote="", comment.char=""),
                    error=function(e) NULL)
    t2n <- tryCatch(read.table(go_term2name_file, header=TRUE, sep="\t",
                               stringsAsFactors=FALSE, check.names=FALSE,
                               quote="", comment.char=""),
                    error=function(e) NULL)
    if (!is.null(t2g) && !is.null(t2n) && nrow(t2g) > 0 && nrow(t2n) > 0) {
      cat(sprintf("GO BP annotations: cache (%d terms, %d gene-term links)\n",
                  nrow(t2n), nrow(t2g)))
      return(list(term2gene=t2g, term2name=t2n))
    }
    file.remove(go_term2gene_file, go_term2name_file)
  }

  cat("GO BP annotations: building cache from OrgDb...\n")
  all_keys <- AnnotationDbi::keys(orgdb_obj, keytype=key_type)
  go_map <- suppressMessages(AnnotationDbi::select(
    orgdb_obj, keys=all_keys, columns=c("GOALL", "ONTOLOGYALL"), keytype=key_type
  ))
  go_map <- go_map[!is.na(go_map$GOALL) & go_map$ONTOLOGYALL == "BP", , drop=FALSE]
  if (!nrow(go_map)) stop(sprintf("No GO BP mappings found for key type '%s'.", key_type))

  t2g <- unique(data.frame(
    TERM=as.character(go_map$GOALL),
    GENE=as.character(go_map[[key_type]]),
    stringsAsFactors=FALSE
  ))
  t2g <- t2g[nzchar(t2g$TERM) & nzchar(t2g$GENE), , drop=FALSE]

  go_ids <- sort(unique(t2g$TERM))

  valid_go_ids <- intersect(
    go_ids,
    AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
  )
  obsolete_go_ids <- setdiff(go_ids, valid_go_ids)
  if (length(obsolete_go_ids) > 0L) {
    cat(sprintf(
      "  Removed %d obsolete/unresolved GO ID(s) before GSEA.\n",
      length(obsolete_go_ids)
    ))
  }
  if (!length(valid_go_ids)) {
    stop("No valid GO BP identifiers remained after validation against GO.db.")
  }

  t2g <- t2g[t2g$TERM %in% valid_go_ids, , drop = FALSE]
  go_ids <- sort(unique(t2g$TERM))
  valid_go_ids <- go_ids
  go_names <- setNames(go_ids, go_ids)
  if (length(valid_go_ids) > 0L) {
    go_desc <- suppressMessages(
      AnnotationDbi::select(
        GO.db::GO.db,
        keys = valid_go_ids,
        columns = "TERM",
        keytype = "GOID"
      )
    )
    go_desc <- go_desc[
      !is.na(go_desc$GOID) & !is.na(go_desc$TERM) & nzchar(go_desc$TERM),
      c("GOID", "TERM"),
      drop = FALSE
    ]
    go_desc <- go_desc[!duplicated(go_desc$GOID), , drop = FALSE]
    go_names[go_desc$GOID] <- go_desc$TERM
  }

  t2n <- data.frame(
    TERM = go_ids,
    NAME = unname(go_names[go_ids]),
    stringsAsFactors = FALSE
  )

  write.table(t2g, go_term2gene_file, sep="\t", row.names=FALSE, quote=FALSE)
  write.table(t2n, go_term2name_file, sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("GO BP annotations: cache saved (%d terms, %d gene-term links)\n",
              nrow(t2n), nrow(t2g)))
  list(term2gene=t2g, term2name=t2n)
}

go_annot <- NULL
if (run_go_bp) {
  go_annot <- tryCatch(
    load_or_build_go_bp(),
    error=function(e) {
      cat(sprintf("  ERROR: GO BP annotations unavailable — %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(go_annot)) stop("GO BP was requested but its annotation data could not be loaded.", call.=FALSE)
}

kegg_term2gene_file <- file.path(kegg_cache, paste0(kegg_org, "_mapping_v2_term2gene.tsv"))
kegg_term2name_file <- file.path(kegg_cache, paste0(kegg_org, "_mapping_v2_term2name.tsv"))

map_kegg_ids_to_input <- function(raw_ids) {
  raw_ids <- unique(trimws(sub("^genes/", "", sub(paste0("^", kegg_org, ":"), "", as.character(raw_ids)))))
  raw_ids <- raw_ids[!is.na(raw_ids) & nzchar(raw_ids)]
  if (!length(raw_ids)) stop("KEGG did not return any usable gene IDs.", call.=FALSE)
  available <- AnnotationDbi::keytypes(orgdb_obj)
  candidates <- intersect(c(key_type, "ENTREZID", "ENSEMBL", "REFSEQ", "SYMBOL", "ORF", "SGD"), available)
  # Convert each plausible KEGG identifier namespace; preserve direct IDs and
  # one-to-many conversions instead of relying on an arbitrary 50% shortcut.
  target_keys <- AnnotationDbi::keys(orgdb_obj, keytype=key_type)
  direct <- intersect(raw_ids, target_keys)
  maps <- list()
  if (length(direct)) maps[["direct"]] <- data.frame(RAW=direct, GENE=direct, stringsAsFactors=FALSE)
  source_stats <- character()
  for (kt in setdiff(candidates, key_type)) {
    source_keys <- tryCatch(AnnotationDbi::keys(orgdb_obj, keytype=kt), error=function(e) character())
    hit <- intersect(raw_ids, source_keys)
    if (!length(hit)) next
    converted <- tryCatch(suppressMessages(AnnotationDbi::select(orgdb_obj, keys=hit,
                                              keytype=kt, columns=key_type)), error=function(e) NULL)
    if (is.null(converted) || !all(c(kt, key_type) %in% names(converted))) next
    converted <- unique(data.frame(RAW=as.character(converted[[kt]]),
                         GENE=as.character(converted[[key_type]]), stringsAsFactors=FALSE))
    converted <- converted[!is.na(converted$RAW) & !is.na(converted$GENE) &
                           nzchar(converted$RAW) & nzchar(converted$GENE) &
                           converted$GENE %in% target_keys, , drop=FALSE]
    if (nrow(converted)) {
      maps[[kt]] <- converted
      source_stats <- c(source_stats, sprintf("%s=%d", kt, length(unique(converted$RAW))))
    }
  }
  if (!length(maps)) stop(sprintf("Unable to map KEGG gene IDs to %s (%s).",key_type,orgdb),call.=FALSE)
  # Use one source namespace consistently. Merging every plausible key type
  # can introduce false mappings when the same raw string denotes distinct genes.
  n_by_source <- vapply(maps, function(x) length(unique(x$RAW)), integer(1))
  best_source <- names(n_by_source)[which.max(n_by_source)]
  mapped <- unique(maps[[best_source]])
  # A direct target match is unambiguous in the requested identifier space;
  # retain it for raw IDs not recognized by the selected source namespace.
  if (best_source != "direct" && "direct" %in% names(maps)) {
    extra <- maps[["direct"]][!maps[["direct"]]$RAW %in% mapped$RAW,,drop=FALSE]
    mapped <- unique(rbind(mapped,extra))
  }
  ambiguous <- sum(table(mapped$RAW) > 1L)
  cat(sprintf("  KEGG IDs mapped to %s: %d/%d raw IDs (source=%s; direct=%d; multiple targets=%d)\n",
       key_type, length(unique(mapped$RAW)),length(raw_ids),best_source,length(direct),ambiguous))
  if (length(unique(mapped$RAW))/length(raw_ids) < 0.50)
    warning("Less than 50% of KEGG raw gene IDs could be mapped; check organism and identifier settings.",call.=FALSE)
  mapped
}

load_or_download_kegg <- function() {
  
  if (file.exists(kegg_term2gene_file) && file.exists(kegg_term2name_file)) {
    t2g <- tryCatch(
      read.table(kegg_term2gene_file, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "", comment.char = ""),
      error = function(e) NULL
    )
    t2n <- tryCatch(
      read.table(kegg_term2name_file, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "", comment.char = ""),
      error = function(e) NULL
    )
    if (!is.null(t2g) && !is.null(t2n) && nrow(t2g) > 0 && nrow(t2n) > 0) {
      # Normalize KEGG IDs consistently in TERM2GENE and TERM2NAME.
      t2g$TERM <- sub("^path:", "", as.character(t2g$TERM))
      t2n$TERM <- sub("^path:", "", as.character(t2n$TERM))
      term_overlap <- intersect(unique(t2g$TERM), unique(t2n$TERM))
      if (length(term_overlap) > 0L) {
        t2g <- t2g[t2g$TERM %in% term_overlap, , drop=FALSE]
        t2n <- t2n[t2n$TERM %in% term_overlap, , drop=FALSE]
        cat(sprintf("KEGG annotations: cache (%d pathways, %d gene-pathway links)\n",
                    nrow(t2n), nrow(t2g)))
        return(list(term2gene = t2g, term2name = t2n))
      }
      cat("  KEGG cache uses incompatible pathway identifiers — rebuilding\n")
    }
    cat("  Cache incomplete or empty — re-downloading\n")
    file.remove(kegg_term2gene_file)
    file.remove(kegg_term2name_file)
  }
  
  cat("KEGG annotations: downloading and building cache...\n")
  
  pw_list <- keggList("pathway", kegg_org)
  t2n <- data.frame(
    TERM = sub("^path:", "", names(pw_list)),
    NAME = sub("^(.*) - ([^-]+)$", "\\1", unname(pw_list), perl=TRUE),
    stringsAsFactors = FALSE
  )
  
  all_links <- keggLink(kegg_org, "pathway")
  raw_links <- data.frame(
    TERM = sub("^path:", "", names(all_links)),
    RAW  = sub(paste0("^", kegg_org, ":"), "", unname(all_links)),
    stringsAsFactors = FALSE
  )
  id_map <- map_kegg_ids_to_input(raw_links$RAW)
  t2g <- merge(raw_links, id_map, by="RAW", all=FALSE)
  t2g <- unique(t2g[, c("TERM", "GENE"), drop=FALSE])
  
  write.table(t2g, kegg_term2gene_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(t2n, kegg_term2name_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("KEGG annotations: cache saved (%d pathways, %d gene-pathway links)\n",
              nrow(t2n), nrow(t2g)))
  
  list(term2gene = t2g, term2name = t2n)
}

kegg_annot <- NULL
if (run_kegg) {
  kegg_annot <- tryCatch(
    load_or_download_kegg(),
    error = function(e) {
      cat(sprintf("  ERROR: KEGG annotations unavailable — %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(kegg_annot)) stop("KEGG was requested but its annotation data could not be loaded.", call.=FALSE)
}

organism_manifest <- data.frame(
  key=c("organism_name","orgdb","key_type","symbol_column","kegg_org"),
  value=c(organism_name,orgdb,key_type,symbol_column,kegg_org),
  stringsAsFactors=FALSE
)
write.table(organism_manifest, file.path(output_dir,"organism_manifest.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

normalize_input_gene_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  if (toupper(key_type) == "ENSEMBL") {
    stripped <- sub("\\.[0-9]+$", "", ids, perl=TRUE)
    attr(stripped, "n_version_suffix_stripped") <- sum(stripped != ids, na.rm=TRUE)
    return(stripped)
  }
  attr(ids, "n_version_suffix_stripped") <- 0L
  ids
}

.duplicate_resolution_rows <- list()

make_ranked <- function(filepath, rank_col) {
  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE,
                   quote = "", comment.char = "")
  if (!"gene" %in% colnames(df))
    stop(sprintf("ERROR: No 'gene' column in %s.", basename(filepath)))
  effective_rank_col <- if (identical(rank_col, "t") && "z.std" %in% colnames(df)) "z.std" else rank_col
  if (!effective_rank_col %in% colnames(df))
    stop(sprintf("ERROR: Column '%s' not found in %s. Available: %s",
                 effective_rank_col, basename(filepath),
                 paste(colnames(df), collapse = ", ")))

  # Annotation and GSEA use original gene IDs, not display suffixes.
  annotation_ids <- if ("original_gene" %in% colnames(df)) df$original_gene else df$gene
  norm_ids <- normalize_input_gene_ids(annotation_ids)
  n_stripped <- attr(norm_ids, "n_version_suffix_stripped")
  rank_values <- suppressWarnings(as.numeric(df[[effective_rank_col]]))
  input_ids <- as.character(df$gene)
  ok <- !is.na(rank_values) & !is.na(norm_ids) & nzchar(norm_ids)
  norm_ids <- norm_ids[ok]
  rank_values <- rank_values[ok]
  input_ids <- input_ids[ok]

  if (anyDuplicated(norm_ids)) {
    # Collapse repeated annotation IDs to one mean ranking statistic and record the resolution.
    groups <- split(seq_along(norm_ids), norm_ids)
    dup_groups <- groups[lengths(groups) > 1L]
    if (length(dup_groups)) {
      contrast_name <- file_path_sans_ext(basename(filepath))
      .duplicate_resolution_rows[[contrast_name]] <<- do.call(rbind, lapply(names(dup_groups), function(id) {
        ii <- dup_groups[[id]]
        data.frame(
          contrast=contrast_name,
          annotation_gene=id,
          n_input_features=length(ii),
          input_feature_ids=paste(input_ids[ii], collapse=" | "),
          rank_metric=effective_rank_col,
          aggregation="arithmetic mean",
          stringsAsFactors=FALSE,
          check.names=FALSE
        )
      }))
      warning(sprintf(
        "Collapsed %d duplicated annotation identifier(s) in %s for GSEA by arithmetic mean of the ranking statistic; individual input rows remain separate in DE/volcano output. See duplicate_gene_resolution.tsv.",
        length(dup_groups), basename(filepath)
      ), call.=FALSE)
    }
    agg_ids <- names(groups)
    agg_values <- vapply(groups, function(ii) mean(rank_values[ii]), numeric(1))
    norm_ids <- agg_ids
    rank_values <- unname(agg_values[agg_ids])
  }
  ranked <- setNames(rank_values, norm_ids)
  ranked <- sort(ranked, decreasing = TRUE)
  attr(ranked, "rank_metric") <- effective_rank_col
  attr(ranked, "n_version_suffix_stripped") <- n_stripped
  ranked
}

run_gsea_go <- function(ranked, label, t2g, t2n, n_retry=3L, sleep_sec=2) {
  last_err <- NULL
  for (attempt in seq_len(n_retry)) {
    if (attempt > 1L) {
      cat(sprintf("    Retry %d/%d after %.0f s...\n", attempt, n_retry, sleep_sec))
      Sys.sleep(sleep_sec)
      sleep_sec <- sleep_sec * 1.5
      gc(verbose=FALSE)
    }
    BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
    res <- tryCatch({
      GSEA(
        geneList=ranked, TERM2GENE=t2g, TERM2NAME=t2n,
        exponent=1, minGSSize=min_gs, maxGSSize=max_gs,
        eps=1e-10, pvalueCutoff=padj_cut, pAdjustMethod=padj_method,
        verbose=FALSE, seed=TRUE, by="fgsea"
      )
    }, error=function(e) {
      last_err <<- conditionMessage(e)
      NULL
    })
    if (!is.null(res)) return(res)
  }
  stop(sprintf("GO BP GSEA failed for contrast '%s' after %d attempts: %s", label, n_retry, last_err), call.=FALSE)
}

run_gsea_kegg <- function(ranked, label, t2g, t2n, n_retry = 3L, sleep_sec = 2) {
  last_err <- NULL
  for (attempt in seq_len(n_retry)) {
    if (attempt > 1L) {
      cat(sprintf("    Retry %d/%d after %.0f s...\n", attempt, n_retry, sleep_sec))
      Sys.sleep(sleep_sec)
      sleep_sec <- sleep_sec * 1.5
      gc(verbose = FALSE)
    }
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    res <- tryCatch({
      GSEA(
        geneList      = ranked,
        TERM2GENE     = t2g,
        TERM2NAME     = t2n,
        exponent      = 1,
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        eps           = 1e-10,
        pvalueCutoff  = padj_cut,
        pAdjustMethod = padj_method,
        verbose       = FALSE,
        seed          = TRUE,
        by            = "fgsea"
      )
    }, error = function(e) {
      last_err <<- conditionMessage(e)
      NULL
    })
    if (!is.null(res)) return(res)
  }
  stop(sprintf("KEGG GSEA failed for contrast '%s' after %d attempts: %s", label, n_retry, last_err), call.=FALSE)
}

contrast_files <- list.files(results_dir, pattern = "\\.txt$", full.names = TRUE)
if (!length(contrast_files))
  stop(sprintf("ERROR: No .txt files found in '%s'.", results_dir))

cat(sprintf("Contrasts to analyze: %d\n", length(contrast_files)))

# Record annotation-ID coverage after ranked-list normalization.
coverage_df <- tryCatch({
  probe <- read.table(contrast_files[[1]], header=TRUE, sep="\t", stringsAsFactors=FALSE,
                      check.names=FALSE, quote="", comment.char="")
  if (!"gene" %in% colnames(probe)) stop("No gene column")
  raw_probe_ids <- trimws(as.character(if ("original_gene" %in% colnames(probe)) probe$original_gene else probe$gene))
  probe_ids <- normalize_input_gene_ids(raw_probe_ids)
  n_stripped <- attr(probe_ids, "n_version_suffix_stripped")
  probe_ids <- unique(probe_ids[!is.na(probe_ids) & nzchar(probe_ids)])
  rows <- list()
  if (run_go_bp && !is.null(go_annot)) {
    mapped <- intersect(probe_ids, unique(as.character(go_annot$term2gene$GENE)))
    rows[["GO_BP"]] <- data.frame(database="GO_BP", input_unique_genes=length(probe_ids),
      mapped_unique_genes=length(mapped), mapped_fraction=if(length(probe_ids)) length(mapped)/length(probe_ids) else NA_real_,
      ensembl_version_suffixes_stripped=n_stripped, stringsAsFactors=FALSE)
  }
  if (run_kegg && !is.null(kegg_annot)) {
    mapped <- intersect(probe_ids, unique(as.character(kegg_annot$term2gene$GENE)))
    rows[["KEGG"]] <- data.frame(database="KEGG", input_unique_genes=length(probe_ids),
      mapped_unique_genes=length(mapped), mapped_fraction=if(length(probe_ids)) length(mapped)/length(probe_ids) else NA_real_,
      ensembl_version_suffixes_stripped=n_stripped, stringsAsFactors=FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}, error=function(e) {
  warning(sprintf("Could not calculate gene-ID mapping coverage: %s", conditionMessage(e)), call.=FALSE)
  NULL
})
if (!is.null(coverage_df)) {
  rownames(coverage_df) <- NULL
  write.table(coverage_df, file.path(output_dir, "gene_id_mapping_coverage.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  for (ii in seq_len(nrow(coverage_df))) {
    cat(sprintf("  %s mapping coverage: %d/%d genes (%.1f%%)%s\n",
                coverage_df$database[ii], coverage_df$mapped_unique_genes[ii],
                coverage_df$input_unique_genes[ii], 100*coverage_df$mapped_fraction[ii],
                if (coverage_df$ensembl_version_suffixes_stripped[ii] > 0)
                  sprintf("; stripped Ensembl version suffixes from %d input IDs", coverage_df$ensembl_version_suffixes_stripped[ii]) else ""))
  }
  cat("\n")
}

summary_rows <- list()
status_file <- file.path(output_dir, "GSEA_status.tsv")
write_status <- function() {
  if (!length(summary_rows)) return(invisible(NULL))
  tab <- do.call(rbind, summary_rows)
  write.table(tab, status_file, sep="\t", row.names=FALSE, quote=FALSE)
  invisible(NULL)
}

for (contrast_index in seq_along(contrast_files)) {
  filepath <- contrast_files[[contrast_index]]
  contrast_name <- file_path_sans_ext(basename(filepath))
  
  ranked <- tryCatch(
    make_ranked(filepath, rank_by),
    error = function(e) stop(sprintf("Cannot build ranking for contrast '%s': %s", contrast_name, conditionMessage(e)), call.=FALSE)
  )
  rank_metric_used <- attr(ranked, "rank_metric")
  
  out_prefix <- file.path(output_dir, contrast_name)
  n_go_sig   <- NA_integer_
  n_kegg_sig <- NA_integer_
  go_status <- if (run_go_bp) "PENDING" else "SKIPPED"
  kegg_status <- if (run_kegg) "PENDING" else "SKIPPED"
  go_message <- ""
  kegg_message <- ""
  
  if (run_go_bp && !is.null(go_annot)) {
    res_go <- tryCatch(
      run_gsea_go(ranked, contrast_name, go_annot$term2gene, go_annot$term2name),
      error=function(e) {
        go_status <<- "FAILED"; go_message <<- conditionMessage(e)
        NULL
      }
    )
    if (identical(go_status, "FAILED")) {
      summary_rows[[contrast_name]] <- data.frame(
        contrast=contrast_name, rank_metric=rank_metric_used, genes_ranked=length(ranked),
        GO_BP_status=go_status, GO_BP_sig=NA_integer_, GO_BP_message=go_message,
        KEGG_status=kegg_status, KEGG_sig=n_kegg_sig, KEGG_message=kegg_message,
        stringsAsFactors=FALSE
      )
      write_status()
      stop(go_message, call.=FALSE)
    }
    go_df <- as.data.frame(res_go)
    go_status <- if (nrow(go_df)) "SUCCESS" else "SUCCESS_EMPTY"
    n_go_sig <- if (nrow(go_df)) sum(go_df$p.adjust < sig_threshold, na.rm=TRUE) else 0L
    tsv_path <- paste0(out_prefix, "_GO_BP.tsv")
    write.table(go_df, tsv_path, sep="\t", row.names=FALSE, quote=FALSE)
  }
  
  if (run_kegg && !is.null(kegg_annot)) {
    res_kegg <- tryCatch(
      run_gsea_kegg(ranked, contrast_name, kegg_annot$term2gene, kegg_annot$term2name),
      error=function(e) {
        kegg_status <<- "FAILED"; kegg_message <<- conditionMessage(e)
        NULL
      }
    )
    if (identical(kegg_status, "FAILED")) {
      summary_rows[[contrast_name]] <- data.frame(
        contrast=contrast_name, rank_metric=rank_metric_used, genes_ranked=length(ranked),
        GO_BP_status=go_status, GO_BP_sig=n_go_sig, GO_BP_message=go_message,
        KEGG_status=kegg_status, KEGG_sig=NA_integer_, KEGG_message=kegg_message,
        stringsAsFactors=FALSE
      )
      write_status()
      stop(kegg_message, call.=FALSE)
    }
    kegg_df <- as.data.frame(res_kegg)
    kegg_status <- if (nrow(kegg_df)) "SUCCESS" else "SUCCESS_EMPTY"
    n_kegg_sig <- if (nrow(kegg_df)) sum(kegg_df$p.adjust < sig_threshold, na.rm=TRUE) else 0L
    tsv_path <- paste0(out_prefix, "_KEGG.tsv")
    write.table(kegg_df, tsv_path, sep="\t", row.names=FALSE, quote=FALSE)
  }
  
  Sys.sleep(0.5)
  gc(verbose = FALSE)
  
  cat(sprintf("  [%d/%d] %s | genes=%d | GO=%s | KEGG=%s\n",
              contrast_index, length(contrast_files), contrast_name, length(ranked),
              ifelse(is.na(n_go_sig), "-", as.character(n_go_sig)),
              ifelse(is.na(n_kegg_sig), "-", as.character(n_kegg_sig))))

  summary_rows[[contrast_name]] <- data.frame(
    contrast=contrast_name,
    rank_metric=rank_metric_used,
    genes_ranked=length(ranked),
    GO_BP_status=go_status,
    GO_BP_sig=n_go_sig,
    GO_BP_message=go_message,
    KEGG_status=kegg_status,
    KEGG_sig=n_kegg_sig,
    KEGG_message=kegg_message,
    stringsAsFactors=FALSE
  )
  write_status()
}

cat(sprintf("GSEA complete: %d contrast(s) | rank=%s | p.adjust<%.4g\n",
            length(summary_rows), paste(unique(vapply(summary_rows, function(x) x$rank_metric[[1]], character(1))), collapse=","), sig_threshold))

if (length(.duplicate_resolution_rows)) {
  dup_tab <- do.call(rbind, .duplicate_resolution_rows)
  rownames(dup_tab) <- NULL
  write.table(dup_tab, file.path(output_dir, "duplicate_gene_resolution.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}
