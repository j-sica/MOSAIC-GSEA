options(warn = 1)  # Print each warning immediately so it is preserved in the stage log.
# Build summary tables

gsea_dir        <- "GSEA_results/"
results_dir     <- "results/"
sig_threshold   <- 0.05


suppressPackageStartupMessages(library(multcompView))

load_tsv <- function(path) {
  df <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      check.names = FALSE,
      fill = TRUE
    ),
    error = function(e) {
      stop(sprintf("Cannot read GSEA file '%s': %s", basename(path),
                   conditionMessage(e)), call. = FALSE)
    }
  )

  if (!"ID" %in% colnames(df)) {
    stop(sprintf("GSEA file '%s' has no ID column. Columns: %s",
                 basename(path), paste(colnames(df), collapse = ", ")),
         call. = FALSE)
  }

  df$ID <- trimws(as.character(df$ID))
  df <- df[!is.na(df$ID) & nzchar(df$ID), , drop = FALSE]
  df <- df[!duplicated(df$ID), , drop = FALSE]
  rownames(df) <- df$ID
  df
}

manifest_path <- file.path(results_dir, "analysis_manifest.tsv")
contrast_manifest_path <- file.path(results_dir, "contrast_manifest.tsv")

if (!file.exists(manifest_path)) {
  cat("WARNING: analysis_manifest.tsv not found — using automatic filename detection\n")
  pipeline_info    <- setNames(character(0), character(0))
  pipeline_contrasts <- data.frame()
} else {
  mi <- read.table(manifest_path, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                   quote="", comment.char="", check.names=FALSE)
  pipeline_info <- setNames(mi$value, mi$key)
  pipeline_contrasts <- if (file.exists(contrast_manifest_path))
    read.table(contrast_manifest_path, header=TRUE, sep="\t",
               stringsAsFactors=FALSE, quote="", comment.char="",
               check.names=FALSE)
  else data.frame()
}

analysis_mode <- if ("analysis_mode" %in% names(pipeline_info))
  pipeline_info[["analysis_mode"]] else "standard"
n_int_factors <- if ("n_interaction_factors" %in% names(pipeline_info))
  as.integer(pipeline_info[["n_interaction_factors"]]) else 0L
contrast_factor_manifest <- if ("contrast_factor" %in% names(pipeline_info))
  pipeline_info[["contrast_factor"]] else "treatment"
reference_level_manifest <- if ("reference_level" %in% names(pipeline_info))
  pipeline_info[["reference_level"]] else NA_character_

cat(sprintf("Mode: %s\n", analysis_mode))

compute_cld <- function(treatments, get_stats) {
  n <- length(treatments)
  if (n == 0) return(setNames(character(0), character(0)))
  if (n == 1) return(setNames("a", treatments))
  safe  <- setNames(paste0("T", seq_along(treatments)), treatments)
  pairs <- combn(treatments, 2, simplify=FALSE)
  stats <- lapply(pairs, function(p) get_stats(p[1], p[2]))
  raw_p <- vapply(stats, function(z) if (is.null(z) || is.na(z$raw)) 1 else as.numeric(z$raw), numeric(1))
  bh_q  <- vapply(stats, function(z) if (is.null(z) || is.na(z$bh)) 1 else as.numeric(z$bh), numeric(1))
  holm_p <- p.adjust(raw_p, method="holm")
  decision_p <- ifelse(bh_q < sig_threshold & holm_p < sig_threshold, 0, 1)
  names(decision_p) <- vapply(pairs, function(p) paste(safe[p[1]], safe[p[2]], sep="-"), character(1))
  cld_obj <- tryCatch(multcompLetters(decision_p, threshold=sig_threshold), error=function(e) NULL)
  if (is.null(cld_obj)) return(setNames(rep("a", n), treatments))
  raw <- cld_obj$Letters
  setNames(unname(raw[safe[treatments]]), treatments)
}

get_gsea_stats <- function(df, gid) {
  if (is.null(df) || !gid %in% rownames(df)) return(list(raw=NA_real_, bh=NA_real_))
  raw_col <- intersect(c("pvalue", "P.Value", "p.value"), colnames(df))[1]
  bh_col  <- intersect(c("p.adjust", "padj", "adj.P.Val"), colnames(df))[1]
  list(
    raw = if (is.na(raw_col)) NA_real_ else as.numeric(df[gid, raw_col]),
    bh  = if (is.na(bh_col)) NA_real_ else as.numeric(df[gid, bh_col])
  )
}

split_condition_components <- function(x) {
  if (grepl("×", x, fixed=TRUE)) strsplit(x, "×", fixed=TRUE)[[1]]
  else as.character(x)
}

refresh_gsea_status_counts <- function() {
  status_file <- file.path(gsea_dir, "GSEA_status.tsv")
  if (!file.exists(status_file)) return(invisible(NULL))
  st <- tryCatch(read.delim(status_file, stringsAsFactors=FALSE, check.names=FALSE,
                            quote="", comment.char=""), error=function(e) NULL)
  if (is.null(st) || !"contrast" %in% colnames(st)) return(invisible(NULL))
  for (type_suffix in c("GO_BP", "KEGG")) {
    sig_col <- paste0(type_suffix, "_sig")
    status_col <- paste0(type_suffix, "_status")
    if (!sig_col %in% colnames(st)) next
    for (ii in seq_len(nrow(st))) {
      if (status_col %in% colnames(st) && !st[[status_col]][ii] %in% c("SUCCESS", "SUCCESS_EMPTY")) next
      f <- file.path(gsea_dir, paste0(st$contrast[ii], "_", type_suffix, ".tsv"))
      if (!file.exists(f)) next
      d <- tryCatch(read.delim(f, stringsAsFactors=FALSE, check.names=FALSE,
                               quote="", comment.char=""), error=function(e) NULL)
      if (is.null(d) || !"p.adjust" %in% colnames(d)) next
      st[[sig_col]][ii] <- sum(as.numeric(d$p.adjust) < sig_threshold, na.rm=TRUE)
    }
  }
  write.table(st, status_file, sep="\t", row.names=FALSE, quote=FALSE)
  invisible(NULL)
}
refresh_gsea_status_counts()

process_type <- function(type_suffix) {
  cat(sprintf("%s: ", type_suffix))
  
  all_files  <- list.files(gsea_dir, pattern=paste0("_", type_suffix, "\\.tsv$"),
                           full.names=TRUE)
  all_files  <- all_files[!grepl("_treatment_summary\\.tsv$", all_files)]
  if (!length(all_files)) { cat("  No GSEA files found\n"); return(invisible(NULL)) }
  
  contrast_names <- gsub(paste0("_", type_suffix, "\\.tsv$"), "",
                         basename(all_files))
  cat(sprintf("%d GSEA file(s)", length(all_files)))
  
  role_by_label <- setNames(rep(NA_character_, length(contrast_names)), contrast_names)
  if (nrow(pipeline_contrasts) && all(c("label", "role") %in% colnames(pipeline_contrasts))) {
    manifest_roles <- setNames(as.character(pipeline_contrasts$role), as.character(pipeline_contrasts$label))
    role_by_label <- manifest_roles[contrast_names]
    names(role_by_label) <- contrast_names
  }
  # Classify by the recorded contrast role, not by the presence of `_vs_`:
  # interaction labels also contain `_vs_` (e.g. `aerobiosis_vs_anaerobiosis`).
  has_inter <- !is.na(role_by_label) & startsWith(role_by_label, "interaction_")
  # For outputs lacking a usable manifest, recognize interaction-label syntax.
  unknown_role <- is.na(role_by_label) | !nzchar(role_by_label)
  has_inter <- has_inter | (unknown_role &
    grepl("×", contrast_names, fixed=TRUE) &
    grepl("__", contrast_names, fixed=TRUE))
  has_vs <- grepl("_vs_", contrast_names, fixed=TRUE) & !has_inter
  
  vs_ref_files   <- all_files[has_vs]
  vs_ref_names   <- contrast_names[has_vs]
  inter_files    <- all_files[has_inter]
  inter_names    <- contrast_names[has_inter]
  pairwise_files <- all_files[has_vs]
  
  is_interaction_mode <- n_int_factors >= 1
  
  n_ifac_detected <- if (length(vs_ref_names) > 0) max(c(0L,
                                                         sapply(vs_ref_names, function(nm) {
                                                           lhs <- strsplit(nm, "_vs_", fixed=TRUE)[[1]][1]
                                                           max(0L, length(split_condition_components(lhs)) - 1L)
                                                         })
  )) else 0L
  if (is_interaction_mode) n_int_factors_eff <- max(n_int_factors, n_ifac_detected) else
    n_int_factors_eff <- 0L

  pairwise_simple_files <- character(0)
  pairwise_simple_names <- character(0)
  pairwise_marginal_files <- character(0)
  pairwise_marginal_names <- character(0)
  if (is_interaction_mode && length(vs_ref_names)) {
    roles <- setNames(rep(NA_character_, length(vs_ref_names)), vs_ref_names)
    if (nrow(pipeline_contrasts) && all(c("label", "role") %in% colnames(pipeline_contrasts))) {
      roles <- setNames(as.character(pipeline_contrasts$role), as.character(pipeline_contrasts$label))[vs_ref_names]
    }
    rhs_first <- vapply(vs_ref_names, function(nm) {
      z <- strsplit(nm, "_vs_", fixed=TRUE)[[1]]
      if (length(z) < 2L) return(NA_character_)
      split_condition_components(z[2])[1]
    }, character(1))
    is_pairwise_marginal <- roles == "pairwise_marginal" | startsWith(vs_ref_names, "marginal__")
    is_pairwise_marginal[is.na(is_pairwise_marginal)] <- FALSE
    pairwise_marginal_files <- vs_ref_files[is_pairwise_marginal]
    pairwise_marginal_names <- sub("^marginal__", "", vs_ref_names[is_pairwise_marginal])

    is_pairwise_simple <- roles == "pairwise_simple"
    is_pairwise_simple[is.na(is_pairwise_simple)] <- !is.na(reference_level_manifest) &
      rhs_first[is.na(is_pairwise_simple)] != reference_level_manifest
    is_pairwise_simple[is.na(is_pairwise_simple)] <- FALSE
    is_pairwise_simple <- is_pairwise_simple & !is_pairwise_marginal
    pairwise_simple_files <- vs_ref_files[is_pairwise_simple]
    pairwise_simple_names <- vs_ref_names[is_pairwise_simple]
    remove_pairwise <- is_pairwise_simple | is_pairwise_marginal
    vs_ref_files <- vs_ref_files[!remove_pairwise]
    vs_ref_names <- vs_ref_names[!remove_pairwise]
  }
  
  if (!is_interaction_mode) {
    ref_level <- detect_reference_classic(vs_ref_names)
    ref_sfx   <- paste0("_vs_", ref_level)
    is_vs_ref <- endsWith(vs_ref_names, ref_sfx)
    
    vs_ref_files_std <- vs_ref_files[is_vs_ref]
    vs_ref_nms_std   <- vs_ref_names[is_vs_ref]
    treatments       <- sub(paste0(ref_sfx, "$"), "", vs_ref_nms_std)
    
    pw_is <- has_vs & !endsWith(contrast_names, ref_sfx)
    pairwise_files <- all_files[pw_is]
    pairwise_names <- contrast_names[pw_is]
    
    
    data_by_treatment <- setNames(lapply(vs_ref_files_std, load_tsv), treatments)
    rows_loaded <- vapply(data_by_treatment, nrow, integer(1))

    all_ids <- Reduce(union, lapply(data_by_treatment, rownames), init=character(0))
    if (!length(all_ids)) {
      stop(sprintf("No %s terms were read from the versus-reference GSEA files.",
                   type_suffix), call.=FALSE)
    }
    all_desc <- character(length(all_ids)); names(all_desc) <- all_ids
    
    for (t in treatments) {
      df <- data_by_treatment[[t]]
      if (is.null(df)) next
      desc_col <- intersect(c("Description","description"), colnames(df))[1]
      if (!is.na(desc_col)) {
        present <- intersect(rownames(df), all_ids)
        all_desc[present] <- df[present, desc_col]
      }
    }
    
    nes_matrix  <- matrix(NA_real_,  nrow=length(all_ids), ncol=length(treatments),
                          dimnames=list(all_ids, treatments))
    padj_matrix <- matrix(NA_real_,  nrow=length(all_ids), ncol=length(treatments),
                          dimnames=list(all_ids, treatments))
    
    for (t in treatments) {
      df <- data_by_treatment[[t]]
      if (is.null(df)) next
      ids_here <- intersect(rownames(df), all_ids)
      nes_col  <- intersect(c("NES","nes"), colnames(df))[1]
      padj_col <- intersect(c("p.adjust","pvalue","padj"), colnames(df))[1]
      if (!is.na(nes_col))  nes_matrix[ids_here, t]  <- df[ids_here, nes_col]
      if (!is.na(padj_col)) padj_matrix[ids_here, t] <- df[ids_here, padj_col]
    }
    
    sig_matrix <- matrix(ifelse(!is.na(padj_matrix) & padj_matrix < sig_threshold, "*", "ns"),
                         nrow=nrow(padj_matrix), ncol=ncol(padj_matrix),
                         dimnames=dimnames(padj_matrix))
    
    pairwise_data <- list()
    for (i in seq_along(pairwise_files)) {
      df <- load_tsv(pairwise_files[i]); nm <- pairwise_names[i]
      if (!is.null(df)) pairwise_data[[nm]] <- df
    }
    
    use_cld <- !is_interaction_mode && length(treatments) >= 2
    cld_list <- if (use_cld) {
      setNames(lapply(all_ids, function(gid) {
        compute_cld(treatments, function(ta, tb) {
          lbl <- paste0(ta, "_vs_", tb)
          if (!lbl %in% names(pairwise_data)) lbl <- paste0(tb, "_vs_", ta)
          get_gsea_stats(pairwise_data[[lbl]], gid)
        })
      }), all_ids)
    } else NULL
    
    out <- data.frame(ID=all_ids, Description=all_desc, stringsAsFactors=FALSE,
                      check.names=FALSE)
    for (t in treatments) {
      out[[paste0(t, "_NES")]]        <- nes_matrix[all_ids, t]
      out[[paste0(t, "_sig_vs_ref")]] <- sig_matrix[all_ids, t]
      if (use_cld) out[[paste0(t, "_CLD")]] <- vapply(all_ids, function(g)
        cld_list[[g]][t], character(1))
    }
    
    if (length(treatments) > 1) {
      treatment_pairs <- combn(treatments, 2, simplify=FALSE)
      out[["treatment_differences"]] <- vapply(all_ids, function(gid) {
        stats <- lapply(treatment_pairs, function(p) {
          lbl <- paste0(p[1], "_vs_", p[2])
          if (!lbl %in% names(pairwise_data)) lbl <- paste0(p[2], "_vs_", p[1])
          get_gsea_stats(pairwise_data[[lbl]], gid)
        })
        raw_p <- vapply(stats, function(z)
          if (is.null(z) || is.na(z$raw)) 1 else as.numeric(z$raw), numeric(1))
        bh_q <- vapply(stats, function(z)
          if (is.null(z) || is.na(z$bh)) 1 else as.numeric(z$bh), numeric(1))
        holm_p <- p.adjust(raw_p, method="holm")
        if (any(bh_q < sig_threshold & holm_p < sig_threshold)) "*" else ""
      }, character(1))
    }
    
    out_file <- file.path(gsea_dir, paste0(type_suffix, "_treatment_summary.tsv"))
    write.table(out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
    cat(sprintf(" -> %d term(s) [%s]\n", nrow(out), basename(out_file)))
    return(invisible(out))
  }
  
  
  parse_condition_label <- function(nm) {
    parts <- strsplit(nm, "_vs_")[[1]]
    if (length(parts) < 2) return(NULL)
    lhs <- parts[1]; rhs <- parts[2]
    list(lhs=lhs, rhs=rhs,
         lhs_parts=split_condition_components(lhs),
         rhs_parts=split_condition_components(rhs))
  }
  
  cond_parsed <- lapply(vs_ref_names, parse_condition_label)
  valid       <- !vapply(cond_parsed, is.null, logical(1))
  cond_files  <- vs_ref_files[valid]
  cond_names  <- vs_ref_names[valid]
  cond_labels <- vapply(cond_parsed[valid], function(p) p$lhs, character(1))
  
  conditions  <- unique(cond_labels)
  
  data_by_cond <- setNames(lapply(conditions, function(cond) {
    idx <- which(cond_labels == cond)
    if (!length(idx)) return(NULL)
    load_tsv(cond_files[idx[1]])
  }), conditions)
  
  all_ids <- Reduce(
    union,
    lapply(data_by_cond, function(df) if (!is.null(df)) rownames(df) else character(0)),
    init = character(0)
  )
  if (!length(all_ids)) {
    stop(sprintf("No %s terms were read from the condition-specific GSEA files.",
                 type_suffix), call.=FALSE)
  }
  all_desc <- character(length(all_ids)); names(all_desc) <- all_ids
  
  for (cond in conditions) {
    df <- data_by_cond[[cond]]
    if (is.null(df)) next
    desc_col <- intersect(c("Description","description"), colnames(df))[1]
    if (!is.na(desc_col)) {
      present <- intersect(rownames(df), all_ids)
      all_desc[present] <- df[present, desc_col]
    }
  }
  
  nes_matrix  <- matrix(NA_real_, nrow=length(all_ids), ncol=length(conditions),
                        dimnames=list(all_ids, conditions))
  padj_matrix <- matrix(NA_real_, nrow=length(all_ids), ncol=length(conditions),
                        dimnames=list(all_ids, conditions))
  
  for (cond in conditions) {
    df <- data_by_cond[[cond]]
    if (is.null(df)) next
    ids_here <- intersect(rownames(df), all_ids)
    nes_col  <- intersect(c("NES","nes"), colnames(df))[1]
    padj_col <- intersect(c("p.adjust","pvalue","padj"), colnames(df))[1]
    if (!is.na(nes_col))  nes_matrix[ids_here, cond]  <- df[ids_here, nes_col]
    if (!is.na(padj_col)) padj_matrix[ids_here, cond] <- df[ids_here, padj_col]
  }
  
  sig_matrix <- matrix(
    ifelse(!is.na(padj_matrix) & padj_matrix < sig_threshold, "*", "ns"),
    nrow=nrow(padj_matrix), ncol=ncol(padj_matrix),
    dimnames=dimnames(padj_matrix))
  
  pairwise_simple_data <- setNames(lapply(pairwise_simple_files, load_tsv),
                                   pairwise_simple_names)
  pairwise_marginal_data <- setNames(lapply(pairwise_marginal_files, load_tsv),
                                     pairwise_marginal_names)
  contrast_levels <- unique(vapply(conditions, function(z)
    split_condition_components(z)[1], character(1)))
  use_cld_interaction <- length(contrast_levels) > 2L &&
    (length(pairwise_simple_data) > 0L || length(pairwise_marginal_data) > 0L)
  cld_stratified <- NULL
  cld_marginal <- NULL
  strata <- unique(vapply(conditions, function(z) {
    pp <- split_condition_components(z)
    paste(pp[-1], collapse="×")
  }, character(1)))

  if (use_cld_interaction && length(pairwise_simple_data) > 0L) {
    cld_stratified <- matrix("", nrow=length(all_ids), ncol=length(conditions),
                             dimnames=list(all_ids, conditions))
    for (stratum_index in seq_along(strata)) {
      stratum <- strata[stratum_index]
      conds_here <- conditions[vapply(conditions, function(z) {
        pp <- split_condition_components(z)
        identical(paste(pp[-1], collapse="×"), stratum)
      }, logical(1))]
      lvls_here <- vapply(conds_here, function(z) split_condition_components(z)[1], character(1))
      if (length(lvls_here) < 2L) next
      prime_suffix <- if (stratum_index == 1L) "" else
        paste(rep("′", stratum_index - 1L), collapse="")
      for (gid in all_ids) {
        letters_here <- compute_cld(lvls_here, function(la, lb) {
          lhs <- paste(c(la, if(nzchar(stratum)) strsplit(stratum, "×", fixed=TRUE)[[1]] else character(0)), collapse="×")
          rhs <- paste(c(lb, if(nzchar(stratum)) strsplit(stratum, "×", fixed=TRUE)[[1]] else character(0)), collapse="×")
          nm <- paste0(lhs, "_vs_", rhs)
          if (!nm %in% names(pairwise_simple_data)) nm <- paste0(rhs, "_vs_", lhs)
          get_gsea_stats(pairwise_simple_data[[nm]], gid)
        })
        for (k in seq_along(conds_here))
          cld_stratified[gid, conds_here[k]] <- paste0(letters_here[lvls_here[k]], prime_suffix)
      }
    }
  }

  if (use_cld_interaction && length(pairwise_marginal_data) > 0L) {
    cld_marginal_levels <- matrix("", nrow=length(all_ids), ncol=length(contrast_levels),
                                  dimnames=list(all_ids, contrast_levels))
    for (gid in all_ids) {
      letters_here <- compute_cld(contrast_levels, function(la, lb) {
        nm <- paste0(la, "_vs_", lb)
        if (!nm %in% names(pairwise_marginal_data)) nm <- paste0(lb, "_vs_", la)
        get_gsea_stats(pairwise_marginal_data[[nm]], gid)
      })
      cld_marginal_levels[gid, contrast_levels] <- letters_here[contrast_levels]
    }
    cld_marginal <- matrix("", nrow=length(all_ids), ncol=length(conditions),
                           dimnames=list(all_ids, conditions))
    for (cond in conditions) {
      lvl <- split_condition_components(cond)[1]
      cld_marginal[, cond] <- cld_marginal_levels[, lvl]
    }
  }

  split_interaction_name <- function(nm) {
    core <- strsplit(nm, "__", fixed=TRUE)[[1]][1]
    if (grepl("×", core, fixed=TRUE)) strsplit(core, "×", fixed=TRUE)[[1]]
    else core
  }
  int_parts <- lapply(inter_names, split_interaction_name)
  int_types <- vapply(int_parts, function(z) max(0L, length(z)-1L), integer(1))

  int_2way_files <- inter_files[int_types == 1L]
  int_2way_names <- inter_names[int_types == 1L]
  int_3way_files <- inter_files[int_types >= 2L]
  int_3way_names <- inter_names[int_types >= 2L]

  load_interaction_stats <- function(files, names) {
    if (!length(files)) return(list())
    setNames(lapply(files, function(f) {
      df <- load_tsv(f); if (is.null(df)) return(NULL)
      padj_col <- intersect(c("p.adjust","pvalue","padj"), colnames(df))[1]
      nes_col  <- intersect(c("NES","nes"), colnames(df))[1]
      if (is.na(padj_col)) return(NULL)
      list(
        padj=setNames(as.numeric(df[[padj_col]]), rownames(df)),
        nes=if (!is.na(nes_col)) setNames(as.numeric(df[[nes_col]]), rownames(df)) else
          setNames(rep(NA_real_, nrow(df)), rownames(df))
      )
    }), names)
  }

  int_2way_stats <- load_interaction_stats(int_2way_files, int_2way_names)
  int_3way_stats <- load_interaction_stats(int_3way_files, int_3way_names)
  int_2way_data <- lapply(int_2way_stats, function(z) if(is.null(z)) NULL else z$padj)
  int_3way_data <- lapply(int_3way_stats, function(z) if(is.null(z)) NULL else z$padj)

  get_one_int_sig <- function(pvec, gid) {
    if (is.null(pvec)) return(FALSE)
    pv <- pvec[gid]
    !is.null(pv) && !is.na(pv) && pv < sig_threshold
  }


  cld_by_condition <- NULL
  if (use_cld_interaction) {
    cld_by_condition <- matrix("", nrow=length(all_ids), ncol=length(conditions),
                               dimnames=list(all_ids, conditions))
    for (gid in all_ids) {
      any_specific_interaction <- any(vapply(int_2way_data, function(pvec)
        get_one_int_sig(pvec, gid), logical(1)))
      source_mat <- if (any_specific_interaction && !is.null(cld_stratified))
        cld_stratified else cld_marginal
      if (is.null(source_mat)) source_mat <- cld_stratified
      if (!is.null(source_mat)) cld_by_condition[gid, ] <- source_mat[gid, conditions]
    }
  }

  out <- data.frame(ID=all_ids, Description=all_desc, stringsAsFactors=FALSE,
                    check.names=FALSE)
  for (cond in conditions) {
    out[[paste0(cond, "_NES")]]        <- nes_matrix[all_ids, cond]
    out[[paste0(cond, "_sig_vs_ref")]] <- sig_matrix[all_ids, cond]
    if (use_cld_interaction)
      out[[paste0(cond, "_CLD")]] <- cld_by_condition[all_ids, cond]
  }

  if (length(int_2way_names)) {
    int_meta <- lapply(int_2way_names, function(nm) {
      parts <- split_interaction_name(nm)
      row_m <- if (nrow(pipeline_contrasts) && "label" %in% names(pipeline_contrasts))
        pipeline_contrasts[pipeline_contrasts$label == nm, , drop=FALSE] else data.frame()
      val <- function(col, fallback=NA_character_) {
        if (nrow(row_m) && col %in% names(row_m)) {
          z <- as.character(row_m[[col]][1])
          if (!is.na(z) && nzchar(z)) return(z)
        }
        fallback
      }
      list(
        name=nm,
        cf=val("contrast_level", if(length(parts)) parts[1] else NA_character_),
        factor=val("interaction_factor", if(length(parts)>=2) parts[2] else NA_character_),
        tested=val("tested_level"),
        reference=val("tested_reference"),
        conditioned_factor=val("conditioned_factor"),
        conditioned_level=val("conditioned_level")
      )
    })
    names(int_meta) <- int_2way_names

    group_key <- vapply(int_meta, function(m)
      paste(m$cf, m$factor, m$tested, m$reference, sep="|||"), character(1))
    int_groups <- split(int_2way_names, group_key)

    for (grp in int_groups) {
      m <- int_meta[[grp[1]]]
      if (is.na(m$cf) || is.na(m$factor)) next
      col_name <- paste0("interaction_x_", m$factor, "__", m$cf,
                         if(!is.na(m$tested) && nzchar(m$tested)) paste0("__", m$tested) else "",
                         if(!is.na(m$reference) && nzchar(m$reference)) paste0("__", m$reference) else "")
      out[[col_name]] <- vapply(all_ids, function(gid) {
        sigs <- vapply(grp, function(nm) get_one_int_sig(int_2way_data[[nm]], gid), logical(1))
        if (!all(sigs)) return("ns")
        if (length(grp) > 1L) {
          nes_vals <- vapply(grp, function(nm) {
            z <- int_2way_stats[[nm]]
            if (is.null(z)) return(NA_real_)
            v <- z$nes[gid]
            if (length(v) && is.finite(v)) as.numeric(v) else NA_real_
          }, numeric(1))
          if (!all(is.finite(nes_vals)) || length(unique(sign(nes_vals))) != 1L) return("ns")
        }
        "*"
      }, character(1))
    }
  }

  if (length(int_3way_names) > 0L) {
    for (nm in int_3way_names) {
      parts <- split_interaction_name(nm)
      cf_lvl <- parts[1]
      facs <- parts[-1]
      col_name <- paste0("interaction_3way__", cf_lvl, "__", paste(facs, collapse="__"))
      pvec <- int_3way_data[[nm]]
      out[[col_name]] <- vapply(all_ids, function(gid) {
        if (get_one_int_sig(pvec, gid)) "*" else "ns"
      }, character(1))
    }
  }

  out_file <- file.path(gsea_dir, paste0(type_suffix, "_treatment_summary.tsv"))
  write.table(out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf(" -> %d term(s) [%s]\n", nrow(out), basename(out_file)))
  invisible(out)
}

detect_reference_classic <- function(contrast_names) {
  rhs <- vapply(contrast_names, function(cn) {
    parts <- strsplit(cn, "_vs_")[[1]]
    if (length(parts) < 2) NA_character_
    else paste(parts[2:length(parts)], collapse="_vs_")
  }, character(1))
  rhs <- rhs[!is.na(rhs)]
  names(which.max(table(rhs)))
}

process_type("GO_BP")
process_type("KEGG")


invisible(NULL)
