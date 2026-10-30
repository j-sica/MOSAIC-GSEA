#!/usr/bin/env Rscript
# Main workflow




options(stringsAsFactors = FALSE, warn = 1)

home_dir <- Sys.getenv("MOSAIC_GSEA_HOME", unset = "")
version_file <- if (nzchar(home_dir)) file.path(home_dir, "VERSION") else "VERSION"
mosaic_version <- if (file.exists(version_file)) trimws(readLines(version_file, n = 1L, warn = FALSE)) else "1.0.0"
if (!nzchar(mosaic_version)) mosaic_version <- "1.0.0"

usage <- function(status = 0L) {
  cat(
paste0("MOSAIC-GSEA ", mosaic_version, "\n\n"),
"Usage:\n",
"  mosaic-gsea -i COUNTS -m METADATA -f FORMULA -c FACTOR -r REFERENCE [options]\n\n",
"Required:\n",
"  -i, --input FILE              Count matrix; first column contains gene IDs\n",
"  -m, --metadata FILE           Metadata; first column contains sample IDs\n",
"  -f, --formula STRING          R model formula, e.g. '~ treatment' or '~ strain*treatment*oxygen'\n",
"  -c, --contrast-factor NAME    Factor whose non-reference levels are contrasted\n",
"  -r, --reference LEVEL         Reference level of the contrast factor\n\n",
"Output:\n",
"  -o, --output DIR              Run directory [./mosaic_gsea_output]\n",
"      --clean BOOL              Remove an existing output directory [false]\n",
"      --redo                     Reuse existing results and rerun from the earliest affected stage\n",
"      --keep-configured BOOL    Retain configured internal scripts [false]\n",
"      --cache-dir DIR           Parent directory for organism annotation data [auto]\n\n",
"Model options:\n",
"      --interaction-factors CSV auto, none, or comma-separated names [auto]\n",
"      --cores N                 BiocParallel workers for model/variance fitting [4]\n",
"      --filter-column NAME      Optional metadata filter column\n",
"      --filter-value VALUE      Optional metadata filter value\n",
"      --variable-types FILE     Optional column types file: factor or number [auto]\n\n",
"Organism options:\n",
"      --organism-name STRING    Scientific name [Saccharomyces cerevisiae]\n",
"      --orgdb PACKAGE           Bioconductor OrgDb package [auto]\n",
"      --key-type TYPE           Input gene ID type [auto: detect using OrgDb]\n",
"      --symbol-column COLUMN    Display gene-name column [GENENAME for yeast, SYMBOL otherwise]\n",
"      --kegg-org CODE           KEGG organism code [auto]\n",
"      --yeast-id-style STYLE    systematic or common [systematic]\n\n",
"GSEA options:\n",
"      --rank-by METRIC          t or logFC [t; dream uses z.std for t-ranking]\n",
"      --min-gs N                Minimum gene-set size [10]\n",
"      --max-gs N                Maximum gene-set size [500]\n",
"      --padj-cut X              GSEA storage cutoff [1.0]\n",
"      --padj-method NAME        Multiple-testing method [BH]\n",
"      --run-go BOOL             Run GO BP [true]\n",
"      --run-kegg BOOL           Run KEGG [true]\n\n",
"Report/filter options:\n",
"      --sig-threshold X         Adjusted-p significance threshold [0.05]\n",
"      --filter-mode MODE        interpretation or all [interpretation]\n",
"      --similarity-reduction BOOL  Reduce redundant GO terms [true]\n",
"      --similarity-threshold X  GO redundancy threshold [0.375]\n",
"      --similarity-metric NAME  jaccard_overlap, jaccard, overlap, dice [jaccard_overlap]\n",
"      --winner-score NAME       score_weighted, max_abs_nes, mean_abs_nes, n_significant [score_weighted]\n",
"      --max-nes-scale X         Fixed NES color limit; omit for automatic\n\n",
"Other:\n",
"  -h, --help                    Show this help\n",
"      --version                 Show MOSAIC-GSEA version\n\n",
"Diagonal/cross-condition contrasts are intentionally not generated.\n",
sep = "")
  quit(save = "no", status = status)
}

parse_bool <- function(x, name) {
  z <- tolower(trimws(x))
  if (z %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (z %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("%s must be true or false; received '%s'", name, x), call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (any(args == "--version")) {
  cat(sprintf("MOSAIC-GSEA %s\n", mosaic_version))
  quit(save = "no", status = 0L)
}
if (!length(args) || any(args %in% c("-h", "--help"))) usage(0L)
cat(sprintf("MOSAIC-GSEA %s\nInitializing...\n", mosaic_version))
flush.console()

short_map <- c(i = "input", m = "metadata", f = "formula", c = "contrast-factor",
               r = "reference", o = "output")
allowed_options <- c(
  "input", "metadata", "formula", "contrast-factor", "reference", "output",
  "clean", "redo", "keep-configured", "cache-dir", "interaction-factors", "interaction-refs",
  "filter-column", "filter-value", "variable-types", "organism-name", "orgdb",
  "key-type", "symbol-column", "kegg-org", "yeast-id-style", "rank-by",
  "min-gs", "max-gs", "padj-cut", "padj-method", "run-go", "run-kegg",
  "sig-threshold", "filter-mode", "similarity-reduction", "similarity-threshold",
  "similarity-metric", "winner-score", "max-nes-scale", "cores"
)
opts <- list()
i <- 1L
while (i <= length(args)) {
  token <- args[[i]]
  if (startsWith(token, "--")) {
    keyval <- substring(token, 3L)
    if (grepl("=", keyval, fixed = TRUE)) {
      pieces <- strsplit(keyval, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    } else {
      key <- keyval
      if (identical(key, "redo") &&
          (i == length(args) || startsWith(args[[i + 1L]], "-"))) {
        value <- "true"
      } else {
        if (i == length(args) || startsWith(args[[i + 1L]], "-"))
          stop(sprintf("Missing value for --%s", key), call. = FALSE)
        i <- i + 1L
        value <- args[[i]]
      }
    }
  } else if (grepl("^-[imfcro]$", token)) {
    short <- substring(token, 2L)
    key <- unname(short_map[[short]])
    if (i == length(args)) stop(sprintf("Missing value for %s", token), call. = FALSE)
    i <- i + 1L
    value <- args[[i]]
  } else {
    stop(sprintf("Unknown argument: %s", token), call. = FALSE)
  }
  if (!key %in% allowed_options)
    stop(sprintf("Unknown option: --%s", key), call. = FALSE)
  if (!is.null(opts[[key]]))
    stop(sprintf("Option supplied more than once: --%s", key), call. = FALSE)
  opts[[key]] <- value
  i <- i + 1L
}

if (!is.null(opts[["interaction-refs"]])) {
  stop("Interaction-factor baselines are selected automatically from the first observed level in the metadata; omit --interaction-refs.", call. = FALSE)
}

required <- c("input", "metadata", "formula", "contrast-factor", "reference")
missing <- required[!vapply(required, function(x) !is.null(opts[[x]]) && nzchar(opts[[x]]), logical(1))]
if (length(missing)) stop(sprintf("Missing required option(s): %s", paste(missing, collapse = ", ")), call. = FALSE)

validate_supported_formula <- function(fstr) {
  # Accept only formula syntax supported by the automatic contrast parser.
  fixed <- gsub("\\([^)]*\\|[^)]*\\)", "", fstr, perl=TRUE)
  bad <- character(0)
  if (grepl("\\^", fixed)) bad <- c(bad, "^")
  if (grepl("/", fixed, fixed=TRUE)) bad <- c(bad, "/")
  if (grepl("%in%", fixed, fixed=TRUE)) bad <- c(bad, "%in%")
  if (grepl("(^|[^A-Za-z0-9_.])(poly|I)\\s*\\(", fixed, perl=TRUE)) bad <- c(bad, "function call")
  if (grepl("[()]", fixed)) bad <- c(bad, "parenthesized/function fixed term")
  if (length(bad)) {
    stop(sprintf(
      "Unsupported model-formula syntax (%s). Use explicit terms joined by +, * or :, with optional simple random-effect terms such as (1|group).",
      paste(unique(bad), collapse=", ")
    ), call.=FALSE)
  }
}
validate_supported_formula(opts[["formula"]])

getopt <- function(name, default = NULL) if (!is.null(opts[[name]])) opts[[name]] else default
abs_path <- function(path, must_exist = FALSE) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(getwd(), path)
  path <- normalizePath(path, winslash = "/", mustWork = must_exist)
  path
}
r_string <- function(x) encodeString(as.character(x), quote = '"')
r_nullable_string <- function(x) if (is.null(x) || !nzchar(x)) "NULL" else r_string(x)
r_vector <- function(x) {
  if (is.null(x)) return("NULL")
  if (!length(x)) return("character(0)")
  paste0("c(", paste(vapply(x, r_string, character(1)), collapse = ", "), ")")
}

script_home <- Sys.getenv("MOSAIC_GSEA_HOME", unset = "")
if (!nzchar(script_home)) {
  file_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- file_args[grep("^--file=", file_args)][1]
  file_arg <- sub("^--file=", "", file_arg)
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  script_home <- dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
}
script_home <- normalizePath(script_home, winslash = "/", mustWork = TRUE)
source(file.path(script_home, "organism_support.R"), local = TRUE)

organism_name <- getopt("organism-name", "Saccharomyces cerevisiae")
orgdb <- resolve_orgdb(organism_name, getopt("orgdb", "auto"))
ensure_orgdb(orgdb, organism_name)
run_kegg_requested <- parse_bool(getopt("run-kegg", "true"), "--run-kegg")
kegg_org <- resolve_kegg_org(organism_name, getopt("kegg-org", "auto"), run_kegg_requested)

cfg <- list(
  input = abs_path(getopt("input"), TRUE),
  metadata = abs_path(getopt("metadata"), TRUE),
  formula = getopt("formula"),
  contrast_factor = getopt("contrast-factor"),
  reference = getopt("reference"),
  output = abs_path(getopt("output", "mosaic_gsea_output"), FALSE),
  cache_dir_option = getopt("cache-dir", ""),
  clean = parse_bool(getopt("clean", "false"), "--clean"),
  keep_configured = parse_bool(getopt("keep-configured", "false"), "--keep-configured"),
  filter_column = getopt("filter-column", ""),
  filter_value = getopt("filter-value", ""),
  variable_types_file = if (nzchar(getopt("variable-types", ""))) abs_path(getopt("variable-types"), TRUE) else "",
  interaction_factors = getopt("interaction-factors", "auto"),
  organism_name = organism_name,
  orgdb = orgdb,
  key_type = getopt("key-type", "auto"),
  symbol_column = getopt("symbol-column", if (orgdb == "org.Sc.sgd.db") "GENENAME" else "SYMBOL"),
  kegg_org = kegg_org,
  yeast_id_style = getopt("yeast-id-style", "systematic"),
  rank_by = getopt("rank-by", "t"),
  min_gs = as.integer(getopt("min-gs", "10")),
  max_gs = as.integer(getopt("max-gs", "500")),
  padj_cut = as.numeric(getopt("padj-cut", "1.0")),
  padj_method = getopt("padj-method", "BH"),
  run_go = parse_bool(getopt("run-go", "true"), "--run-go"),
  run_kegg = run_kegg_requested,
  sig_threshold = as.numeric(getopt("sig-threshold", "0.05")),
  filter_mode = getopt("filter-mode", "interpretation"),
  similarity_reduction = parse_bool(getopt("similarity-reduction", "true"), "--similarity-reduction"),
  similarity_threshold = as.numeric(getopt("similarity-threshold", "0.375")),
  similarity_metric = getopt("similarity-metric", "jaccard_overlap"),
  winner_score = getopt("winner-score", "score_weighted"),
  max_nes_scale = getopt("max-nes-scale", ""),
  cores = suppressWarnings(as.integer(getopt("cores", "4")))
)
redo_requested <- parse_bool(getopt("redo", "false"), "--redo")

`%notin%` <- function(x, table) !(x %in% table)

if (cfg$yeast_id_style %notin% c("systematic", "common")) stop("--yeast-id-style must be systematic or common", call. = FALSE)

# Resolve identifier namespace from actual count-matrix IDs, not from a fixed
# OrgDb central key. An explicit --key-type always takes precedence.
orgdb_obj <- getExportedValue(cfg$orgdb, cfg$orgdb)
available_keytypes <- AnnotationDbi::keytypes(orgdb_obj)
if (!identical(tolower(cfg$key_type), "auto") && !cfg$key_type %in% available_keytypes) {
  stop(sprintf("--key-type '%s' is not supported by %s. Available: %s",
               cfg$key_type, cfg$orgdb, paste(available_keytypes, collapse=", ")), call.=FALSE)
}
read_input_id_sample <- function(path, limit=Inf, scan_lines=-1L) {
  con <- if (grepl("\\.gz$", path, ignore.case=TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add=TRUE)
  lines <- readLines(con, n=scan_lines, warn=FALSE)
  if (length(lines) < 2L) stop("Count matrix has no gene rows", call.=FALSE)
  ids <- vapply(strsplit(lines[-1L], "\t", fixed=TRUE), function(x) x[[1L]], character(1))
  ids <- unique(trimws(ids[nzchar(trimws(ids))]))
  if (is.finite(limit) && length(ids) > limit) ids <- ids[unique(round(seq(1, length(ids), length.out=limit)))]
  ids
}
input_id_sample <- read_input_id_sample(cfg$input)
keytypes_to_test <- intersect(c("ORF", "COMMON", "ENSEMBL", "ENTREZID", "SYMBOL", "SGD", "TAIR", "REFSEQ", "ALIAS", "ENSEMBLTRANS", "ENSEMBLPROT", "FLYBASE", "WORMBASE", "MGI", "UNIPROT"), available_keytypes)
if (!length(keytypes_to_test)) keytypes_to_test <- available_keytypes
id_coverage <- function(kt) {
  ids <- input_id_sample
  if (kt == "ENSEMBL") ids <- sub("\\.[0-9]+$", "", ids)
  known <- tryCatch(AnnotationDbi::keys(orgdb_obj, keytype=kt), error=function(e) character())
  if (!length(ids)) return(0)
  sum(ids %in% known) / length(ids)
}
if (identical(tolower(cfg$key_type), "auto")) {
  cover <- vapply(keytypes_to_test, id_coverage, numeric(1))
  cat("Input gene ID coverage by namespace: ", paste(sprintf("%s %.1f%%", names(sort(cover, decreasing=TRUE)), 100*sort(cover, decreasing=TRUE)), collapse=", "), "\n", sep="")
  best <- max(cover)
  if (!is.finite(best) || best < 0.50) {
    stop(sprintf("Cannot safely detect gene ID type for %s: best OrgDb coverage is %.1f%%. Specify --key-type explicitly.",
                 cfg$orgdb, 100*best), call.=FALSE)
  }
  near <- names(cover)[cover >= best - 0.02]
  # SGD identifiers may occur in both ORF and ENSEMBL: prefer the native ORF key.
  preferred <- if (identical(cfg$orgdb, "org.Sc.sgd.db")) {
    if (identical(cfg$yeast_id_style, "common")) "COMMON" else "ORF"
  } else ""
  if (nzchar(preferred) && preferred %in% near) {
    cfg$key_type <- preferred
  } else if (length(near) == 1L) {
    cfg$key_type <- near[[1L]]
  } else {
    stop(sprintf("Ambiguous input ID namespace (%s). Specify --key-type explicitly.",
                 paste(sprintf("%s %.1f%%", near, 100*cover[near]), collapse=", ")), call.=FALSE)
  }
  cat(sprintf("Gene ID type auto-detected: %s (%.1f%% matched in %s; %d input IDs inspected)\n",
              cfg$key_type, 100*cover[[cfg$key_type]], cfg$orgdb, length(input_id_sample)))
  if (best < 0.80) warning(sprintf("Only %.1f%% of input IDs match the best namespace. Check whether the selected OrgDb and gene IDs are compatible.", 100*best), call.=FALSE)
} else {
  coverage <- id_coverage(cfg$key_type)
  cat(sprintf("Gene ID type explicitly selected: %s (%.1f%% sample coverage in %s)\n",
              cfg$key_type, 100*coverage, cfg$orgdb))
  if (coverage < 0.80) warning(sprintf("Only %.1f%% of input genes match --key-type %s in %s; check the identifier namespace.",
                               100*coverage, cfg$key_type, cfg$orgdb), call.=FALSE)
}
if (identical(cfg$orgdb, "org.Sc.sgd.db") && is.null(opts[["symbol-column"]])) {
  cfg$symbol_column <- if (identical(cfg$key_type, "COMMON")) "COMMON" else "GENENAME"
}
if (!cfg$symbol_column %in% AnnotationDbi::columns(orgdb_obj)) {
  warning(sprintf("Display-name column %s is unavailable in %s; gene identifiers will be displayed instead.", cfg$symbol_column, cfg$orgdb), call.=FALSE)
  cfg$symbol_column <- cfg$key_type
}
if (!cfg$rank_by %in% c("t", "logFC")) stop("--rank-by must be t or logFC", call. = FALSE)
if (!cfg$filter_mode %in% c("interpretation", "all")) stop("--filter-mode must be interpretation or all", call. = FALSE)
if (!cfg$similarity_metric %in% c("jaccard_overlap", "combined", "jaccard", "overlap", "dice")) stop("Unsupported --similarity-metric", call. = FALSE)
if (!cfg$winner_score %in% c("score_weighted", "max_abs_nes", "mean_abs_nes", "n_significant")) stop("Unsupported --winner-score", call. = FALSE)
if (is.na(cfg$min_gs) || is.na(cfg$max_gs) || cfg$min_gs < 1L || cfg$max_gs < cfg$min_gs) stop("Invalid gene-set size limits", call. = FALSE)
if (!is.finite(cfg$padj_cut) || cfg$padj_cut <= 0 || cfg$padj_cut > 1) stop("--padj-cut must be greater than 0 and at most 1", call. = FALSE)
if (!cfg$padj_method %in% p.adjust.methods) stop(sprintf("--padj-method must be one of: %s", paste(p.adjust.methods, collapse=", ")), call. = FALSE)
if (!is.finite(cfg$sig_threshold) || cfg$sig_threshold <= 0 || cfg$sig_threshold >= 1) stop("--sig-threshold must be between 0 and 1", call. = FALSE)
if (cfg$padj_cut < cfg$sig_threshold) warning(sprintf("--padj-cut (%.4g) is below --sig-threshold (%.4g): GSEA will discard terms before MOSAIC-GSEA can evaluate them at the requested significance threshold.", cfg$padj_cut, cfg$sig_threshold), call. = FALSE)
if (!is.finite(cfg$similarity_threshold) || cfg$similarity_threshold < 0 || cfg$similarity_threshold > 1) stop("--similarity-threshold must be between 0 and 1", call. = FALSE)
if (nzchar(cfg$max_nes_scale)) {
  max_nes_num <- suppressWarnings(as.numeric(cfg$max_nes_scale))
  if (!is.finite(max_nes_num) || max_nes_num <= 0) stop("--max-nes-scale must be a positive number", call. = FALSE)
}
if (xor(nzchar(cfg$filter_column), nzchar(cfg$filter_value))) stop("--filter-column and --filter-value must be supplied together", call. = FALSE)
if (!cfg$run_go && !cfg$run_kegg) stop("At least one of --run-go or --run-kegg must be true", call. = FALSE)
cores_raw <- getopt("cores", "4")
if (!grepl("^[1-9][0-9]*$", cores_raw) || is.na(cfg$cores) || cfg$cores < 1L) stop("--cores must be an integer greater than or equal to 1", call. = FALSE)


if (cfg$clean && redo_requested) {
  stop("--clean and --redo cannot be used together.", call. = FALSE)
}
if (dir.exists(cfg$output)) {
  nonempty <- length(list.files(cfg$output, all.files = TRUE, no.. = TRUE)) > 0L
  if (nonempty && !cfg$clean && !redo_requested) {
    stop(sprintf(
      "Output directory already exists and is not empty: %s\nUse --redo to reuse valid upstream results, --clean true to restart, or choose another -o directory.",
      cfg$output
    ), call. = FALSE)
  }
  if (cfg$clean) unlink(cfg$output, recursive = TRUE, force = TRUE)
}

organism_slug <- tolower(cfg$organism_name)
organism_slug <- gsub("[^a-z0-9]+", "_", organism_slug)
organism_slug <- gsub("^_+|_+$", "", organism_slug)
annotation_dir_name <- paste0(organism_slug, "_GO_and_KEGG_data")

base_dir <- script_home

cache_parent <- if (nzchar(cfg$cache_dir_option)) {
  abs_path(cfg$cache_dir_option, FALSE)
} else {
  candidate_parents <- unique(c(base_dir, getwd(), dirname(cfg$output)))
  existing <- candidate_parents[dir.exists(file.path(candidate_parents, annotation_dir_name))]
  if (length(existing)) existing[[1]] else base_dir
}
cache_root <- file.path(cache_parent, annotation_dir_name)

dirs <- list(
  root = cfg$output,
  de = file.path(cfg$output, "01_differential_expression"),
  gsea = file.path(cfg$output, "02_gsea"),
  final = file.path(cfg$output, "03_final_reports"),
  cache = cache_root,
  go_cache = file.path(cache_root, "go", cfg$orgdb, cfg$key_type),
  kegg_cache = file.path(cache_root, "kegg", cfg$kegg_org, cfg$key_type),
  logs = file.path(cfg$output, "logs"),
  config = file.path(cfg$output, "config"),
  internal = file.path(cfg$output, ".pipeline_internal")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
sgd_cache_file <- file.path(cache_root, paste0("sgd_pathway_cache_", cfg$kegg_org, ".rds"))
if (identical(cfg$orgdb, "org.Sc.sgd.db") && identical(cfg$kegg_org, "sce") && !file.exists(sgd_cache_file)) {
  bundled_sgd_cache <- file.path(base_dir, annotation_dir_name, "sgd_pathway_cache_sce.rds")
  if (file.exists(bundled_sgd_cache)) {
    dir.create(dirname(sgd_cache_file), recursive=TRUE, showWarnings=FALSE)
    copied <- file.copy(bundled_sgd_cache, sgd_cache_file, overwrite=FALSE)
    if (!isTRUE(copied) && !file.exists(sgd_cache_file)) {
      warning("Bundled SGD biochemical pathway cache could not be copied to the active annotation cache.", call.=FALSE)
    }
  }
}


# Rebuild malformed SGD pathway caches before report generation.
sgd_cache_is_valid <- function(path) {
  if (!file.exists(path)) return(FALSE)
  x <- tryCatch(readRDS(path), error=function(e) NULL)
  if (!is.list(x) || is.null(x$pathways) || !is.list(x$pathways)) return(FALSE)
  nm <- unique(unlist(lapply(x$pathways, function(g) {
    if (is.null(g) || !length(g)) return(character(0))
    vapply(g, function(e) {
      if (is.list(e) && !is.null(e$name) && length(e$name)) as.character(e$name[[1]]) else ""
    }, character(1))
  }), use.names=FALSE))
  nm <- trimws(nm[nzchar(trimws(nm))])
  if (length(nm) < 10L) return(FALSE)
  ref_like <- grepl("(?:PMID|SGD_REF):", nm, ignore.case=TRUE, perl=TRUE)
  mean(ref_like) < 0.10
}

if (identical(cfg$orgdb, "org.Sc.sgd.db") && identical(cfg$kegg_org, "sce")) {
  needs_sgd_rebuild <- !sgd_cache_is_valid(sgd_cache_file)
  if (isTRUE(needs_sgd_rebuild)) {
    if (file.exists(sgd_cache_file)) {
      cat("Detected an invalid/outdated SGD pathway cache (for example PMID strings used as pathway names). Rebuilding it...\n")
      unlink(sgd_cache_file, force=TRUE)
    } else {
      cat("No valid SGD biochemical pathway cache found. Building it from the current SGD bulk table...\n")
    }
    builder <- file.path(base_dir, "build_sgd_pathway_cache.R")
    rscript_bin <- file.path(R.home("bin"), "Rscript")
    rebuilt <- FALSE
    if (file.exists(builder) && file.exists(rscript_bin)) {
      status <- tryCatch(
        system2(
          rscript_bin,
          args=c(shQuote(builder), shQuote(sgd_cache_file)),
          stdout="", stderr=""
        ),
        error=function(e) 1L
      )
      rebuilt <- identical(as.integer(status), 0L) && sgd_cache_is_valid(sgd_cache_file)
    }
    if (!isTRUE(rebuilt)) {
      if (file.exists(sgd_cache_file) && !sgd_cache_is_valid(sgd_cache_file)) unlink(sgd_cache_file, force=TRUE)
      warning(
        paste0(
          "Could not build a valid SGD biochemical pathway cache. ",
          "The analysis will continue, but the SGD Pathways panel will be unavailable. ",
          "When internet access is available, run: Rscript build_sgd_pathway_cache.R"
        ),
        call.=FALSE
      )
    }
  }
}

template_dir <- file.path(base_dir, "templates")
required_templates <- c("step1_de.R", "step2_gsea.R", "step3_summary.R", "step4_reports.R")
if (!all(file.exists(file.path(template_dir, required_templates)))) stop("Pipeline templates are missing next to the CLI script.", call. = FALSE)

code_files <- c(
  launcher = file.path(base_dir, "mosaic-gsea"),
  driver = file.path(base_dir, "mosaic_gsea.R"),
  organism_support = file.path(base_dir, "organism_support.R"),
  step1 = file.path(template_dir, "step1_de.R"),
  step2 = file.path(template_dir, "step2_gsea.R"),
  step3 = file.path(template_dir, "step3_summary.R"),
  step4 = file.path(template_dir, "step4_reports.R")
)
code_stage <- c(launcher=1L, driver=1L, organism_support=1L, step1=1L, step2=2L, step3=3L, step4=4L)
if (!all(file.exists(code_files))) stop("One or more pipeline source files are missing.", call. = FALSE)
current_code_signatures <- data.frame(
  component = names(code_files),
  stage = unname(code_stage[names(code_files)]),
  md5 = unname(tools::md5sum(code_files)),
  stringsAsFactors = FALSE
)
parameter_values <- vapply(cfg, function(x) paste(x, collapse = ","), character(1))
parameter_values <- c(software_version=mosaic_version, parameter_values,
                      cache_root_resolved=normalizePath(cache_root, winslash="/", mustWork=FALSE))

hash_directory <- function(path) {
  if (!dir.exists(path)) return(NA_character_)
  files <- list.files(path, recursive=TRUE, full.names=TRUE, all.files=TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) return("EMPTY")
  root_norm <- normalizePath(path, winslash="/", mustWork=TRUE)
  file_norm <- normalizePath(files, winslash="/", mustWork=TRUE)
  rel <- substring(file_norm, nchar(root_norm) + 2L)
  ord <- order(rel)
  rel <- rel[ord]; file_norm <- file_norm[ord]
  rows <- paste(rel, unname(tools::md5sum(file_norm)), sep="\t")
  tmp <- tempfile("mosaic_signature_", fileext=".txt")
  on.exit(unlink(tmp), add=TRUE)
  writeLines(rows, tmp, useBytes=TRUE)
  unname(tools::md5sum(tmp))
}

build_data_signatures <- function() {
  components <- c("input_counts", "metadata")
  stages <- c(1L, 1L)
  paths <- c(cfg$input, cfg$metadata)
  md5 <- unname(tools::md5sum(paths))
  if (nzchar(cfg$variable_types_file)) {
    components <- c(components, "variable_types")
    stages <- c(stages, 1L)
    paths <- c(paths, cfg$variable_types_file)
    md5 <- c(md5, unname(tools::md5sum(cfg$variable_types_file)))
  }
  components <- c(components, "annotation_cache")
  stages <- c(stages, 2L)
  paths <- c(paths, normalizePath(cache_root, winslash="/", mustWork=FALSE))
  md5 <- c(md5, hash_directory(cache_root))
  data.frame(component=components, stage=stages, path=paths, md5=md5, stringsAsFactors=FALSE)
}

current_data_signatures <- build_data_signatures()

build_environment_signature <- function() {
  packages <- unique(c(
    "limma", "edgeR", "variancePartition", "BiocParallel", "clusterProfiler",
    "AnnotationDbi", "GO.db", "KEGGREST", "ComplexHeatmap", "GOSemSim",
    "fgsea", "circlize", "multcompView", "openxlsx", "writexl", "jsonlite",
    "performance", cfg$orgdb
  ))
  versions <- vapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) return("NOT_INSTALLED")
    as.character(utils::packageVersion(pkg))
  }, character(1))
  bioc_version <- if (requireNamespace("BiocManager", quietly=TRUE)) {
    as.character(BiocManager::version())
  } else "NOT_INSTALLED"
  data.frame(
    component=c("R", "Bioconductor", packages),
    value=c(paste(R.version$major, R.version$minor, sep="."), bioc_version, unname(versions)),
    stringsAsFactors=FALSE
  )
}

current_environment_signature <- build_environment_signature()

redo_start_step <- 1L
redo_changed <- character(0)
redo_code_changed <- character(0)
redo_data_changed <- character(0)
if (redo_requested) {
  old_parameter_file <- file.path(dirs$config, "run_parameters.tsv")
  old_values <- character(0)
  if (file.exists(old_parameter_file)) {
    old_df <- tryCatch(
      read.delim(old_parameter_file, stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "", comment.char = ""),
      error = function(e) NULL
    )
    if (!is.null(old_df) && all(c("parameter", "value") %in% names(old_df))) {
      old_values <- setNames(as.character(old_df$value), as.character(old_df$parameter))
    }
  }

  old_code_file <- file.path(dirs$config, "code_signatures.tsv")
  old_code <- NULL
  if (file.exists(old_code_file)) {
    old_code <- tryCatch(
      read.delim(old_code_file, stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "", comment.char = ""),
      error = function(e) NULL
    )
    if (is.null(old_code) || !all(c("component", "stage", "md5") %in% names(old_code))) old_code <- NULL
  }

  old_data_file <- file.path(dirs$config, "data_signatures.tsv")
  old_data <- NULL
  if (file.exists(old_data_file)) {
    old_data <- tryCatch(
      read.delim(old_data_file, stringsAsFactors = FALSE, check.names = FALSE,
                 quote = "", comment.char = ""),
      error = function(e) NULL
    )
    if (is.null(old_data) || !all(c("component", "stage", "path", "md5") %in% names(old_data))) old_data <- NULL
  }

  old_environment_file <- file.path(dirs$config, "environment_signature.tsv")
  old_environment <- NULL
  if (file.exists(old_environment_file)) {
    old_environment <- tryCatch(
      read.delim(old_environment_file, stringsAsFactors=FALSE, check.names=FALSE,
                 quote="", comment.char=""),
      error=function(e) NULL
    )
    if (is.null(old_environment) || !all(c("component", "value") %in% names(old_environment))) old_environment <- NULL
  }

  current_values <- parameter_values
  common_names <- union(names(current_values), names(old_values))
  changed <- common_names[vapply(common_names, function(nm) {
    old <- if (nm %in% names(old_values)) old_values[[nm]] else NA_character_
    cur <- if (nm %in% names(current_values)) current_values[[nm]] else NA_character_
    !identical(old, cur)
  }, logical(1))]
  changed <- setdiff(changed, c("clean", "keep_configured", "output"))
  redo_changed <- changed

  stage1_parameters <- c(
    "input", "metadata", "formula", "contrast_factor", "reference",
    "filter_column", "filter_value", "variable_types_file", "interaction_factors", "cores"
  )
  stage2_parameters <- c(
    "cache_dir_option", "cache_root_resolved", "organism_name", "orgdb", "key_type", "symbol_column",
    "kegg_org", "yeast_id_style", "rank_by", "min_gs", "max_gs", "padj_cut",
    "padj_method", "run_go", "run_kegg"
  )
  stage3_parameters <- c("sig_threshold")
  stage4_parameters <- c(
    "filter_mode", "similarity_reduction", "similarity_threshold",
    "similarity_metric", "winner_score", "max_nes_scale"
  )

  if (!length(old_values)) {
    redo_start_step <- 1L
    message("--redo: run parameters were not found; restarting from Step 1.")
  } else if (any(changed %in% stage1_parameters)) {
    redo_start_step <- 1L
  } else if (any(changed %in% stage2_parameters)) {
    redo_start_step <- 2L
  } else if (any(changed %in% stage3_parameters)) {
    redo_start_step <- 3L
  } else {
    redo_start_step <- 4L
  }

  code_start_step <- 4L
  if (is.null(old_code)) {
    redo_code_changed <- current_code_signatures$component
    code_start_step <- 1L
  } else {
    old_md5 <- setNames(as.character(old_code$md5), as.character(old_code$component))
    current_md5 <- setNames(as.character(current_code_signatures$md5), as.character(current_code_signatures$component))
    all_components <- union(names(current_md5), names(old_md5))
    redo_code_changed <- all_components[vapply(all_components, function(nm) {
      old <- if (nm %in% names(old_md5)) old_md5[[nm]] else NA_character_
      cur <- if (nm %in% names(current_md5)) current_md5[[nm]] else NA_character_
      !identical(old, cur)
    }, logical(1))]
    if (length(redo_code_changed)) {
      affected <- current_code_signatures$stage[match(redo_code_changed, current_code_signatures$component)]
      affected[is.na(affected)] <- 1L
      code_start_step <- min(affected)
    }
  }
  redo_start_step <- min(redo_start_step, code_start_step)

  data_start_step <- 4L
  if (is.null(old_data)) {
    redo_data_changed <- current_data_signatures$component
    data_start_step <- 1L
  } else {
    old_sig <- setNames(as.character(old_data$md5), as.character(old_data$component))
    cur_sig <- setNames(as.character(current_data_signatures$md5), as.character(current_data_signatures$component))
    all_data_components <- union(names(cur_sig), names(old_sig))
    changed_data_components <- all_data_components[vapply(all_data_components, function(nm) {
      old <- if (nm %in% names(old_sig)) old_sig[[nm]] else NA_character_
      cur <- if (nm %in% names(cur_sig)) cur_sig[[nm]] else NA_character_
      !identical(old, cur)
    }, logical(1))]
    redo_data_changed <- changed_data_components
    if (length(changed_data_components)) {
      stage_lookup <- setNames(current_data_signatures$stage, current_data_signatures$component)
      old_stage_lookup <- setNames(old_data$stage, old_data$component)
      affected <- vapply(changed_data_components, function(nm) {
        if (nm %in% names(stage_lookup)) stage_lookup[[nm]]
        else if (nm %in% names(old_stage_lookup)) old_stage_lookup[[nm]]
        else 1L
      }, integer(1))
      data_start_step <- min(affected)
    }
  }
  redo_start_step <- min(redo_start_step, data_start_step)

  environment_changed <- character(0)
  if (is.null(old_environment)) {
    environment_changed <- current_environment_signature$component
    redo_start_step <- 1L
  } else {
    old_env <- setNames(as.character(old_environment$value), as.character(old_environment$component))
    cur_env <- setNames(as.character(current_environment_signature$value), as.character(current_environment_signature$component))
    all_env_components <- union(names(cur_env), names(old_env))
    environment_changed <- all_env_components[vapply(all_env_components, function(nm) {
      old <- if (nm %in% names(old_env)) old_env[[nm]] else NA_character_
      cur <- if (nm %in% names(cur_env)) cur_env[[nm]] else NA_character_
      !identical(old, cur)
    }, logical(1))]
    if (length(environment_changed)) redo_start_step <- 1L
  }

  de_ready <- file.exists(file.path(dirs$de, "contrast_manifest.tsv")) &&
    length(list.files(dirs$de, pattern = "[.]txt$", full.names = TRUE)) > 0L
  gsea_status_file <- file.path(dirs$gsea, "GSEA_status.tsv")
  gsea_ready <- FALSE
  if (de_ready && file.exists(gsea_status_file)) {
    cm_redo <- tryCatch(read.delim(file.path(dirs$de, "contrast_manifest.tsv"), stringsAsFactors=FALSE,
                                   check.names=FALSE, quote="", comment.char=""), error=function(e) NULL)
    gs_redo <- tryCatch(read.delim(gsea_status_file, stringsAsFactors=FALSE,
                                   check.names=FALSE, quote="", comment.char=""), error=function(e) NULL)
    if (!is.null(cm_redo) && !is.null(gs_redo) && "label" %in% names(cm_redo) && "contrast" %in% names(gs_redo)) {
      expected_contrasts <- unique(as.character(cm_redo$label))
      ok_status <- c("SUCCESS", "SUCCESS_EMPTY")
      status_complete <- all(expected_contrasts %in% as.character(gs_redo$contrast))
      type_complete <- TRUE
      if (status_complete && cfg$run_go) {
        type_complete <- type_complete && "GO_BP_status" %in% names(gs_redo) &&
          all(gs_redo$GO_BP_status[match(expected_contrasts, gs_redo$contrast)] %in% ok_status) &&
          all(file.exists(file.path(dirs$gsea, paste0(expected_contrasts, "_GO_BP.tsv"))))
      }
      if (status_complete && cfg$run_kegg) {
        type_complete <- type_complete && "KEGG_status" %in% names(gs_redo) &&
          all(gs_redo$KEGG_status[match(expected_contrasts, gs_redo$contrast)] %in% ok_status) &&
          all(file.exists(file.path(dirs$gsea, paste0(expected_contrasts, "_KEGG.tsv"))))
      }
      gsea_ready <- status_complete && type_complete
    }
  }
  requested_summaries <- character(0)
  if (cfg$run_go) requested_summaries <- c(requested_summaries, "GO_BP_treatment_summary.tsv")
  if (cfg$run_kegg) requested_summaries <- c(requested_summaries, "KEGG_treatment_summary.tsv")
  summaries_ready <- length(requested_summaries) > 0L &&
    all(file.exists(file.path(dirs$gsea, requested_summaries)))

  if (redo_start_step > 1L && !de_ready) redo_start_step <- 1L
  if (redo_start_step > 2L && !gsea_ready) redo_start_step <- 2L
  if (redo_start_step > 3L && !summaries_ready) redo_start_step <- 3L

  changed_label <- if (length(redo_changed)) paste(redo_changed, collapse = ", ") else "none"
  code_changed_label <- if (length(redo_code_changed)) paste(redo_code_changed, collapse = ", ") else "none"
  data_changed_label <- if (length(redo_data_changed)) paste(redo_data_changed, collapse = ", ") else "none"
  cat(sprintf("\n--redo enabled: changed parameter(s): %s\n", changed_label))
  cat(sprintf("--redo enabled: changed pipeline component(s): %s\n", code_changed_label))
  environment_changed_label <- if (length(environment_changed)) paste(environment_changed, collapse=", ") else "none"
  cat(sprintf("--redo enabled: changed input/cache component(s): %s\n", data_changed_label))
  cat(sprintf("--redo enabled: changed R/Bioconductor component(s): %s\n", environment_changed_label))
  cat(sprintf("Resuming from Step %d of 4.\n", redo_start_step))

  if (redo_start_step <= 1L) {
    unlink(c(dirs$de, dirs$gsea, dirs$final), recursive = TRUE, force = TRUE)
    unlink(list.files(dirs$logs, pattern = "^(01|02|03|04)_", full.names = TRUE), force = TRUE)
  } else if (redo_start_step == 2L) {
    unlink(c(dirs$gsea, dirs$final), recursive = TRUE, force = TRUE)
    unlink(list.files(dirs$logs, pattern = "^(02|03|04)_", full.names = TRUE), force = TRUE)
  } else if (redo_start_step == 3L) {
    unlink(file.path(dirs$gsea, c("GO_BP_treatment_summary.tsv", "KEGG_treatment_summary.tsv")), force = TRUE)
    unlink(dirs$final, recursive = TRUE, force = TRUE)
    unlink(list.files(dirs$logs, pattern = "^(03|04)_", full.names = TRUE), force = TRUE)
  } else {
    unlink(dirs$final, recursive = TRUE, force = TRUE)
    unlink(file.path(dirs$logs, "04_final_reports.log"), force = TRUE)
  }
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

writeLines(c(
  "# Organism settings",
  sprintf("organism_name <- %s", r_string(cfg$organism_name)),
  sprintf("orgdb <- %s", r_string(cfg$orgdb)),
  sprintf("kegg_org <- %s", r_string(cfg$kegg_org)),
  sprintf("yeast_gene_id_style <- %s", r_string(cfg$yeast_id_style)),
  sprintf("key_type <- %s", r_string(cfg$key_type)),
  sprintf("symbol_column <- %s", r_string(cfg$symbol_column)),
  "is_sgd_yeast <- identical(orgdb, \"org.Sc.sgd.db\")",
  "map_gene_labels <- !identical(key_type, symbol_column)",
  "enable_sgd_gene_links <- is_sgd_yeast",
  "enable_sgd_pathway_lookup <- is_sgd_yeast"
), file.path(dirs$config, "organism_config.R"), useBytes = TRUE)

replace_assignment <- function(text, variable, expression) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  pat <- paste0("^[[:space:]]*", variable, "[[:space:]]*<-")
  idx <- grep(pat, lines)
  if (!length(idx)) {
    stop(sprintf("Could not configure variable '%s' in template", variable), call. = FALSE)
  }
  lines[idx[1]] <- sprintf("%s <- %s", variable, expression)
  paste(lines, collapse = "\n")
}
configure <- function(src, dest, assignments, literal_replacements = list()) {
  txt <- paste(readLines(src, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (nm in names(assignments)) txt <- replace_assignment(txt, nm, assignments[[nm]])
  if (length(literal_replacements)) {
    for (old in names(literal_replacements)) txt <- gsub(old, literal_replacements[[old]], txt, fixed = TRUE)
  }
  writeLines(txt, dest, useBytes = TRUE)
}

int_factor_expr <- if (tolower(cfg$interaction_factors) == "auto") {
  "NULL"
} else if (tolower(cfg$interaction_factors) %in% c("none", "no", "false", "0", "")) {
  "character(0)"
} else {
  r_vector(trimws(strsplit(cfg$interaction_factors, ",", fixed = TRUE)[[1]]))
}
max_nes_expr <- if (!nzchar(cfg$max_nes_scale)) "NULL" else as.character(as.numeric(cfg$max_nes_scale))

configured <- c(
  step1 = file.path(dirs$internal, "step1_de_configured.R"),
  step2 = file.path(dirs$internal, "step2_gsea_configured.R"),
  step3 = file.path(dirs$internal, "step3_summary_configured.R"),
  step4 = file.path(dirs$internal, "step4_reports_configured.R")
)

configure(file.path(template_dir, "step1_de.R"), configured[["step1"]], list(
  counts_file = r_string(cfg$input), metadata_file = r_string(cfg$metadata),
  formula_str = r_string(cfg$formula), contrast_factor = r_string(cfg$contrast_factor),
  reference_level = r_string(cfg$reference), output_dir = r_string(dirs$de),
  filter_column = r_nullable_string(cfg$filter_column), filter_value = r_nullable_string(cfg$filter_value),
  variable_types_file = r_nullable_string(cfg$variable_types_file),
  interaction_factors = int_factor_expr, sig_threshold = as.character(cfg$sig_threshold),
  cores = as.character(cfg$cores)
))

configure(file.path(template_dir, "step2_gsea.R"), configured[["step2"]], list(
  results_dir = r_string(dirs$de), output_dir = r_string(dirs$gsea),
  kegg_cache = r_string(dirs$kegg_cache), go_cache = r_string(dirs$go_cache),
  rank_by = r_string(cfg$rank_by), organism_config_file = r_string(file.path(dirs$config, "organism_config.R")),
  min_gs = as.character(cfg$min_gs), max_gs = as.character(cfg$max_gs),
  padj_cut = as.character(cfg$padj_cut), padj_method = r_string(cfg$padj_method),
  sig_threshold = as.character(cfg$sig_threshold),
  run_go_bp = if (cfg$run_go) "TRUE" else "FALSE", run_kegg = if (cfg$run_kegg) "TRUE" else "FALSE"
))

configure(file.path(template_dir, "step3_summary.R"), configured[["step3"]], list(
  gsea_dir = r_string(dirs$gsea), results_dir = r_string(dirs$de),
  sig_threshold = as.character(cfg$sig_threshold)
))

configure(file.path(template_dir, "step4_reports.R"), configured[["step4"]], list(
  gsea_dir = r_string(dirs$gsea), results_dir = r_string(dirs$de), output_dir = r_string(dirs$final),
  organism_config_file = r_string(file.path(dirs$config, "organism_config.R")),
  reference_level = r_string(cfg$reference), contrast_factor = r_string(cfg$contrast_factor),
  sig_threshold = as.character(cfg$sig_threshold),
  do_similarity_reduction = if (cfg$similarity_reduction) "TRUE" else "FALSE",
  similarity_threshold = as.character(cfg$similarity_threshold), similarity_metric = r_string(cfg$similarity_metric),
  max_nes_scale = max_nes_expr, filter_mode = r_string(cfg$filter_mode), winner_score = r_string(cfg$winner_score),
  sgd_cache_file = r_string(sgd_cache_file)
), literal_replacements = setNames(list(r_string(file.path(dirs$kegg_cache, paste0(cfg$kegg_org, "_term2gene.tsv")))),
                                    'file.path("kegg_cache", paste0(kegg_org, "_term2gene.tsv"))'))

manifest <- data.frame(
  parameter = names(parameter_values),
  value = unname(parameter_values),
  stringsAsFactors = FALSE
)
pending_parameter_file <- file.path(dirs$config, "run_parameters.pending.tsv")
write.table(manifest, pending_parameter_file, sep = "\t", quote = FALSE, row.names = FALSE)
pending_code_file <- file.path(dirs$config, "code_signatures.pending.tsv")
write.table(current_code_signatures, pending_code_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(normalizePath(dirs$cache, winslash="/", mustWork=FALSE), file.path(dirs$config, "CACHE_LOCATION.txt"))

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript was not found in PATH. Activate the supplied Conda environment.", call. = FALSE)

run_step <- function(number, label, script, working_dir = dirs$root) {
  log_file <- normalizePath(file.path(dirs$logs, sprintf("%02d_%s.log", number, label)), winslash = "/", mustWork = FALSE)
  script_abs <- normalizePath(script, winslash = "/", mustWork = TRUE)
  working_dir <- normalizePath(working_dir, winslash = "/", mustWork = TRUE)
  cat(sprintf("\n[%d/4] %s\n", number, label))
  flush.console()

  runner <- tempfile(sprintf("mosaic_stage_%02d_", number), fileext = ".sh")
  runner_lines <- c(
    "#!/usr/bin/env bash",
    "set -o pipefail",
    sprintf("cd %s", shQuote(working_dir)),
    sprintf("%s %s 2>&1 | tee %s", shQuote(rscript), shQuote(script_abs), shQuote(log_file)),
    "exit ${PIPESTATUS[0]}"
  )
  writeLines(runner_lines, runner, useBytes = TRUE)
  Sys.chmod(runner, mode = "0755")
  on.exit(unlink(runner), add = TRUE)

  status <- system2("bash", args = shQuote(runner), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed. Full log: %s", label, log_file), call. = FALSE)
  }
  cat("  done\n")
}

find_summary_files <- function() {
  list.files(dirs$root, pattern = "^(GO_BP|KEGG)_treatment_summary[.]tsv$",
             recursive = TRUE, full.names = TRUE)
}

ensure_summary_files <- function() {
  requested <- character(0)
  if (cfg$run_go) requested <- c(requested, "GO_BP_treatment_summary.tsv")
  if (cfg$run_kegg) requested <- c(requested, "KEGG_treatment_summary.tsv")
  expected <- file.path(dirs$gsea, requested)
  missing <- expected[!file.exists(expected)]
  if (length(missing)) {
    discovered <- find_summary_files()
    for (target in missing) {
      hit <- discovered[basename(discovered) == basename(target)]
      if (length(hit)) file.copy(hit[[1]], target, overwrite = TRUE)
    }
  }
  missing <- expected[!file.exists(expected)]
  if (length(missing)) {
    gsea_files <- list.files(dirs$gsea, recursive = TRUE, full.names = TRUE)
    preview <- if (length(gsea_files)) paste(head(gsea_files, 30L), collapse = "\n  ") else "<none>"
    stop(sprintf(
      "Summary step completed but required summary file(s) are missing: %s\nFiles currently in 02_gsea/:\n  %s\nCheck logs/03_summary.log.",
      paste(basename(missing), collapse = ", "), preview
    ), call. = FALSE)
  }
}

start_time <- Sys.time()
status <- tryCatch({
  start_step <- if (redo_requested) redo_start_step else 1L
  if (start_step <= 1L) run_step(1L, "differential_expression", configured[["step1"]], dirs$root)
  if (start_step <= 2L) run_step(2L, "gsea", configured[["step2"]], dirs$root)
  if (start_step <= 3L) run_step(3L, "summary", configured[["step3"]], dirs$gsea)
  ensure_summary_files()
  run_step(4L, "final_reports", configured[["step4"]], dirs$gsea)
  TRUE
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(dirs$logs, "PIPELINE_FAILED.txt"))
  message("\nPIPELINE FAILED: ", conditionMessage(e))
  FALSE
})

if (!status) quit(save = "no", status = 1L)

outputs <- list.files(dirs$final, recursive = TRUE, full.names = TRUE)
report_outputs <- outputs[grepl("[.](png|html|xlsx)$", outputs, ignore.case = TRUE)]
if (!length(report_outputs)) {
  warning("Pipeline completed successfully, but no PNG, HTML, or XLSX enrichment report was produced. This can occur when no terms remain after filtering; see logs/04_final_reports.log for details.", call. = FALSE)
}
invisible(file.copy(pending_parameter_file, file.path(dirs$config, "run_parameters.tsv"), overwrite = TRUE))
unlink(pending_parameter_file, force = TRUE)
invisible(file.copy(pending_code_file, file.path(dirs$config, "code_signatures.tsv"), overwrite = TRUE))
unlink(pending_code_file, force = TRUE)
write.table(build_data_signatures(), file.path(dirs$config, "data_signatures.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(build_environment_signature(), file.path(dirs$config, "environment_signature.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(dirs$config, "sessionInfo.txt"), useBytes=TRUE)
summary_file <- file.path(dirs$root, "README_RESULTS.txt")
writeLines(c(
  sprintf("MOSAIC-GSEA %s completed successfully.", mosaic_version),
  sprintf("Started: %s", format(start_time)),
  sprintf("Finished: %s", format(Sys.time())),
  "",
  "Directory layout:",
  "  01_differential_expression/  limma/dream contrasts, manifests, and gene-wise variance fractions",
  "  02_gsea/                     GO BP and KEGG GSEA tables and summaries",
  "  03_final_reports/            final enrichment reports and variance-partition summary plot",
  sprintf("  %s/", annotation_dir_name),
  "                               reusable GO, KEGG, and SGD pathway data",
  "  config/                      organism and run parameters",
  "  logs/                        one log per pipeline stage",
  "",
  "Final files:",
  paste0("  ", basename(outputs))
), summary_file)

if (!cfg$keep_configured) unlink(dirs$internal, recursive = TRUE, force = TRUE)
cat("\nCompleted successfully.\n")
cat(sprintf("Reports : %s\n", dirs$final))
cat(sprintf("Summary : %s\n", summary_file))
cat(sprintf("Logs    : %s\n", dirs$logs))
