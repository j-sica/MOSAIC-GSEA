# Organism setup

normalize_species_name <- function(x) {
  x <- trimws(as.character(x))
  gsub("[[:space:]]+", " ", x)
}

ensure_biocmanager <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
}

known_orgdb <- c(
  "saccharomyces cerevisiae" = "org.Sc.sgd.db",
  "homo sapiens" = "org.Hs.eg.db",
  "mus musculus" = "org.Mm.eg.db",
  "rattus norvegicus" = "org.Rn.eg.db",
  "drosophila melanogaster" = "org.Dm.eg.db",
  "caenorhabditis elegans" = "org.Ce.eg.db",
  "danio rerio" = "org.Dr.eg.db",
  "arabidopsis thaliana" = "org.At.tair.db"
)

known_kegg <- c(
  "saccharomyces cerevisiae" = "sce",
  "homo sapiens" = "hsa",
  "mus musculus" = "mmu",
  "rattus norvegicus" = "rno",
  "drosophila melanogaster" = "dme",
  "caenorhabditis elegans" = "cel",
  "danio rerio" = "dre",
  "arabidopsis thaliana" = "ath"
)

resolve_kegg_org <- function(organism_name, requested = "auto", run_kegg = TRUE) {
  requested <- trimws(as.character(requested))
  if (nzchar(requested) && tolower(requested) != "auto") return(requested)
  known <- unname(known_kegg[tolower(normalize_species_name(organism_name))])
  if (length(known) == 1L && !is.na(known)) return(known)
  if (!isTRUE(run_kegg)) return("")
  stop(sprintf(
    "KEGG code for '%s' is not in the built-in table. Use --kegg-org CODE, or set --run-kegg false.",
    organism_name
  ), call. = FALSE)
}

orgdb_candidates <- function(organism_name) {
  parts <- strsplit(normalize_species_name(organism_name), " ", fixed = TRUE)[[1]]
  if (length(parts) < 2L) return(character(0))
  code <- paste0(toupper(substr(parts[[1]], 1L, 1L)), tolower(substr(parts[[2]], 1L, 1L)))
  pattern <- paste0("^org\\.", code, "\\..*\\.db$")
  installed <- rownames(installed.packages())
  local_hits <- installed[grepl(pattern, installed)]
  if (length(local_hits)) return(local_hits)

  ensure_biocmanager()
  repos <- BiocManager::repositories()
  ap <- tryCatch(
    available.packages(repos = repos, type = "source"),
    error = function(e) NULL
  )
  if (is.null(ap)) return(character(0))
  pkgs <- rownames(ap)
  pkgs[grepl(pattern, pkgs)]
}

resolve_orgdb <- function(organism_name, requested = "auto") {
  organism_name <- normalize_species_name(organism_name)
  requested <- trimws(as.character(requested))
  if (nzchar(requested) && tolower(requested) != "auto") return(requested)

  known <- unname(known_orgdb[tolower(organism_name)])
  if (length(known) == 1L && !is.na(known)) return(known)

  candidates <- orgdb_candidates(organism_name)
  if (length(candidates) == 1L) return(candidates[[1]])
  if (!length(candidates)) {
    stop(sprintf(
      "No Bioconductor OrgDb package was found automatically for '%s'. Use --orgdb PACKAGE if one exists.",
      organism_name
    ), call. = FALSE)
  }
  stop(sprintf(
    "More than one OrgDb package may match '%s': %s. Choose one with --orgdb PACKAGE.",
    organism_name, paste(candidates, collapse = ", ")
  ), call. = FALSE)
}

ensure_orgdb <- function(package, organism_name = NULL) {
  ensure_biocmanager()

  if (!requireNamespace(package, quietly = TRUE)) {
    cat(sprintf("OrgDb %s is not installed. Installing it now...\n", package))
    ok <- tryCatch({
      BiocManager::install(package, ask = FALSE, update = FALSE)
      TRUE
    }, error = function(e) {
      message(conditionMessage(e))
      FALSE
    })
    if (!ok || !requireNamespace(package, quietly = TRUE)) {
      stop(sprintf(
        "OrgDb package '%s' could not be installed. Check that it exists for your Bioconductor release.",
        package
      ), call. = FALSE)
    }
  }

  suppressPackageStartupMessages(library(package, character.only = TRUE))
  obj <- getExportedValue(package, package)

  if (!is.null(organism_name) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    md <- tryCatch(AnnotationDbi::metadata(obj), error = function(e) NULL)
    if (!is.null(md) && all(c("name", "value") %in% names(md))) {
      hit <- md$value[toupper(md$name) == "ORGANISM"]
      if (length(hit)) {
        expected <- tolower(normalize_species_name(organism_name))
        observed <- tolower(normalize_species_name(hit[[1]]))
        if (nzchar(expected) && nzchar(observed) && !identical(expected, observed)) {
          stop(sprintf(
            "OrgDb '%s' is for '%s', not '%s'.",
            package, hit[[1]], organism_name
          ), call. = FALSE)
        }
      }
    }
  }

  invisible(obj)
}
