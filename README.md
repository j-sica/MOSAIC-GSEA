# MOSAIC-GSEA

**Version 1.0.0**

MOSAIC-GSEA is a command-line workflow for analysing functional responses in multifactorial bulk RNA-seq experiments. From a count matrix, sample metadata and a model formula, it performs differential expression analysis, builds model-based contrasts and runs ranked GSEA for GO Biological Process and KEGG pathways. Results are presented in summary heatmaps, interactive HTML reports and Excel workbooks.

Fixed-effect models use limma-voom; mixed models use dream. Depending on the design, the workflow reports treatment effects, pairwise comparisons and two- or three-factor interactions. It also provides gene-wise variance summaries and tools for exploring leading-edge genes.

## Quick start

MOSAIC-GSEA runs on 64-bit Linux. Clone the complete GitHub repository:

```bash
git clone https://github.com/j-sica/MOSAIC-GSEA.git
cd MOSAIC-GSEA
```

Create and activate the supplied Conda environment:

```bash
conda env create -f environment.yml
conda activate mosaic-gsea
chmod +x mosaic-gsea
```

Run the first example:

```bash
./mosaic-gsea \
  -i test_data/feature_counts_ex_1.txt \
  -m test_data/metadata_ex1.tsv \
  -f '~ strain + treatment' \
  -c treatment \
  -r YPD \
  -o example_1
```

Here, `-f` defines the R model, `-c` selects the factor to compare, `-r` selects its reference level, and `-o` sets the output directory. For the complete command-line interface, run `./mosaic-gsea --help`.

### Run from any directory

To call `mosaic-gsea` without specifying its installation path, add the package directory to your Bash `PATH` (change the path if necessary):

```bash
printf '\nexport PATH="$HOME/MOSAIC-GSEA:$PATH"\n' >> ~/.bashrc
source ~/.bashrc
mosaic-gsea --version
```

You can then launch an analysis from a project directory:

```bash
mosaic-gsea -i counts.tsv -m metadata.tsv -f '~ treatment' -c treatment -r control -o results
```

Relative input and output paths refer to the directory where you run the command. MOSAIC-GSEA finds its scripts, templates and bundled annotations relative to the installation directory, including when launched through `PATH`. By default, it prefers the installed annotation cache; use `--cache-dir DIR` to choose another cache parent directory.

## Input files

Provide two tabular files:

- **Count matrix:** gene or feature identifiers in the first column, followed by one column per sample with raw counts.
- **Metadata:** sample identifiers in the first column, followed by the variables referenced in the model formula. Sample names must match the count matrix.

Model variables cannot contain missing values. Metadata names and factor levels must not have surrounding spaces or contain the reserved strings `_vs_`, `_to_`, `__`, `×` or filename-reserved characters; an ordinary `x` is allowed.

Use unique gene identifiers whenever possible. Duplicate count-matrix labels can be analysed as separate rows, but identifiers mapping to the same annotation are combined into a single GSEA ranking entry. The original names and any such resolution are documented in the output files.

### Metadata types

Text and integer-only metadata columns are treated as factors by default; columns containing decimal numbers are treated as numeric. To override this, supply a tab- or space-separated file with `--variable-types`, for example:

```text
variable type
age number
batch factor
time number
```

Continuous variables can be included as covariates, but the primary contrast factor (`-c`) must be categorical. Effective metadata types are saved in `01_differential_expression/metadata_variable_types.tsv`.

## Example designs

The package includes data for three common analyses.

**1. Additive model — treatment effects adjusted for strain**

```bash
mosaic-gsea -i test_data/feature_counts_ex_1.txt -m test_data/metadata_ex1.tsv \
  -f '~ strain + treatment' -c treatment -r YPD -o example_1
```

**2. Two-factor interaction — does oxygen change the treatment response?**

```bash
mosaic-gsea -i test_data/feature_counts_ethanol_ex2-3.txt -m test_data/metadata_ex2-3.tsv \
  -f '~ strain + treatment * oxygen' -c treatment -r YPD -o example_2
```

**3. Three-factor interaction — does the treatment response depend on strain and oxygen together?**

```bash
mosaic-gsea -i test_data/feature_counts_ethanol_ex2-3.txt -m test_data/metadata_ex2-3.tsv \
  -f '~ strain * treatment * oxygen' -c treatment -r YPD -o example_3
```

Interacting categorical factors are detected automatically. Three-factor interaction reporting requires a complete 2 × 2 × 2 design. More complex interaction structures are not supported by the automatic reporting workflow.

MOSAIC-GSEA automatically constructs the relevant model-based contrasts and records their definitions and coefficients in the analysis files.

## Results

Each analysis creates:

```text
results/
├── 01_differential_expression/
├── 02_gsea/
├── 03_final_reports/
├── config/
├── logs/
└── README_RESULTS.txt
```

The final reports contain GO BP and KEGG heatmaps, interactive HTML files and Excel workbooks. The HTML reports let you inspect enrichment statistics, leading-edge genes, volcano plots and overlaps between conditions. For *S. cerevisiae*, relevant SGD gene and biochemical-pathway links are also shown.

MOSAIC-GSEA filters lowly expressed genes with `edgeR::filterByExpr()` and normalizes library composition using TMM. Fixed-effect models are fitted with limma-voom and mixed models with dream. Ranked GSEA uses the corresponding signed gene-level statistic (by default, moderated `t` for limma or standardized `z` for dream), with Benjamini–Hochberg adjustment **within each contrast**.

### Interpretability and GO redundancy

The optional **model-aware interpretability filter** keeps significantly enriched terms whose patterns across conditions are supported by the corresponding model contrasts. For instance, if a GO term is enriched in ethanol under aerobiosis but not anaerobiosis, an apparent difference is not treated as an interaction unless the model supports that difference. The filter uses statistical significance tests of the relevant model-derived contrasts to assess whether enrichment patterns are supported; excluded terms are retained in the Excel reports.

For GO Biological Process, **leading-edge redundancy reduction** groups terms whose enrichment is driven by overlapping sets of genes. Similarity is computed from leading-edge genes across conditions (default: combined Jaccard and overlap similarity; threshold 0.375). One representative is displayed per group in the heatmap, with the other terms still accessible in the reports. Neither step changes the underlying GSEA statistics.

Gene-wise variance estimates are reported separately from enrichment results. For fixed-effect models, the report uses **hierarchical Shapley attribution of voom-weighted model R-squared**: all fixed variables (including categorical factors with multiple levels) are compared through nested models with the same samples and weights. An interaction enters a nested model only when all its participating variables are present; its explanatory contribution is thereby distributed among those variables rather than shown as a separate interaction component. The resulting factor contributions sum to the full weighted model R-squared, with the remaining fraction reported as residual variance, independently of factor reference levels. This exact calculation enumerates subsets of up to 10 fixed variables. Mixed models with categorical fixed effects continue to use aggregate fixed, random and residual components based on marginal and conditional R-squared. These measures are descriptive, not additional significance tests.

Logs include warnings from each processing step. `README_RESULTS.txt`, the model and contrast manifests, and the saved settings document the analysis. To rerun after changing inputs or settings, use the same output directory with `--redo`; valid upstream steps are reused when possible.

## Organisms and annotations

The default organism is *Saccharomyces cerevisiae*. GO analysis requires a compatible Bioconductor `OrgDb` package; MOSAIC-GSEA selects one for supported organisms or accepts `--orgdb` when manual selection is needed. KEGG uses the corresponding organism code, which can be set with `--kegg-org`; use `--run-kegg false` when KEGG is unavailable or should not be queried.

Gene identifiers must match the annotation database. **`--key-type auto` is the default**: MOSAIC-GSEA samples input gene IDs and tests their actual coverage against supported identifier namespaces in the selected `OrgDb`. In yeast, if `ORF` and `ENSEMBL` have similarly high coverage of systematic SGD IDs (for example, `YBR160W`, also known as `CDC28`), `ORF` is preferred. Common names such as `CDC28` can be detected as `COMMON`. Reports normally display familiar gene names where an appropriate annotation is available. `--yeast-id-style systematic|common` remains supported as a tie-breaking preference for backward compatibility, **not** as a setting that overrides an explicitly supplied `--key-type`.

Specify `--key-type ORF`, `--key-type ENSEMBL`, `--key-type ENTREZID`, `--key-type SYMBOL` or any other key type supported by the selected `OrgDb` to override auto-detection. If auto-detection has low coverage (below 50%) or similarly good matches in multiple identifier namespaces without a yeast-specific tie-break, the program stops and asks for an explicit choice. An explicit choice is always respected, but input coverage below 80% generates a warning. Auto-detection checks all unique input gene IDs; it cannot guarantee all input IDs have annotations. Mapping coverage is also reported by the GSEA step. **Changing the identifier namespace can change GO/KEGG enrichment results**: re-run analyses after changing the key type. GO and KEGG mapping caches are isolated by key type; the bundled yeast GO cache is for the legacy `ENSEMBL` namespace, so the first `ORF` run generates and caches its own GO annotations using `org.Sc.sgd.db`.

The package includes an ENSEMBL-keyed legacy GO annotation cache and an SGD biochemical-pathway cache for *S. cerevisiae*; the default auto-detected `ORF` setting builds a separate GO cache on first use. **KEGG data are not bundled:** on the first run with KEGG enabled, MOSAIC-GSEA downloads the gene–pathway associations and pathway names using the KEGG API, then saves them in a local cache for subsequent runs. **An internet connection is required for this first KEGG run**, for new organisms or identifier types, or if the cache is deleted. KEGG access is subject to its provider's terms (academic access through the public API; commercial users need an appropriate KEGG license). To run GO analysis without accessing KEGG, use `--run-kegg false`.

Valid local annotation caches are reused. If the SGD cache is missing or invalid, MOSAIC-GSEA can rebuild it, potentially requiring access to SGD. SGD-specific links and pathway lookups are omitted for other organisms. Annotation availability varies across species; MOSAIC-GSEA does not transfer annotations between organisms.

## Useful options and limitations

- `--interaction-factors auto` detects relevant interactions; use `none` or an explicit list to override detection.
- `--filter-mode`, `--similarity-metric` and `--similarity-threshold` control interpretation and GO redundancy reporting.
- `--rank-by`, `--min-gs`, `--max-gs`, `--padj-method` and `--sig-threshold` control ranking and enrichment settings.
- `--cores N` sets parallel workers for mixed-model and variance fitting (default: 4). If fewer cores are detected, the program warns and uses the available number.
- `--key-type auto` detects input gene IDs from the selected OrgDb; explicit `--key-type TYPE` overrides it.
- `--cache-dir DIR` selects a different annotation-cache location; `--redo` resumes a changed analysis.

The workflow requires estimable models with residual degrees of freedom. It does not support no-intercept or saturated designs, arbitrary formula functions, or continuous variables as the primary contrast. GO enrichment relies on organism annotations, not de novo gene function prediction. KEGG resources are subject to the provider's access and licensing terms.

For further options, run `mosaic-gsea --help`.

## License

MOSAIC-GSEA source code is released under **GPL-3.0** (see `LICENSE`), including permission for commercial use of the software itself. External annotation data remain subject to their providers' licenses: Gene Ontology Consortium data are available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (source: [Gene Ontology Consortium](https://geneontology.org/)); SGD annotations are credited to the [Saccharomyces Genome Database](https://www.yeastgenome.org/). KEGG data are **not included** in this distribution. Use of the KEGG API and KEGG-derived data is subject to [KEGG's terms](https://www.genome.jp/kegg/legal.html); non-academic users must obtain the appropriate KEGG license. The GPL license for MOSAIC-GSEA does not grant rights to KEGG data.

### Identifier conversion and diagnostic messages

MOSAIC-GSEA checks every input gene ID against the available identifier namespaces of the selected OrgDb and reports coverage by namespace. Auto-detection requires at least 50% coverage, warns when coverage is below 80%, and requests an explicit `--key-type` when competing namespaces cannot be distinguished safely. Explicit key types are validated against OrgDb and warned when coverage is below 80%. KEGG identifiers are mapped using direct and available OrgDb-mediated mappings, reporting ambiguous mappings and conversion coverage. The versioned KEGG identifier-mapping cache avoids silently reusing legacy conversions. Existing enrichment results should be regenerated if identifier mappings change.
