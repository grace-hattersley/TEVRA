# TEVRA workflow

``` r

library(TEVRA)
library(SingleCellExperiment)

set.seed(0)
```

## Load data

TEVRA accepts a `SingleCellExperiment` object that has been
pre-processed, including reference alignment, transcript quantification,
QC and cell annotation.

``` r

data(tevra_example)
tevra_example
#> class: SingleCellExperiment 
#> dim: 2000 3430 
#> metadata(1): graphs
#> assays(1): logcounts
#> rownames(2000): Pcmtd1 Sgk3 ... CAAA01147332.1 AC149090.1
#> rowData names(0):
#> colnames(3430): AAACCTGAGTTACCCA-1 AAACCTGCAATGTTGC-1 ...
#>   TTTGTCATCACCGTAA-1 TTTGTCATCGCTAGCG-1
#> colData names(18): orig.ident nCount_RNA ... CTaa_BCR_IGH CTaa_BCR_IGL
#> reducedDimNames(2): PCA UMAP
#> mainExpName: RNA
#> altExpNames(0):
```

## Defining cellular aggregates

Communities of transcriptionally similar cells can be computed upstream
of the TEVRA workflow using any preferred method e.g., *miloR*,
*metacell* or *supercell*. If discrete communities have been defined,
these should be stored as a column in `colData(sce_object)`. If
overlapping neighbourhoods were computed e.g., via *miloR*, these should
be stored as a binary `nhoods` matrix.

We recommend utilising purpose-built softwares for community detection
such as the methods listed above for more robust cell aggregation.
However, TEVRA can also do this natively with create_communities(), a
custom function for community detection, by Leiden partitioning of the
shared-nearest-neighbour graph. The key parameters are:

- `resolution`: determines the granularity of clustering
- `min_cells`: informs the minimum number of cells per community. Cells
  belonging to communities with fewer cells are reassigned to the
  community occupied by the majority of the cell’s nearest neighbours

To establish communities of the desired size, different resolutions can
be evaluated via
[`test_resolution()`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
and visualised via
[`plot_resolution_results()`](https://grace-hattersley.github.io/TEVRA/reference/plot_resolution_results.md)

``` r

res_results <- test_resolution(tevra_example,
                               graph.name = "RNA_snn",
                               resolutions = c(0.5,1,2.5,5,10,20),
                               min_cells = 5)
#> Testing resolutions...
#> Resolution: 0.5
#> Resolution: 1
#> Resolution: 2.5
#> Resolution: 5
#> Resolution: 10
#> Resolution: 20
#> 
#> Resolution testing complete!
head(res_results)
#>   resolution n_communities mean_cells_per_community sd_cells_per_community
#> 1        0.5             9                381.11111             374.370954
#> 2        1.0            14                245.00000              87.088814
#> 3        2.5            27                127.03704              53.868416
#> 4        5.0            52                 65.96154              23.790723
#> 5       10.0            88                 38.97727              15.577007
#> 6       20.0           178                 19.26966               8.606612
#>   median_cells_per_community min_cells_per_community max_cells_per_community
#> 1                      168.0                      83                    1081
#> 2                      237.0                      83                     370
#> 3                      134.0                      27                     244
#> 4                       66.5                      22                     124
#> 5                       37.0                       9                      82
#> 6                       18.0                       5                      47
#>   modularity
#> 1  0.7002931
#> 2  0.7295705
#> 3  0.7035432
#> 4  0.6404798
#> 5  0.5851922
#> 6  0.4821389
```

``` r

plot_resolution_results(res_results, 
                        show_error_bars = TRUE)
```

![](TEVRA_files/figure-html/plot_resolution_results-1.png)

Once the optimal resolution has been selected, communities can be
generated via
[`create_communities()`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
and community sizes evaluated via
[`plot_nhood_sizes()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_sizes.md)

``` r

communities <- create_communities(tevra_example,
                                  graph.name = "RNA_snn",
                                  resolution = 5,
                                  min_cells = 5,
                                  feature.col = "nFeature_RNA",
                                  count.col = "nCount_RNA")
#> [1] "Created 50 communities"
#> [1] "Mean cells per community: 68.6"
communities
#> $sce
#> class: SingleCellExperiment 
#> dim: 2000 3430 
#> metadata(1): graphs
#> assays(1): logcounts
#> rownames(2000): Pcmtd1 Sgk3 ... CAAA01147332.1 AC149090.1
#> rowData names(0):
#> colnames(3430): AAACCTGAGTTACCCA-1 AAACCTGCAATGTTGC-1 ...
#>   TTTGTCATCACCGTAA-1 TTTGTCATCGCTAGCG-1
#> colData names(19): orig.ident nCount_RNA ... CTaa_BCR_IGL community_id
#> reducedDimNames(2): PCA UMAP
#> mainExpName: RNA
#> altExpNames(0):
#> 
#> $communities_summary
#> # A tibble: 50 × 4
#>    community_id n_cells mean_detected mean_total_counts
#>           <dbl>   <int>         <dbl>             <dbl>
#>  1            1      99         1989.             4771.
#>  2            2      49         1989.             4218.
#>  3            3      83         1366.             2199.
#>  4            4      84         2003.             4965.
#>  5            5      80         2139.             5694.
#>  6            6      58         1881.             4352.
#>  7            7      90         2056.             4813.
#>  8            8      74         1813.             4176.
#>  9            9      72         2072.             4912.
#> 10           10      54         1935.             4383.
#> # ℹ 40 more rows

sce_obj <- communities$sce
```

``` r

plot_nhood_sizes(communities$communities_summary)
```

![](TEVRA_files/figure-html/plot_nhood_sizes-1.png)

A community-level UMAP can be generated via aggregating (mean or median)
the consistuent cells’ individual UMAP coordinates per community

``` r

sce_obj <- nhood_umap(sce_obj,
                      nhood_col = "community_id",
                      aggregate_fun = "median",
                      result_name   = "nhood_umap")
#> Neighbourhood UMAP centroids (median) for 50 neighbourhoods stored in metadata(object)[['nhood_umap']].

head(metadata(sce_obj)$nhood_umap)
#>   nhood nhood_size median_UMAP_1 median_UMAP_2
#> 1     1         99    -0.1635164    3.39397471
#> 2     2         49    -6.7208307   -0.48542459
#> 3     3         83    -2.8592288   -2.98583253
#> 4     4         84     5.6871135   -0.40307231
#> 5     5         80     6.3954184   -0.02601494
#> 6     6         58     2.3985686    1.21709789
```

## Aggregate-level analysis

Once cells are aggregated into transcriptionally-defined clusters, these
can be fed into TEVRA’s aggregate-level analysis functions.

### Cell type annotation

Each aggregate can be annotated according to the most common identity
among its constituent cells. Specify one or more `colData` columns in
`annotation_cols`. For each column, the most frequent value across the
aggregate’s cells becomes the aggregate’s annotation, provided it
accounts for a proportion of cells above the defined `threshold`.

``` r

sce_obj <- nhood_annotation(object = sce_obj,
                            annotation_cols = "Bcelltype",
                            nhood_col = "community_id",
                            threshold = 0.5,
                            result_name = "nhood_annotation")
#> Annotation column 'Bcelltype':
#>   49 / 50 neighbourhoods assigned a dominant identity (threshold = 0.50).
#>   1 / 50 neighbourhoods labelled 'mixed' (no identity exceeded threshold = 0.50).
#>   Dominant identity breakdown:
#>     Follicular 1: 12 neighbourhood(s)
#>     Follicular 2: 11 neighbourhood(s)
#>     Follicular 3: 8 neighbourhood(s)
#>     Transitional: 7 neighbourhood(s)
#>     Follicular 4: 5 neighbourhood(s)
#>     B1: 2 neighbourhood(s)
#>     MZ: 2 neighbourhood(s)
#>     Follicular 5: 1 neighbourhood(s)
#>     Plasma: 1 neighbourhood(s)
#> Stored annotation results (50 neighbourhoods x 1 annotation column(s)) in metadata(object)[['nhood_annotation']].

head(metadata(sce_obj)$nhood_annotation)
#>   nhood nhood_size    Bcelltype Bcelltype_dominant_prop
#> 1     1         99 Follicular 4               0.9191919
#> 2     2         49 Transitional               1.0000000
#> 3     3         83       Plasma               1.0000000
#> 4     4         84           B1               0.9761905
#> 5     5         80           B1               0.9750000
#> 6     6         58 Follicular 1               0.7068966
```

To investigate the proportion of cells per aggregate with each identity
without assigning an annotation, run
[`nhood_annotation_props()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md)

### Isotype annotation

[`nhood_annotation()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md)
can additionally be used to annotate isotype usage (e.g., IGHD, IGHM,
IGHA, IGHG) per cellular aggregate.

### Clonal diversity

The clonal repertoire of each aggregate is summarised using
[`nhood_diversity()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md).
This function incorporates a number of different methods to summarise
clonal diversity. The most simple is `clones_per_contig`, which is
always reported in the output and informs the proportion of recovered
contigs that represent a unique clone. The other available metrics are
listed below and must be selected by the user to be calculated:

- `shannon`:
- `norm.entropy`
- `inv.simpson`
- `gini.simpson`
- `chao1`
- `ace`
- `gini`
- `standardised_n_clones`

If clones have not been defined, the cdr3 sequence is used as a
surrogate, although this is not recommended for B cells due to the
potential for significant sequence modification within clones via
somatic hypermutation. Metrics are not computed for aggregates with less
than `min_n_contigs` as diversity estimates are unstable at low sample
sizes. Rarefaction or coverage- based standardisation is recommended due
to the sensitivity of diversity estimates to sampling depth

``` r

sce_obj <- nhood_diversity(sce_obj,
                           nhood_col = "community_id",
                           clone_id_col = "CTstrict_BCR",
                           cdr3_col = "CTnt_BCR",
                           min_n_contigs = 10,
                           metrics = c("shannon", "inv.simpson", "gini.simpson", "chao1"),
                           downsample_method = "rarefaction",
                           n_boots = 100,
                           feature_scale = FALSE,
                           result_name = "nhood_diversity")
#> Warning in nhood_diversity(sce_obj, nhood_col = "community_id", clone_id_col =
#> "CTstrict_BCR", : Rarefaction is not compatible with chao1 or ace estimators.
#> These metrics will be computed on observed counts without rarefaction.
#> _diversity_sd will be NA for these metrics.
#> Contig summary:
#> Cells with contig: 1515 / 3430
#> Neighbourhoods with >= 10 cells with contig: 49 / 50
#> Rarefaction depth: 27
#> Neighbourhoods successfully rarefied (n_contigs > rarefaction depth): 24 / 49
#> Stored diversity results (50 neighbourhoods x 4 metrics) in metadata(object)[['nhood_diversity']].

head(metadata(sce_obj)$nhood_diversity)
#>   nhood nhood_size n_contigs contig_capture_rate n_clones clone_per_contig
#> 1     1         99        43           0.4343434       41        0.9534884
#> 2    10         54        25           0.4629630       23        0.9200000
#> 3    11         68        17           0.2500000       16        0.9411765
#> 4    12         90        26           0.2888889       26        1.0000000
#> 5    13        119        41           0.3445378       37        0.9024390
#> 6    14         87        39           0.4482759       36        0.9230769
#>     coverage diversity_reliable shannon_diversity shannon_diversity_sd
#> 1 0.09302326               TRUE          3.254248           0.03476596
#> 2 0.16000000               TRUE          3.107972                   NA
#> 3 0.11764706               TRUE          2.751667                   NA
#> 4 0.00000000               TRUE          3.258097                   NA
#> 5 0.19512195               TRUE          3.208038           0.04911501
#> 6 0.15384615               TRUE          3.221388           0.04401398
#>   inv.simpson_diversity inv.simpson_diversity_sd gini.simpson_diversity
#> 1              25.29050                 1.110711              0.9605761
#> 2              21.55172                       NA              0.9536000
#> 3              15.21053                       NA              0.9342561
#> 4              26.00000                       NA              0.9615385
#> 5              23.94501                 1.595969              0.9577503
#> 6              24.72817                 1.296804              0.9587929
#>   gini.simpson_diversity_sd chao1_diversity chao1_diversity_sd
#> 1               0.001975085           288.0                 NA
#> 2                        NA            93.0                 NA
#> 3                        NA            68.5                 NA
#> 4                        NA            26.0                 NA
#> 5               0.002571843           142.6                 NA
#> 6               0.002152818           168.0                 NA
```

### Inter-clonal cdr3 convergence

TEVRA provides a novel inter-clonal CDR3 convergence metric computed via
[`nhood_cdr3_convergence()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md).
This reports the similarity of the CDR3 amino acid sequences between
distinct clonotypes, within each aggregate. An aggregate with high CDR3
convergence represents one in which clonotypes are converging onto a
similar CDR3 sequence, and therefore likely respond to the same epitope.
Alternatively, an aggregate with low CDR3 convergence suggests
clonotypes respond to multiple antigens or multiple epitopes of the same
antigen with dissimilar sequences. The function computes this for both
heavy and light chain CDR3 sequences when `chain` = “both”. For more
details on specifically how this metric is computed, refer to the
section at the end of this vignette.

The key metrics are:

- `receptor_type`: for TCRs, similarity scores are computed between each
  unique CDR3 sequence as these represent unique clones. For BCRs,
  somatic hypermutation enables CDR3 sequences of the same clone to have
  different sequences. Hence, a representative sequence is first
  selected per clone as the CDR3 sequence that is the most similar to
  all other CDR3 sequences of that clone
- `min_n_contigs`: the minimum number of contigs required per aggregate.
  Must exceed `min_clones`
- `min_clones`: the minimum number of clones required per aggregate
- `distance_method`: CDR3 sequences can differ by either insertions,
  deletions or substitutions. All insertions or deletions incur a cost
  of 1. Using Levenshtein distances, all substitutions also incur a cost
  of 1. Using blosum- weighted distances, the cost of substitutions is
  weighted according to the biological significance of this substitution
  (see the end of the vignette for more details). A custom substitution
  matrix can also be supplied with unique substitution weightings,
  ranging from -1 (least costly substitution) to 0 (most costly
  substitution)
- `summary_metric`: inter-clonal convergence scores are summarised per
  aggregate either globally by providing the mean or median score or
  locally. Local convergence scores report the mean distance of each
  CDR3 sequence to its k nearest neighbours, averaged per aggregate (see
  the end of the vignette for more details).
- `downsample`: if TRUE, matrices of clonal convergence scores are
  downsampled to the median clonal count over `n_boots` iterations to
  control for variation in clone number across aggregates

``` r

sce_obj <- nhood_cdr3_convergence(object = sce_obj,
                                  nhood_col ="community_id",
                                  receptor_type = "BCR",
                                  chain = "both",
                                  aa_col = c(VDJ = "CTaa_BCR_IGH", VJ = "CTaa_BCR_IGL"),
                                  clone_id_col = "CTstrict_BCR",
                                  min_n_contigs = 10,
                                  min_clones = NULL,
                                  distance_method = "levenshtein", 
                                  summary_metric = "converge_score",
                                  downsample = TRUE,
                                  n_boots = 100,
                                  BPPARAM = BiocParallel::SerialParam(),
                                  result_name = "nhood_cdr3_convergence") 
#> Extracting neighbourhood membership...
#> Filtering to cells with complete contig annotation...
#> Pre-pass: counting clones per neighbourhood...
#>   VDJ: 33 / 50 neighbourhoods will be processed (downsampling depth = 16 representatives) | 17 skipped (< 10 cells with contig annotation)
#>   VJ: 49 / 50 neighbourhoods will be processed (downsampling depth = 26 representatives) | 1 skipped (< 10 cells with contig annotation)
#> Computing CDR3 distances for 49 valid neighbourhoods (of 50 total)...
#>   Processing neighbourhood 5 / 50 (10%)...
#>   Processing neighbourhood 10 / 50 (20%)...
#>   Processing neighbourhood 15 / 50 (30%)...
#>   Processing neighbourhood 20 / 50 (40%)...
#>   Processing neighbourhood 25 / 50 (50%)...
#>   Processing neighbourhood 30 / 50 (60%)...
#>   Processing neighbourhood 35 / 50 (70%)...
#>   Processing neighbourhood 40 / 50 (80%)...
#>   Processing neighbourhood 45 / 50 (90%)...
#>   Processing neighbourhood 50 / 50 (100%)...
#> Collating results...
#> 
#> CDR3 convergence score computation complete.
#>   Receptor type  : BCR
#>   Distance method: levenshtein
#>   Chain(s)       : VDJ, VJ
#>   Summary metric(s): converge_score
#>   VDJ: 33 / 50 neighbourhoods with reliable convergence estimates (0 with NA converge_score due to n_clones <= k)
#>   VJ: 49 / 50 neighbourhoods with reliable convergence estimates (0 with NA converge_score due to n_clones <= k)
#> Stored CDR3 convergence results (50 neighbourhoods x 14 columns) in metadata(object)[['nhood_cdr3_convergence']].

head(metadata(sce_obj)$nhood_cdr3_convergence)
#>   nhood nhood_size VDJ_n_contigs VDJ_contig_capture_rate VDJ_n_clones
#> 1     1         99            24               0.2424242           23
#> 2     2         49            24               0.4897959           23
#> 3     3         83            11               0.1325301           10
#> 4     4         84            28               0.3333333           23
#> 5     5         80            26               0.3250000           25
#> 6     6         58             9               0.1551724            8
#>   VDJ_converge_score_reliable CDR3_converge_score_VDJ
#> 1                        TRUE                0.578668
#> 2                        TRUE                0.576261
#> 3                        TRUE                0.550000
#> 4                        TRUE                0.593832
#> 5                        TRUE                0.583904
#> 6                       FALSE                      NA
#>   CDR3_converge_score_sd_VDJ VJ_n_contigs VJ_contig_capture_rate VJ_n_clones
#> 1                 0.02072197           43              0.4343434          41
#> 2                 0.01761828           29              0.5918367          27
#> 3                         NA           12              0.1445783          11
#> 4                 0.01817885           42              0.5000000          35
#> 5                 0.01668564           40              0.5000000          36
#> 6                         NA           24              0.4137931          23
#>   VJ_converge_score_reliable CDR3_converge_score_VJ CDR3_converge_score_sd_VJ
#> 1                       TRUE               0.691739               0.017611755
#> 2                       TRUE               0.634082               0.003144673
#> 3                       TRUE               0.635000                        NA
#> 4                       TRUE               0.647581               0.012022132
#> 5                       TRUE               0.641388               0.015655762
#> 6                       TRUE               0.656900                        NA
```

To investigate the raw convergence matrices without neighbourhood-level
summarisation, run
[`nhood_cdr3_convergence_matrices()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence_matrices.md).

### Somatic hypermutation

BCR mutation frequency or counts can be calculated per aggregate using
[`nhood_shm()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md).
This function returns either the mean or median mutation metric and its
standard deviation or interquartile range per aggregate.

### Merging multiple aggregate-level analysese

The above functions compute aggregate-level statistics for respective
analyses.
[`merge_nhood_results()`](https://grace-hattersley.github.io/TEVRA/reference/merge_nhood_results.md)
easily aggregates all neighbourhood-level analyses into a single
dataframe, so that each nhood has multiple AIRR metrics that can be
jointly analysed. For example, following merging, users may cluster
neighbourhoods using multiple parameters derived from the above
analyses, such that V(D)J metrics can serve an independent axis for cell
or aggregate identity.

``` r

sce_obj <- merge_nhood_results(object = sce_obj,
                               dataframe_names = c("nhood_umap", "nhood_annotation", 
                                "nhood_diversity", "nhood_cdr3_convergence"),
                               merge_key = "nhood",
                               size_col = "nhood_size",
                               verbose = TRUE,
                               result_name = "merged_nhood_analysis")
#> ======================================================================
#> Merging Neighbourhood Analysis Results
#> ======================================================================
#> Checking available data (4 requested):
#>   nhood_umap: ✓ Found
#>   nhood_annotation: ✓ Found
#>   nhood_diversity: ✓ Found
#>   nhood_cdr3_convergence: ✓ Found
#> Input dimensions:
#>   nhood_umap                          50 rows × 4 columns
#>   nhood_annotation                    50 rows × 4 columns
#>   nhood_diversity                     50 rows × 16 columns
#>   nhood_cdr3_convergence              50 rows × 14 columns
#> Validating neighbourhood consistency...
#> ----------------------------------------------------------------------
#>   ✓ 'nhood_umap' vs 'nhood_annotation': all 50 neighbourhoods shared
#>   ✓ 'nhood_umap' vs 'nhood_annotation': size column ('nhood_size') IDENTICAL across all 50 shared neighbourhoods
#>   ✓ 'nhood_umap' vs 'nhood_diversity': all 50 neighbourhoods shared
#>   ✓ 'nhood_umap' vs 'nhood_diversity': size column ('nhood_size') IDENTICAL across all 50 shared neighbourhoods
#>   ✓ 'nhood_umap' vs 'nhood_cdr3_convergence': all 50 neighbourhoods shared
#>   ✓ 'nhood_umap' vs 'nhood_cdr3_convergence': size column ('nhood_size') IDENTICAL across all 50 shared neighbourhoods
#> ----------------------------------------------------------------------
#> Performing sequential merge on 'nhood' column...
#>   Step 1: 'nhood_umap'  +  'nhood_annotation'
#>   Checking 1 common column(s) between 'nhood_umap' and 'nhood_annotation':
#>     ✓ nhood_size                     IDENTICAL
#>   Removing 1 redundant column(s) from 'nhood_annotation': nhood_size
#>     → 50 rows × 6 columns
#>   Step 2: 'nhood_umap + nhood_annotation'  +  'nhood_diversity'
#>   Checking 1 common column(s) between 'nhood_umap + nhood_annotation' and 'nhood_diversity':
#>     ✓ nhood_size                     IDENTICAL
#>   Removing 1 redundant column(s) from 'nhood_diversity': nhood_size
#>     → 50 rows × 20 columns
#>   Step 3: 'nhood_umap + nhood_annotation + nhood_diversity'  +  'nhood_cdr3_convergence'
#>   Checking 1 common column(s) between 'nhood_umap + nhood_annotation + nhood_diversity' and 'nhood_cdr3_convergence':
#>     ✓ nhood_size                     IDENTICAL
#>   Removing 1 redundant column(s) from 'nhood_cdr3_convergence': nhood_size
#>     → 50 rows × 32 columns
#> Stored merged results (50 rows × 32 columns) in metadata(object)[['merged_nhood_analysis']].
```

## Plot

TEVRA offers a number of visualisation functions. Each neighbourhood is
visualised as a point.

### Neighbourhood UMAP

The aggregated UMAP can be visualised via
[`plot_nhood_umap()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_umap.md).
This function enables, the aggregates to be plotted and to be optionally
overlayed onto the original cell-level UMAP from which the aggregated
coordinates are derived. If the user has also computed aggregate
similarity weightings e.g., as obtained via *metacell*, these can also
be plotted.

``` r

cell_umap_coords <- as.data.frame(reducedDims(sce_obj)$UMAP)
summary(row.names(cell_umap_coords) == row.names(colData(sce_obj)))
#>    Mode    TRUE 
#> logical    3430
cell_umap_coords$community_id <- colData(sce_obj)$community_id

plot_nhood_umap(metadata(sce_obj)$merged_nhood_analysis,
                nhood_col = "nhood",
                x_col = "median_UMAP_1",
                y_col = "median_UMAP_2",
                fill_by = "Bcelltype",
                show_cells = T,
                cell_umap_coords = cell_umap_coords,
                cell_x_col = "umap_1",
                cell_y_col = "umap_2",
                cell_colour_by = "nhood",
                cell_nhood_col = "community_id")
```

![](TEVRA_files/figure-html/plot_nhood_umap-1.png)

### Capture efficiency

Visualises the efficiency of contig capture per neighbourhood. Follows
[`nhood_diversity()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md)

- `capture_rate` demonstrates the proportion of cells per aggregate with
  a recovered contig as a jitter plot with bix plot overlayed. This
  enables visual investigation of the rate of contig recovery per cell
  type
- `count_vs_size` is a scatter plot with a line of best fit overlayed,
  demonstrating the relationship between aggregate size and the number
  of contigs recovered. A positive, linear relationship is expected,
  with the number of contigs increasing with neighbourhood size.

``` r

plot_capture_efficiency(nhood_diversity_df = metadata(sce_obj)$merged_nhood_analysis, 
                        group_by = "Bcelltype") 
#> $capture_rate
```

![](TEVRA_files/figure-html/plot_capture_efficiency-1.png)

    #> 
    #> $count_vs_size
    #> `geom_smooth()` using formula = 'y ~ x'

![](TEVRA_files/figure-html/plot_capture_efficiency-2.png)

### Beeswarm

[`plot_nhood_beeswarm()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_beeswarm.md)
visualises neighbourhood analysis as a beeswarm plot. A plot is provided
per metric listed in `plot_cols`, grouped by the column named in
`group_by`, with further faceting by the metric listed in `facet_by`.
The reliability column that corresponds to the metric to be plotted
should be provided in `reliability_col` to ensure that aggregates that
were excluded from the relevant analysis are not visualised. The user
can choose via `summary_metric` whether to report the mean or median and
whether this should be for all aggregates plotted or per group.

``` r

plot_nhood_beeswarm(metadata(sce_obj)$merged_nhood_analysis,
                    plot_cols = c("clone_per_contig" , "inv.simpson_diversity"),
                    reliability_col = "diversity_reliable",
                    group_by = "Bcelltype",
                    summary_metric = "group_median",
                    colour_by = "Bcelltype",
                    size_by = "nhood_size")
#> Orientation inferred to be along y-axis; override with
#> `position_quasirandom(orientation = 'x')`
```

![](TEVRA_files/figure-html/plot_nhood_beeswarm-1.png)

### Scatter plot

Two metrics can be compared via `plot_nhood_scatter`. The plot can be
facetted by naming an annotation of interest in `facet_by`. To ensure
unreliable aggregates are not plotted, both of the relevant reliability
columns corresponding to both metrics requested are provided.

``` r

plot_nhood_scatter(metadata(sce_obj)$merged_nhood_analysis,
                   x_col = "inv.simpson_diversity",
                   y_col = "CDR3_converge_score_VDJ",
                   x_reliability_col = "diversity_reliable",
                   y_reliability_col = "VDJ_converge_score_reliable",
                   colour_by = "Bcelltype",
                   size_by = "nhood_size")
```

![](TEVRA_files/figure-html/plot_nhood_scatter-1.png)

## Additional details on inter-clonal CDR3 convergence

Pairwise distances are computed via *pwalign*’s `stringDist()` function
across all unique CDR3 amino acid sequences per neighbourhood exceeding
min_n_contigs and optionally min_clones. CDR3 sequences can differ via
insertions, deletions or substitutions. In all cases, insertions and
deletions are assigned a cost of 1. For Levenshtein distances, all
substitutions also incur a cost of 1, hence the maximum cost is the
length of the longest sequence (all amino acids substituted, insertion
of the extra amino acids necessary to account for sequence length
variation). Substitutions can also be weighted according to their
biological significance via blosum-weighted scores. The raw BLOSUM
substitution matrix is converted to a cost matrix via \\c(i,j) =
(s(i,i) + s(j,j))/2 - s(i,j)\\. This cost matrix is then normalised to
\[0,1\] by dividing each cost by the maximum off-diagonal cost. The
normalised cost matrix is then negated such that the most costly amino
acid substitutions score -1 and the least significant substitutions
score 0. Overall, the maximum distances computed using blosum-weighted
substitutions is also the length of the longest sequence, which would be
acheived if all amino acids substitutions were the most costly
substitution. If the user wishes to use an alternative method of
substitution weighting, they must provide their own substitution matrix
which matches the same logic (bounded \[-1,0\], more similar amino acids
have a higher score).

![Transformation of BLOSUM62 substitution score
matrix](figures/blosum.png)

Transformation of BLOSUM62 substitution score matrix

Pairwise distances are computed between all CDR3 amino acid sequences
per aggregate using the above method. However, for BCRs, some of these
unique sequences belong to the same clone. In this scenario, the CDR3
sequence with the lowest average pairwise distance to all other
sequences in the clone is selected. All other sequences are discarded
from the matrix, such that it contains only one sequence per clone.

The distance matrix is then summarised to obtain a singular convergence
score representing the similarity of CDR3 sequences of different clones
per aggregate. This can be a global convergence score, defined as the
mean or median pairwise distance across all representative CDR3
sequences, or a local score. Local convergence scores are computed as
the mean distance of each CDR3 sequence to its k nearest neighbours,
averaged across all representatives in the aggregate. This metric is
more sensitive to repertoire structure than a global average. For
example, a neighbourhood containing clones that respond to two distinct
epitopes would contain two clusters of highly similar CDR3s that are
dissimilar from one another. Here, global convergence would be high due
to the difference between these two antigens, whilst local convergence
would reflect sequence similarity within these two groups. By default, k
is set adaptively at n/4,where n is the number of representatives, with
a minimum k of 3; a fixed k can also be supplied. Local convergence is
not calculated when the number of representatives does not exceed k.

Finally, all distance-based metrics are converted to convergence scores
via \\1 - distance\\.
