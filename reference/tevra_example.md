# Example dataset: murine splenic B cells

A subset of the 10x Genomics murine splenocyte dataset with paired gene
expression and BCR sequencing, used in the TEVRA tutorial.

## Usage

``` r
tevra_example
```

## Format

A `SingleCellExperiment` with 2,000 highly variable genes and 3,430
cells:

- assays:

  `logcounts`: log-normalised expression

- reducedDims:

  `PCA`, `UMAP`

- colData:

  Per-cell metadata with 19 columns:

  `orig.ident`

  :   Tissue of origin (spleen)

  `nCount_RNA`

  :   Total UMI counts per cell

  `nFeature_RNA`

  :   Number of genes detected per cell

  `percent.mt`

  :   Percentage of counts from mitochondrial genes

  `celltype`

  :   Broad cell type annotation (B cell)

  `S.Score`, `G2M.Score`

  :   Cell cycle phase scores from Seurat's
      [`CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html)

  `Phase`

  :   Assigned cell cycle phase (G1, S or G2M)

  `Bcelltype`

  :   B cell subset annotation

  `CTgene_BCR`

  :   V(D)J and constant gene usage of the heavy and light chains

  `CTnt_BCR`

  :   CDR3 nucleotide sequences of the heavy and light chains

  `CTaa_BCR`

  :   CDR3 amino acid sequences of the heavy and light chains

  `CTstrict_BCR`

  :   Strict clonotype identifier from scRepertoire, used as the default
      clone definition

  `clonalProportion_BCR`

  :   Proportion of cells in the sample belonging to the cell's clone

  `clonalFrequency_BCR`

  :   Number of cells in the sample belonging to the cell's clone

  `cloneSize_BCR`

  :   Binned clonal expansion category

  `CTaa_BCR_IGH`

  :   Heavy chain CDR3 amino acid sequence

  `CTaa_BCR_IGL`

  :   Light chain CDR3 amino acid sequence

- metadata:

  List of additional objects:

  `graphs`

  :   Cell-cell graphs computed on the full gene set with Seurat's
      [`FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html):

      `RNA_snn`

      :   Shared nearest neighbour graph

## Source

Dataset from 10x Genomics
(<https://www.10xgenomics.com/datasets/integrated-gex-and-vdj-analysis-of-connect-generated-library-from-mouse-splenocytes-2-standard-6-0-1>),
licensed under CC BY 4.0
(<https://creativecommons.org/licenses/by/4.0/>). Modified from the
original: cell types were assigned with Seurat, clonal identity was
computed with scRepertoire, and the object was converted to a
`SingleCellExperiment` and restricted to the top 2,000 highly variable
genes.
