#' Example dataset: murine splenic B cells
#'
#' A subset of the 10x Genomics murine splenocyte dataset with paired
#' gene expression and BCR sequencing, used in the TEVRA tutorial.
#'
#' @format A \code{SingleCellExperiment} with 2,000 highly variable genes
#'   and 3,430 cells:
#' \describe{
#'   \item{assays}{\code{logcounts}: log-normalised expression}
#'   \item{reducedDims}{\code{PCA}, \code{UMAP}}
#'   \item{colData}{Per-cell metadata with 19 columns:
#'     \describe{
#'       \item{\code{orig.ident}}{Tissue of origin (spleen)}
#'       \item{\code{nCount_RNA}}{Total UMI counts per cell}
#'       \item{\code{nFeature_RNA}}{Number of genes detected per cell}
#'       \item{\code{percent.mt}}{Percentage of counts from mitochondrial genes}
#'       \item{\code{celltype}}{Broad cell type annotation (B cell)}
#'       \item{\code{S.Score}, \code{G2M.Score}}{Cell cycle phase scores from
#'         Seurat's \code{CellCycleScoring()}}
#'       \item{\code{Phase}}{Assigned cell cycle phase (G1, S or G2M)}
#'       \item{\code{Bcelltype}}{B cell subset annotation}
#'       \item{\code{CTgene_BCR}}{V(D)J and constant gene usage of the heavy and
#'         light chains}
#'       \item{\code{CTnt_BCR}}{CDR3 nucleotide sequences of the heavy and
#'         light chains}
#'       \item{\code{CTaa_BCR}}{CDR3 amino acid sequences of the heavy and
#'         light chains}
#'       \item{\code{CTstrict_BCR}}{Strict clonotype identifier from
#'         scRepertoire, used as the default clone definition}
#'       \item{\code{clonalProportion_BCR}}{Proportion of cells in the sample
#'         belonging to the cell's clone}
#'       \item{\code{clonalFrequency_BCR}}{Number of cells in the sample
#'         belonging to the cell's clone}
#'       \item{\code{cloneSize_BCR}}{Binned clonal expansion category}
#'       \item{\code{CTaa_BCR_IGH}}{Heavy chain CDR3 amino acid sequence}
#'       \item{\code{CTaa_BCR_IGL}}{Light chain CDR3 amino acid sequence}
#'     }
#'   }
#'   \item{metadata}{List of additional objects:
#'     \describe{
#'       \item{\code{graphs}}{Cell-cell graphs computed on the full gene set
#'         with Seurat's \code{FindNeighbors()}:
#'         \describe{
#'           \item{\code{RNA_snn}}{Shared nearest neighbour graph}
#'         }
#'       }
#'     }
#'   }
#' }
#' @source Dataset from 10x Genomics
#'   (\url{https://www.10xgenomics.com/datasets/integrated-gex-and-vdj-analysis-of-connect-generated-library-from-mouse-splenocytes-2-standard-6-0-1}),
#'   licensed under CC BY 4.0
#'   (\url{https://creativecommons.org/licenses/by/4.0/}). Modified from the
#'   original: cell types were assigned with Seurat, clonal identity was
#'   computed with scRepertoire, and the object was converted to a
#'   \code{SingleCellExperiment} and restricted to the top 2,000 highly
#'   variable genes.
"tevra_example"
