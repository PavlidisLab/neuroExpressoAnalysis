#' The list of studies included in the neuroexpresso database.
#' 
#' The list of studies included in the neuroexpresso database. Includes studies that are not included in the 
#' gene selection and the neuroexpresso website for the sake of completion. These studies can be filtered
#' out using ShinyNames or PyramidalDeep fields
#' 
#' @format
#' \itemize{
#'  \item GSE: GSE identifier of the dataset
#'  \item samples: GSM identifiers or file names of the samples included
#'  \item MajorType: Does the study represent glial or neuronal cells
#'  \item Neurotransmitter: Grouping of studies based on neurotransmitter release
#'  \item ShinyNames: Cell type names used in the neuroexpresso.org
#'  \item PyramidalDeep: internal shorter cell type names for cell types used in the study
#'  \item Description: Description provided by the original authors
#'  \item Age: Average age of the samples. Adult denotes an uncertain adult mouse while negative
#'   numbers denote embryonic days
#'  \item Region: Brain regions that the samples were included in for the study
#'  \item Platform: Microarray platform that was used for the samples
#'  \item Reference: Source of the dataset
#'  \item PMID: PubMed id of the paper
#'  \item Normalize: Samples that are used when merging GPL339 and GPL1261 platform
#'  \item Normalize2: Samples that are used when using only GPL1261 platform
#'  \item Note: Notes about the quality of the sample or about the changes made on the naming or metadata
#'  }
"n_expressoStudies"



'n_expressoStudies'

'n_expressoSamples'

