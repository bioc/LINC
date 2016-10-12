\name{getbio-methods}
\docType{methods}
\alias{getbio-methods}
\alias{getbio}
\alias{getbio,LINCcluster-method}
\title{
 Find Enriched Biological Terms For A Cluster
}
\description{
The function provides an interface to the \code{clusterProfiler} package. For each query in a cluster it seeks the biological terms that can be associated with the co-expressed genes, respectively. The input for \code{getbio} is of the class \code{'LINCcluster'}. 
}
\section{Methods}{
\describe{

\item{\code{signature(cluster = "LINCcluster")}}{
(see details)
}
}}

\usage{
getbio(cluster,
       enrichFun   = 'enrichGO',
       ont = "BP",
       ...)
}

\arguments{
  \item{cluster}{
 a \code{'LINCcluster'}. The number of co-expressed genes has to be sufficient.
}
  \item{enrichFun}{
a function given as character string which will derive significant biological terms based on the set of co-expressed genes from a gene annotation resource. Supported functions are:  \code{c("enrichGO", "enrichPathway", "enrichDO")}  
}
  \item{ont}{
a subontology, only used for \code{enrichFun = 'enrichGO'}. This has to be one of \code{"MF", "BP", "CC"}.
}
  \item{...}{
further arguments, mainly for functions from \code{clusterProfiler}
}
}

\details{
In contrast to the function \code{singlelinc} here, a group of queries, those present in the input cluster, will be analyzed for enriched biological terms. The annotation function can be one of \code{c("enrichGO", "enrichPathway", "enrichDO")}. [1]
The gene system of the input object has to be translated for the enrichment function in case genes are not given as Entrez ids. The function \code{clusterProfiler:bitr} [2] will be used in order to translate gene ids.
}

\value{
an object of the class 'LINCmatrix' (S4) with 6 \code{Slots}
  \item{results}{a \code{list} containing the identified enriched biological terms plus their respective p-values}
    \item{assignment}{a \code{character} vector of protein-coding genes}
  \item{correlation}{ a \code{list} of \code{cormatrix}, the correlation of non-coding to protein-coding genes and 
  \code{lnctolnc}, the correlation of non-coding to non-coding genes}
  \item{expression}{the original expression matrix}
  \item{history}{a storage environment of important methods, objects and parameters used to create the object}
  \item{linCenvir}{a storage environment ensuring the compatibility to other objects of the \code{LINC} class}
}

\references{
[1] Yu G, Wang L, Han Y and He Q (2012). "clusterProfiler: an R package for comparing biological themes among gene clusters." OMICS: A Journal of Integrative Biology, 16(5), pp. 284-287. 
(https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
}
\author{
Manuel Goepferich
}

\seealso{
\code{clusterlinc} ; 
\code{singlelinc} ;
}
\examples{
data(BRAIN_EXPR)

## Find the enriched cellular components for each query in the cluster 
crbl_cc <- getbio(crbl_cluster, ont = "CC")
 plotlinc(crbl_cc)
}  
\keyword{ ~getbio }
\keyword{methods}
