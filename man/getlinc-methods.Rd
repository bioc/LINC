\name{getlinc-methods}
\docType{methods}
\alias{getlinc-methods}
\alias{getlinc,ANY,character-method}
\alias{getlinc}


\title{
 Subsetting for \code{LINC} objects
}

\description{
\code{getlinc} is a function to derive substructures from \code{LINC} objects.
}
\section{Methods}{
\describe{

\item{\code{signature(input = "ANY", subject = "character")}}{
(see details)
}
}}

\usage{
getlinc(input,
        subject = "queries")
}

\arguments{
  \item{input}{
 a \code{'LINC'} object, for instance \code{LINCmatrix}
}
  \item{subject}{
has to be one of \code{c("queries", "geneSystem", "results", "history", "customID")}
}
}


\value{
the respective substructure or information
  }

\author{
Manuel Goepferich
}

\seealso{
\code{linc} ; 
}
\examples{
data(BRAIN_EXPR)

# getlinc() is used to accesss information
getlinc(crbl_cluster, subject = "geneSystem")
getlinc(crbl_cluster, subject = "queries")
getlinc(crbl_cluster, subject = "customID")

}  
\keyword{ ~getlinc }

\keyword{methods}
