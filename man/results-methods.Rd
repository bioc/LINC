\name{results-methods}
\docType{methods}
\alias{results-methods}
\alias{results}
\alias{results,LINCbio-method}
\alias{results,LINCcluster-method}
\alias{results,LINCmatrix-method}
\alias{results,LINCsingle-method}
\title{Methods for the Getter Function \code{results} in Package \pkg{LINC}}
\description{
Access the slot \code{results} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
results slot
}

\item{\code{signature(x = "LINCcluster")}}{
results slot
}

\item{\code{signature(x = "LINCmatrix")}}{
results slot
}

\item{\code{signature(x = "LINCsingle")}}{
results slot
}
}}

\arguments{
  \item{x}{
 a \code{'LINC'} object, for instance \code{LINCmatrix}
}
}

\value{
the respective substructure or information
}

\author{
Manuel Goepferich
}

\examples{
data(BRAIN_EXPR)

results(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~results }
