\name{correlation-methods}
\docType{methods}
\alias{correlation-methods}
\alias{correlation,LINCbio-method}
\alias{correlation,LINCcluster-method}
\alias{correlation,LINCmatrix-method}
\alias{correlation,LINCsingle-method}
\title{Methods for the Getter Function \code{correlation} in Package \pkg{LINC}}
\description{
Access the slot \code{correlation} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
correlation slot
}

\item{\code{signature(x = "LINCcluster")}}{
correlation slot
}

\item{\code{signature(x = "LINCmatrix")}}{
correlation slot
}

\item{\code{signature(x = "LINCsingle")}}{
correlation slot
}
}}

\arguments{
  \item{x}{
 a \code{'LINC'} object, for instance \code{LINCmatrix}
}

\value{
the respective substructure or information
}

\author{
Manuel Goepferich
}

\examples{
data(BRAIN_EXPR)

correlation(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~correlation }
