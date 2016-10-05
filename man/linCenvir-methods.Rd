\name{linCenvir-methods}
\docType{methods}
\alias{linCenvir-methods}
\alias{linCenvir,LINCbio-method}
\alias{linCenvir,LINCcluster-method}
\alias{linCenvir,LINCmatrix-method}
\alias{linCenvir,LINCsingle-method}
\title{Methods for the Getter Function \code{linCenvir} in Package \pkg{LINC}}
\description{
Access the slot \code{linCenvir} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
linCenvir slot
}

\item{\code{signature(x = "LINCcluster")}}{
linCenvir slot
}

\item{\code{signature(x = "LINCmatrix")}}{
linCenvir slot
}

\item{\code{signature(x = "LINCsingle")}}{
linCenvir slot
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

linCenvir(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~linCenvir }
