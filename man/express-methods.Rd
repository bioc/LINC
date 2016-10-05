\name{express-methods}
\docType{methods}
\alias{express-methods}
\alias{express,LINCbio-method}
\alias{express,LINCcluster-method}
\alias{express,LINCmatrix-method}
\alias{express,LINCsingle-method}
\title{Methods for the Getter Function \code{express} in Package \pkg{LINC}}
\description{
Access the slot \code{expression} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
expression slot
}

\item{\code{signature(x = "LINCcluster")}}{
expression slot
}

\item{\code{signature(x = "LINCmatrix")}}{
expression slot
}

\item{\code{signature(x = "LINCsingle")}}{
expression slot
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

express(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~express }
