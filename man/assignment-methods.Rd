\name{assignment-methods}
\docType{methods}
\alias{assignment}
\alias{assignment-methods}
\alias{assignment,LINCbio-method}
\alias{assignment,LINCcluster-method}
\alias{assignment,LINCmatrix-method}
\alias{assignment,LINCsingle-method}
\title{Methods for the Getter Function \code{assignment} in Package \pkg{LINC}}
\description{
Access the slot \code{assignment} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
assignment slot
}

\item{\code{signature(x = "LINCcluster")}}{
assignment slot
}

\item{\code{signature(x = "LINCmatrix")}}{
assignment slot
}

\item{\code{signature(x = "LINCsingle")}}{
assignment slot
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

assignment(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~assignment }
