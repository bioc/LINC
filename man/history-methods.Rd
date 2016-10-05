\name{history-methods}
\docType{methods}
\alias{history}
\alias{history-methods}
\alias{history,LINCbio-method}
\alias{history,LINCcluster-method}
\alias{history,LINCmatrix-method}
\alias{history,LINCsingle-method}
\title{Methods for the Getter Function \code{history} in Package \pkg{LINC}}
\description{
Access the slot \code{history} of an object of class \code{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(x = "LINCbio")}}{
history slot
}

\item{\code{signature(x = "LINCcluster")}}{
history slot
}

\item{\code{signature(x = "LINCmatrix")}}{
history slot
}

\item{\code{signature(x = "LINCsingle")}}{
history slot
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

history(crbl_cluster)

} 

\keyword{methods}
\keyword{ ~history }
