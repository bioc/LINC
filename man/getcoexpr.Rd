\name{getcoexpr}
\alias{getcoexpr}

\title{
 Get IDs For Co-Expressed Genes from The 'LINC' Class
}
\description{
\code{getcoexpr} provides access to co-expressed genes of a query in \code{'LINC'} objects.
}
\usage{
getcoexpr(input,
          query   = NULL,
          keyType = NULL)
}

\arguments{
  \item{input}{
a \code{'LINCcluster'} or \code{'LINCsingle'} object
}
  \item{query}{
for a \code{'LINCcluster'} the name of the query gene
}
   \item{query}{
The \code{'keyType'} (gene system) for the output.
}
}

\value{
a vector containing the co-expressed genes for a query
}

\author{
Manuel Goepferich
}


\examples{
data(BRAIN_EXPR)
# Get the co-expressed genes for the query gene '55384' alias MEG3
getcoexpr(crbl_cluster, query = "55384")

# The co-expressed genes can also be returned as gene symbols.
getcoexpr(crbl_cluster, query = "55384", keyType = 'SYMBOL')
}  
\keyword{ ~getcoexpr }
