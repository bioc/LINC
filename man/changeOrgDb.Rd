\name{changeOrgDb}
\alias{changeOrgDb}
\title{
Change the Gene Annotation / Model Organism
}
\description{
The standard gene annotation in \pkg{LINC} is \code{"org.Hs.eg.db"}. This function is only relevant in case the input gene expression matrix is from another organism than Homo Sapiens.
}
\usage{
changeOrgDb(object, OrgDb)
}
\arguments{
  \item{object}{
an instance of the \code{LINC} class
}
  \item{OrgDb}{
has to be one of:
                    \code{c("org.Ag.eg.db",
                            "org.At.tair.db",
                            "org.Bt.eg.db",
                            "org.Cf.eg.db",
                            "org.Ce.eg.db",
                            "org.Gg.eg.db",
                            "org.Pt.eg.db",
                            "org.Sco.eg.db", 
                            "org.EcK12.eg.db",
                            "org.EcSakai.eg.db",
                            "org.Dm.eg.db",
                            "org.Tgondii.eg.db",
                            "org.Hs.eg.db",
                            "org.Pf.plasmo.db",
                            "org.Mm.eg.db",
                            "org.Ss.eg.db",
                            "org.Rn.eg.db",
                            "org.Mmu.eg.db",
                            "org.Xl.eg.db",
                            "org.Sc.sgd.db",
                            "org.Dr.eg.db",
                            "anopheles",   
                            "arabidopsis", 
                            "bovine",      
                            "canine",      
                            "celegans",    
                            "chicken",    
                            "chimp",       
                            "coelicolor",  
                            "ecolik12",    
                            "ecsakai",     
                            "fly",         
                            "gondii",      
                            "human",       
                            "malaria",     
                            "mouse",       
                            "pig",         
                            "rat",         
                            "rhesus",      
                            "xenopus",    
                            "yeast",       
                            "zebrafish")}
}
}

\value{
object with changed gene annotation hook
}

\author{
Manuel Goepferich
}
\note{
Please, also consider the documentation of \pkg{clusterProfiler} and \pkg{ReactomePA}.
}

\examples{
data(BRAIN_EXPR)

# change the used gene annotation, here from "human" to "mouse"
murine_matrix <- changeOrgDb(crbl_matrix, OrgDb = 'org.Mm.eg.db')
}  
\keyword{ ~changeOrgDb }