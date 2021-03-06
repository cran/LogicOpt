% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l.partybans.1-data.R
\docType{data}
\name{l.partybans.1}
\alias{l.partybans.1}
\title{Logicopt truth table created from "partybans.csv" dataset}
\format{R data frame table}
\source{
compass.org website
}
\usage{
data(l.partybans.1)
}
\description{
l.partybans.1 is an logicopt compatible truth table generated 
from the QCA dataset "partybans.csv" where output "PB" is 1.
}
\examples{
\dontrun{
# Read raw QCA dataset from csv file
inpath <- system.file("extdata/raw_qca/partybans.csv", package="LogicOpt")
partybans <- read.csv(inpath,row.names=1,na="")

# Load QCA package 
library(QCAGUI) 

# Create the QCA truth table 
q.partybans.1 <- truthTable(partybans, conditions = c("C","F","T","R","V"), outcome = "PB{1}")

# Create the logicopt truth table
l.partybans.1 <- QCAtt2LOtt(q.partybans.1)
}

# Load up logicopt truth table
data(l.partybans.1)

# Optimize logicopt truth table and print results 
partybans1 <- logicopt(l.partybans.1,5,1,find_dc=TRUE,mode="multi-min")
print_multi_tt(partybans1,eqn=TRUE,n_in=5,n_out=1,QCA=TRUE)
}
\keyword{QCA}
\keyword{datasets}

