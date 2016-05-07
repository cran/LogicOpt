\name{LogicOpt-package}

\alias{LogicOpt-package}

\alias{LogicOpt}

\docType{package}

\title{
LogicOpt: A Package for Logic Optimization of Truth Tables using Espresso 
}

\description{
The package \pkg{LogicOpt} provides access to powerful logic minimization
algorithms that operate on a sum-of-products truth table. The core algorithms
are built on Espresso Version 2.3 developed at UC Berkeley for digital logic
synthesis purposes.  The espresso C code has been extended extensively 
to integrate within the R framework and support additional logic optimization 
use cases that operate on R data frame tables.  The primary interface to the 
package is through the \code{\link{logicopt}} function.  

New algorithms have been developed that leverage the espresso core routines
and data structures to support logic optimization use cases used by Qualitative 
Comparative Analysis (QCA).  See the "primes", "multi-min", and "multi-full" 
options for the "mode" paramater in function \code{\link{logicopt}}. These modes 
along with the paramater find_dc=TRUE provide features similar to those
provided by \code{QCAGUI::eqmcc} and \code{QCApro::eQMC} but the espresso-based
routines here are able to handle much larger functions and run in less 
runtime.  Integration to QCA is through the \code{QCAGUI::truthTable} or  
\code{QCApro::truthTable} functions provided in QCA packages \pkg{QCAGUI} 
and \pkg{QCApro} respectively (other QCA packages may also work but have not
been tested).  These functions provide the necessary functionality
to convert a raw QCA dataset into a QCA truth table (see documentation in those
packages for more details).  The function \code{\link{QCAtt2LOtt}} then provides 
the last step to get to an expresso compatible truth table.  See dataset 
\code{\link{l.represent.1}} for an example flow from a raw QCA dataset to 
final optimized QCA results.

This package also contains a Genetic Programming use case where a maze 
navigation program for a robot is optimized.  See the paper referenced below
as well as the example truth table \code{\link{l.robot1}} for more details. The 
\code{\link{logicopt}} option mode="echo" was developed to support the abiilty 
to view the input truth table prior to optimization.  

}

\details{
    \tabular{ll}{
    Package: \tab LogicOpt\cr
    Type:    \tab Package\cr
    Version: \tab 0.1.0\cr
    Date:    \tab 2016-04-15\cr
    License: \tab GPL (>= 2)\cr
    }
}

\author{

\bold{Authors}:\cr
William Stiehl\cr
\email{wwstiehl@gmail.com}

\bold{Maintainer}:\cr
William Stiehl
}

\references{

Brayton, Robert King; Hachtel, Gary D.; McMullen, Curtis T.; Sangiovanni-Vincentelli, Alberto L.. (1984), Logic Minimization Algorithms for VLSI Synthesis, Kluwer Academic Publishers, ISBN 0-89838-164-9

Rudell, Richard L. (1986-06-05), "Multiple-Valued Logic Minimization for PLA Synthesis" Memorandum No. UCB/ERL M86-65 (Berkeley) \url{http://www.eecs.berkeley.edu/Pubs/TechRpts/1986/ERL-86-65.pdf}

Dusa, Adrian (2016). QCAGUI: Modern functions for Qualitative Comparative
Analysis. R Package URL: \url{http://cran.r-project.org/package=QCAGUI}

Thiem, Alrik. 2016. Professional Functionality for Performing and
Evaluating Qualitative Comparative Analysis. R Package Version 1.1-0.
\url{http://www.alrik-thiem.net/software}

Keane, A.J. 2015. "Genetic Programming, Logic Design and Case-Based Reasoning for
Obstacle Avoidance."  Learning and Intelligent Optimization: 9th International 
Conference; pp.104-118

}

