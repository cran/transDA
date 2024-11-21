\name{summary.tda}
\alias{summary.tda}
\alias{print.summary.tda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing discriminant analysis results from function \code{'tda'}
}
\description{
Summary for the results of discriminant analysis obtained by function \code{'tda'}
}
\usage{
\method{summary}{tda}(object, \dots)
\method{print}{summary.tda}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{'tda'} produced by function \code{'tda'}}
  \item{x}{An object of class \code{'summary.tda'} obtained by function \code{'summary.tda'}}
  \item{\dots}{Additional arguments that are passed to or from other methods}
}
% \details{%%  ~~ If necessary, more details than the description above ~~}
\value{
Function \code{'summary.tda'} calculates and returns a list of summary statistics for the model estimated by function \code{'tda'}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{%% ~put references to the literature/web site here ~}
% \author{%%  ~~who you are~~}
%% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tda}}
}
\examples{
set.seed(123)
model <- tda(x = iris[,1:4], max_k = 2, ID = iris$Species, trans = FALSE)
summary(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
