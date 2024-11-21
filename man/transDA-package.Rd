\name{transDA-package}
\alias{transDA-package}
\alias{transDA}
\docType{package}
\title{
Transformation Discriminant Analysis
}
\description{
Transformation and mixture discriminant analysis.
}
\details{
\tabular{ll}{
Package: \tab transDA\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2024-6-10\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}
Function \code{'tda'} is the main function of the package that allows running traditional linear and quadratic discriminant procedures as well as more flexible options such as mixture discriminant analysis, transformation discriminant analysis, or their combination.

Function \code{'predict.tda'} makes classification predictions for a provided testing data set.

}
\author{
Jing Li and Yana Melnykov

Maintainer: Jing Li \email{jli178@crimson.ua.edu}
}
\references{
Hastie, T., & Tibshirani, R. (1996). \emph{Discriminant analysis by Gaussian mixtures}. Journal of the Royal Statistical Society Series B: Statistical Methodology, 58(1), 155-176.
  
Yana Melnykov & Jing Li. (2024)  \emph{Transformation Discriminant Analysis.} (under review)

}
% \references{Not ready yet}
\keyword{Transformation}
\keyword{mixture}
\keyword{discriminant}
\seealso{
  \code{\link{tda}}, 
  \code{\link{predict.tda}}
%%  ~~ Optional links to other man pages, e.g. ~~
%%  ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
set.seed(123)
data(iris)
MDA <- tda(x = iris[,1:4], max_k = 2, ID = iris$Species, trans = FALSE)
summary(MDA)
}
