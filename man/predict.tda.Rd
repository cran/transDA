\name{predict.tda}
\alias{predict.tda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classification of new observations based on the model obtained by function \code{'tda'}
}
\description{
Predicts class memberships for new observations based on the model returned by function \code{'tda'}.
}
\usage{
\method{predict}{tda}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{'tda'} produced by function \code{\link{tda}}
}
  \item{newdata}{
A data frame or matrix with data to be classified
}
  \item{\dots}{
Additional arguments that are passed to or from other methods
}
}
%% \details{%%  ~~ If necessary, more details than the description above ~~}
\value{
  \item{classification}{A factor containing predicted class labels for \code{newdata}}
  \item{Z}{A matrix of posterior probabilities for observations from \code{newdata}}
}
% \references{%% ~put references to the literature/web site here ~}
% \author{%%  ~~who you are~~}
% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{tda}}}
\examples{
     set.seed(123)
     split_indices <- sample(nrow(iris), size = round(0.8 * nrow(iris)), replace = FALSE)
     split_indices <- sort(split_indices)
     traindata <- iris[split_indices,]
     testdata <- iris[-split_indices,]
     MDA <- tda(x = traindata[,1:4], max_k = 2, ID = traindata$Species, trans = FALSE)
     pred <- predict(MDA, testdata[,-5]); pred
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
