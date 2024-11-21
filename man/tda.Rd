\name{tda}
\alias{tda}
\alias{print.tda}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Transformation Discriminant Analysis}

\description{
Implements discriminant analysis methods including traditional linear (LDA), quadratic (QDA), transformation (TDA), mixture (MDA) discriminant analysis, and their combinations such as TQDA or TLMDA. The user chooses a specific method by specifying options for common or varying transformation parameters as well as covariance matrices.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tda(x, max_k, ID, trans = TRUE, common_lambda = FALSE,
                 common_sigma = FALSE, iter = 50, subgroup = NULL, 
                 tol= 0.001, lambda0 = 0.015)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A frame or matrix containing a training data set}
  \item{max_k}{The maximum number of mixture components within each class to be fitted}
  \item{ID}{A variable containing class memberships for all observations}
  \item{trans}{A transformation indicator: \code{'trans = TRUE'} if transformation is needed, \code{'trans = FALSE'} if transformation is not needed}
  \item{common_lambda}{A parameter that regulates transformations. If \code{'common_lambda = TRUE'}, each mixture component or class has the same transformation parameter vector. If \code{'common_lambda = FALSE'}, each component or class has a different transformation vector}
  %\item{prop}{REMOVE THIS OPTION as we discussed}
  \item{common_sigma}{A homoscedasticity parameter: if \code{'common_sigma = TRUE'}, all subgroups across all classes have a common covariance matrix, if \code{'common_sigma = FALSE'}, groups have varying covariance matrices}
  \item{iter}{A maximum number of iterations of the EM algorithm; the default value is 50}
  \item{subgroup}{A vector containing the number of mixture components per each class to be fitted}
  \item{tol}{Tolerance level for a stopping critetion based on the relative difference in two consecutive log-likelihood values}
  \item{lambda0}{Starting value for transformation parameters}
}
%% \details{%%  ~~ If necessary, more details than the description above ~~}
\value{
 \item{BIC}{Values of the Bayesian Information Criterion calculated for each evaluated model}
\item{subprior }{Estimated component priors for each class}
\item{mu }{Estimated component means for each class}
\item{sigma }{Estimated component covariance matrices for each group}
\item{lambda }{Estimated transformation parameters}
%\item{True_ID }{REMOVE THIS AS REDUNDANT} 
\item{loglik }{The log-likelihood value for the model with the lowest BIC}
\item{pred_ID }{Estimated classification of observations in the training data set}
\item{prior }{Estimated class priors}
\item{misclassification_rate }{Misclassification rate for the training data set}
\item{ARI}{Adjusted Rand index value}
\item{Z}{Matrix of posterior probabilities for the training data set}
%\item{test_data_pred_ID}{predict ID on testing data.}
%\item{test_data_true_ID}{true ID of test dataset}
%\item{test_data_misclassification_rate}{misclassification rate of testing dataset}
%\item{test_data_ARI }{Adjust rand index on test data.}
}
%% \references{%% ~put references to the literature/web site here ~}
%% \author{%%  ~~who you are~~}
%% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{summary.tda}}, 
  \code{\link{predict.tda}}
}
\examples{

set.seed(123)
# Example 1:
MDA <- tda(x = iris[,1:4], max_k = 2,ID = iris$Species, trans = FALSE)
print(MDA)
summary(MDA)

# Example 2:
LDA <- tda(x = iris[,1:4], max_k = 1, ID = iris$Species, trans = FALSE,
        common_sigma = TRUE)
print(LDA)
summary(LDA)

# Example 3:
QDA <-  tda(x = iris[,1:4], subgroup = c(1, 1, 1), ID = iris$Species, 
        trans = FALSE, common_sigma = FALSE)
print(QDA)
summary(QDA)

# Example 4:
TQDA <- tda(x = iris[,1:4], subgroup = c(1, 1, 1), ID = iris$Species, 
        trans = TRUE, common_sigma = FALSE, common_lambda = TRUE)
print(TQDA)
summary(TQDA)
}
\keyword{transformation}
\keyword{mixture model}
\keyword{discriminant analysis}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
