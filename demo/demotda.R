
library(transDA)

# function that test MDA
MDA <- tda(x = iris[,1:4], ID = iris$Species, trans = FALSE, max_k = 2)
print(MDA)
summary(MDA)

# function that test LDA
LDA <- tda(x = iris[,1:4], ID = iris$Species, max_k = 1, trans = FALSE, common_sigma = TRUE)
print(LDA)
summary(LDA)
# function that test QDA
QDA <-  tda(x = iris[,1:4], ID = iris$Species, subgroup = c(1,1,1), trans = FALSE, common_sigma = FALSE)
print(QDA)
summary(QDA)
# function that test TQDA
TQDA <- tda(x = iris[,1:4], ID = iris$Species, subgroup = c(1,1,1), trans = TRUE, common_sigma = FALSE)
print(TQDA)
summary(TQDA)