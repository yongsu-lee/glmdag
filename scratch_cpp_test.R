sourceCpp("./expm_test.cpp")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

set.seed(1)
x = rnorm(25)
X = matrix(x, 5, 5)
X
expm_test(X)
