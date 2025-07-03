
library(devtools)
library(usethis)

document()

check()

# 1. Get sensitivity for avg all-cause rates,
# e.g. using sen_arriaga_instantaneous(), this returns a vector.
# let's call this 'sen'

# 2. get COD rate matrix, with causes in columns and ages in rows for men and women, take the element-wise difference. let's call this 'delta'

# 3. then the decomposition becomes sen * delta, here exploiting the R artifact that a vector times a matrix is interpreted as a column vector being element-wise multiplied into the rows of the matrix, like so:

# sen <- runif(10)
# delta <- matrix(runif(80), ncol = 8)
# 
# cc <- sen * delta
# cc




