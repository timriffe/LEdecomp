# library(DemoDecomp)
#  a <- 0.001
#  b <- 0.07
#  age <- 0:100
#  mx <- a * exp(age * b)
#  sex = 't'
#  perturb = 1e-5
# sen_horiuchi <- function(mx, age = 0:(length(mx)-1), sex = 't', perturb = 1e-5){
#   mx1 <- mx * (1 / (1 - perturb))
#   mx2 <- mx * ((1 - perturb) / 1)
#   
#   cc <- horiuchi(mx_to_e0, mx1, mx2, age = age, sex = sex, N = 20)
#   sen <- cc / (mx2 - mx1)
#   sen
# }

# plot(age, sen_horiuchi(mx, age) - sen_arriaga_instantaneous(mx))
# plot(age, sen_horiuchi(mx, age) - numDeriv::grad(mx_to_e0,mx,age = age, sex = sex))
# plot(age[-1], sen_horiuchi(mx, age)[-1] - sen_e0_mx_lt(mx,age,sex)[-1])

a <-.001
b <-0.07
x<-0:100
mx1 <- a* exp(x*b)
mx2 <- a/2 * exp(x*b)
mx_to_e0 <- function(mx){
  sum(exp(cumsum(-mx)))
}
mx_to_e0(mx2) - mx_to_e0(mx1)
cc1 <- DemoDecomp::horiuchi(mx_to_e0,mx1,mx2,N=100)
cc2 <- DemoDecomp::horiuchi(mx_to_e0,mx2,mx1,N=100)
sum(cc1)
sum(cc2)
all(abs(cc1 + cc2) < 1e-12)
