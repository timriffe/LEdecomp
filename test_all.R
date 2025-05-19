a <- .001
b <- .07
x <- 0:100
mx1 <- a * exp(x * b)
mx2 <- a/2 * exp(x * b)

mx1_c <- matrix(NA, 101, 3, dimnames = list(x, c("c1", "c2", "c3")))
mx2_c <- matrix(NA, 101, 3, dimnames = list(x, c("c1", "c2", "c3")))

mx1_c[,1] <- mx1/2
mx1_c[,2] <- mx1/3
mx1_c[,3] <- mx1/6

mx2_c[,1] <- mx2/4
mx2_c[,2] <- mx2/2
mx2_c[,3] <- mx2/4

rowSums(mx1_c) - mx1
rowSums(mx2_c) - mx2

#devtools::load_all()
library(DemoTools)
library(dplyr)
library(data.table)

#Discrete Life Expectancy -- ALL CAUSES
mx_to_e0(mx2,age=x,sex="t") - mx_to_e0(mx1,age=x,sex="t")  # 9.501863

methods <- c(
"lifetable", "arriaga", "arriaga_sym",
"sen_arriaga", "sen_arriaga_sym",
"sen_arriaga_inst", "sen_arriaga_inst2",
"chandrasekaran_ii","chandrasekaran_ii_sym",
"sen_chandrasekaran_ii", "sen_chandrasekaran_ii_inst",
"sen_chandrasekaran_ii_inst2", "sen_chandrasekaran_ii_sym",
"chandrasekaran_iii","chandrasekaran_iii_sym",
"sen_chandrasekaran_iii", "sen_chandrasekaran_iii_inst",
"sen_chandrasekaran_iii_inst2", "sen_chandrasekaran_iii_sym",
"lopez_ruzicka", "lopez_ruzicka_sym",
"sen_lopez_ruzicka", "sen_lopez_ruzicka_sym",
"sen_lopez_ruzicka_inst", "sen_lopez_ruzicka_inst2",
"horiuchi", "stepwise", "numerical")

# default settings:
D_defaults <- rep(0,length(methods))
      names(D_defaults) <-   methods
for (m in methods){
  D_defaults[[m]] <- LEdecomp(mx1,mx2,age=x,sex1="t",sex2="t",method=m,opt=FALSE)
}
