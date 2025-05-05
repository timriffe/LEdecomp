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
dle <- sen_e0_mx_lt((mx2+mx1)/2,age=x,sex='t',closeout=TRUE)
dle_2 <- LEdecomp(mx1, mx2, age = x, sex1 = 't', closeout = TRUE, method = "lifetable")

dle_2$sens - dle

dle_2$LEdecomp - dle*(mx2 - mx1)
sum(dle*(mx2 - mx1)); sum(dle_2$LEdecomp)

#Discrete Life Expectancy -- BY CAUSES
dle_2cau <- LEdecomp(mx1_c, mx2_c, age = x,
                     sex1 = 't', closeout = TRUE, method = "lifetable")
sum(dle_2$LEdecomp)
sum(dle_2cau$LEdecomp)

#Arriaga Decomposition
cc1 <- arriaga(mx1, mx2, age = x)
cc1_2 <- LEdecomp(mx1, mx2, age = x, method = "arriaga")
cc1_2
plot.LEdecomp(cc1_2)
grepl(cc1_2$method,pattern="sen")
cc1 - cc1_2$LEdecomp
sum(cc1); sum(cc1_2$LEdecomp)

cc2 <- arriaga(mx2, mx1, age = x)
cc2_2 <- LEdecomp(mx2, mx1, age = x, method = "arriaga")

cc2 - cc2_2$LEdecomp
sum(cc2); sum(cc2_2$LEdecomp)

#Arriaga Decomposition BY CAUSES
cc1_cau2 <- LEdecomp(mx1 = mx1_c, mx2 = mx2_c, age = x, method = "arriaga")
cc1_cau2$LEdecomp
cc1_cau2
plot(cc1_cau2)

sum(cc1_2$LEdecomp); sum(cc1_cau2$LEdecomp)

cc2_cau2 <- LEdecomp(mx2_c, mx1_c, age = x, method = "arriaga")

sum(cc2_2$LEdecomp); sum(cc2_cau2$LEdecomp)

#Arriaga Symmetric Decomposition
cc_sym <- arriaga_sym(mx1, mx2, age = age)
cc_sym2 <- LEdecomp(mx1, mx2, age = x, method = "arriaga_sym")

cc_sym - cc_sym2$LEdecomp
sum(cc_sym); sum(cc_sym2$LEdecomp)

#Arriaga Symmetric Decomposition BY CAUSES
cc_sym_cau <- LEdecomp(mx1_causes, mx2_causes, age = x, method = "arriaga_sym")

sum(cc_sym2$LEdecomp); sum(cc_sym_cau$LEdecomp)


plot(x, cc1, type = "l", col = "blue")
lines(x, -cc2, type = "l", col = "red")
lines(x, cc_sym, type = "l", col = "green")

#Arriaga sensitivities
s_arri1 <- sen_arriaga(mx1, mx2, age = x)
s_arri1_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_arriaga")

s_arri1 - s_arri1_2$sens
s_arri1*(mx2 - mx1) - s_arri1_2$LEdecomp

sum(s_arri1*(mx2 - mx1)); sum(s_arri1_2$LEdecomp)

s_arri2 <- sen_arriaga(mx2, mx1, age = x)
s_arri2_2 <- LEdecomp(mx2, mx1, age = x, method = "sen_arriaga")

s_arri2 - s_arri2_2$sens
s_arri2*(mx1 - mx2) - s_arri2_2$LEdecomp

sum(s_arri2 * (mx1 - mx2)); sum(s_arri2_2$LEdecomp)

#Arriaga sensitivities BY CAUSES
s_arri_cau_2 <- LEdecomp(mx1_c, mx2_c, age = x, method = "sen_arriaga")

sum(s_arri1_2$LEdecomp); sum(s_arri_cau_2$LEdecomp)

s_arri2_cau_2 <- LEdecomp(mx2_c, mx1_c, age = x, method = "sen_arriaga")

sum(s_arri2_2$LEdecomp); sum(s_arri2_cau_2$LEdecomp)

#Arriaga-symmetrical sensitivities
s_arri_sym <- sen_arriaga_sym(mx1, mx2, age = x)
s_arri_sym_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_arriaga_sym")

s_arri_sym - s_arri_sym_2$sens
s_arri_sym*(mx2 - mx1) - s_arri_sym_2$LEdecomp

sum(s_arri_sym * (mx2 - mx1)); sum(s_arri_sym_2$LEdecomp)

#Arriaga-symmetrical sensitivities BY CAUSES
s_arri_sym_cau <- LEdecomp(mx1_c, mx2_c, age = x, method = "sen_arriaga_sym")

sum(s_arri_sym_2$LEdecomp); sum(s_arri_sym_cau$LEdecomp)

#Arriaga-Instantaneous sensitivities
mx <- (mx1 + mx2)/2
s_arri_inst <- sen_arriaga_instantaneous(mx, age = x, sex = 't', perturb = 1e-6,
                                         closeout = TRUE)
s_arri_inst_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_arriaga_inst")

s_arri_inst - s_arri_inst_2$sens
s_arri_inst*(mx2 - mx1) - s_arri_inst_2$LEdecomp

sum(s_arri_inst * (mx2 - mx1)); sum(s_arri_inst_2$LEdecomp)

#Arriaga-Instantaneous sensitivities BY CAUSES
s_arri_inst_cau <- LEdecomp(mx1_causes, mx2_causes, age = x, method = "sen_arriaga_inst")

sum(s_arri_inst_2$LEdecomp); sum(s_arri_inst_cau$LEdecomp)

#Arriaga-Instantaneous2 sensitivities
mx <- (mx1 + mx2)/2
s_arri_inst2 <- sen_arriaga_instantaneous2(mx, age = x, sex = 't', perturb = 1e-6,
                                         closeout = TRUE)
s_arri_inst2_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_arriaga_inst2")

s_arri_inst2 - s_arri_inst2_2$sens
s_arri_inst2*(mx2 - mx1) - s_arri_inst2_2$LEdecomp

sum(s_arri_inst2 * (mx2 - mx1)); sum(s_arri_inst2_2$LEdecomp)

#Arriaga-Instantaneous2 sensitivities BY CAUSES
s_arri_inst2_cau <- LEdecomp(mx1_c, mx2_c, age = x, method = "sen_arriaga_inst2")

sum(s_arri_inst2_2$LEdecomp); sum(s_arri_inst2_cau$LEdecomp)


plot(x, cc1, type = "l", col = "blue")
lines(x, -cc2, type = "l", col = "red")
lines(x, cc_sym, type = "l", col = "green")
lines(x, s_arri1 * (mx2 - mx1), lty = 2, col = "purple")
lines(x, s_arri2 * (mx2 - mx1), lty = 2, col = "purple")
lines(x, s_arri_inst * (mx2 - mx1), lty = 2, col = "yellow")

#Chandrasekaran Approach II is symetric
cc_chan1 <- chandrasekaran_II(mx1 = mx1, mx2 = mx2, age = x)
cc_chan1_1 <- LEdecomp(mx1, mx2, age = x,
                       method = "chandrasekaran_ii")

cc_chan1 - cc_chan1_1$LEdecomp

sum(cc_chan1) - sum(cc_chan1_1$LEdecomp)

cc_chan1_2 <- LEdecomp(mx2, mx1, age = x,
                       method = "chandrasekaran_ii")

cc_chan1_1$LEdecomp + cc_chan1_2$LEdecomp
sum(cc_chan1_1$LEdecomp); sum(cc_chan1_2$LEdecomp)

#Chandrasekaran Approach II BY CAUSES
cc_chan1_cau <- LEdecomp(mx1 = mx1_c, mx2 = mx2_c, age = x,
                       method = "chandrasekaran_ii")

sum(cc_chan1_1$LEdecomp); sum(cc_chan1_cau$LEdecomp)

#Chandrasekaran sensitivity Approach II is symetric
cc_chan_sen <- LEdecomp(mx1, mx2, age = x,
                        method = "sen_chandrasekaran_ii")
cc_chan_sen$sens
sum(cc_chan_sen$LEdecomp)

#Chandrasekaran sensitivity Approach II is symetric BY CAUSES
cc_chan_sen_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                        method = "sen_chandrasekaran_ii")

sum(cc_chan_sen$LEdecomp); sum(cc_chan_sen_cau$LEdecomp)



#Chandrasekaran sensitivity inst Approach II is symetric
cc_chan_sen_inst <- LEdecomp(mx1, mx2, age = x,
                        method = "sen_chandrasekaran_ii_inst")
cc_chan_sen_inst$sens
sum(cc_chan_sen$LEdecomp)

#Chandrasekaran sensitivity inst Approach II is symetric BY CAUSES
cc_chan_sen_inst_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                             method = "sen_chandrasekaran_ii_inst")
sum(cc_chan_sen_inst$LEdecomp); sum(cc_chan_sen_inst_cau$LEdecomp)


#Chandrasekaran sensitivity inst2 Approach II is symetric
sum(sen_chandrasekaran_II_instantaneous2((mx1+mx2)/2, age = x,)*(mx2 - mx1))

cc_chan_sen_inst2 <- LEdecomp(mx1, mx2, age = x,
                        method = "sen_chandrasekaran_ii_inst2")
cc_chan_sen_inst2$sens
sum(cc_chan_sen_inst2$LEdecomp)

#Chandrasekaran sensitivity inst2 Approach II is symetric BY CAUSES
cc_chan_sen_inst2_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                              method = "sen_chandrasekaran_ii_inst2")

sum(cc_chan_sen_inst2$LEdecomp); sum(cc_chan_sen_inst2_cau$LEdecomp)


#Approach III is symetric
cc_chand2 <- chandrasekaran_III(mx1 = mx1, mx2 = mx2, age = x)
cc_chand2_1 <- LEdecomp(mx1 = mx1, mx2 = mx2, age = x,
                        method = "chandrasekaran_iii")

cc_chand2 - cc_chand2_1$LEdecomp
sum(cc_chand2); sum(cc_chand2_1$LEdecomp)

cc_chand2_2 <- LEdecomp(mx1 = mx2, mx2 = mx1, age = x,
                        method = "chandrasekaran_iii")

cc_chand2_1$LEdecomp + cc_chand2_2$LEdecomp

sum(cc_chand2_1$LEdecomp); sum(cc_chand2_2$LEdecomp)

#Approach III is symetric BY CAUSES
cc_chand2_cau <- LEdecomp(mx1 = mx1_c, mx2 = mx2_c, age = x,
                        method = "chandrasekaran_iii")

sum(cc_chand2_1$LEdecomp); sum(cc_chand2_cau$LEdecomp)


#Chandrasekaran sensitivity Approach III is symetric
cc_chan2_sen <- LEdecomp(mx1, mx2, age = x,
                        method = "sen_chandrasekaran_iii")
cc_chan2_sen$sens
sum(cc_chan2_sen$LEdecomp)

#Chandrasekaran sensitivity Approach III is symetric BY CAUSES
cc_chan2_sen_cua <- LEdecomp(mx1_c, mx2_c, age = x,
                         method = "sen_chandrasekaran_iii")
sum(cc_chan2_sen$LEdecomp); sum(cc_chan2_sen_cua$LEdecomp)

#Chandrasekaran sensitivity inst Approach II is symetric
cc_chan2_sen_inst <- LEdecomp(mx1, mx2, age = x,
                             method = "sen_chandrasekaran_ii_inst")
cc_chan2_sen_inst$sens
sum(cc_chan2_sen$LEdecomp)

#Chandrasekaran sensitivity inst Approach II is symetric BY CAUSES
cc_chan2_sen_inst_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                              method = "sen_chandrasekaran_ii_inst")
sum(cc_chan2_sen_inst$LEdecomp); sum(cc_chan2_sen_inst_cau$LEdecomp)


#Chandrasekaran sensitivity inst2 Approach II is symetric
cc_chan2_sen_inst2 <- LEdecomp(mx1, mx2, age = x,
                              method = "sen_chandrasekaran_ii_inst2")
cc_chan2_sen_inst2$sens
sum(cc_chan2_sen_inst2$LEdecomp)
#Chandrasekaran sensitivity inst2 Approach II is symetric BY CAUSES
cc_chan2_sen_inst2_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                               method = "sen_chandrasekaran_ii_inst2")
sum(cc_chan2_sen_inst2$LEdecomp); sum(cc_chan2_sen_inst2_cau$LEdecomp)


plot(x, cc_chand2_1$LEdecomp, type = "l", col = "blue")
lines(x, cc_chan2_sen$LEdecomp, type = "l", col = "red")
lines(x, cc_chan2_sen_inst$LEdecomp, type = "l", col = "green")
lines(x, cc_chan2_sen_inst2$LEdecomp, type = "l", col = "purple")

#Lopez-Ruzicka
cc_lr1 <- lopez_ruzicka(mx1 = mx1, mx2 = mx2, age = x)
cc_lr1_2 <- LEdecomp(mx1 = mx1, mx2 = mx2, age = x,
                     method = "lopez_ruzicka")
cc_lr1 - cc_lr1_2$LEdecomp
sum(cc_lr1) - sum(cc_lr1_2$LEdecomp)

cc_lr2 <- lopez_ruzicka(mx1 = mx2, mx2 = mx1, age = x)
cc_lr2_2 <- LEdecomp(mx1 = mx2, mx2 = mx1, age = x,
                     method = "lopez_ruzicka")
cc_lr2 - cc_lr2_2$LEdecomp
sum(cc_lr2) - sum(cc_lr2_2$LEdecomp)

cc_lr1_2$LEdecomp - cc_lr2_2$LEdecomp

plot(x, cc_lr1_2$LEdecomp, type = "l", col = "blue")
lines(x, -cc_lr2_2$LEdecomp, col = "red")

sum(cc_lr1); sum(cc_lr2)

#Lopez-Ruzicka BY CAUSES
cc_lr1_cc <- LEdecomp(mx1 = mx1_c, mx2 = mx2_c, age = x,
                     method = "lopez_ruzicka")
cc_lr2_cc <- LEdecomp(mx1 = mx2_c, mx2 = mx1_c, age = x,
                     method = "lopez_ruzicka")

sum(cc_lr1_2$LEdecomp); sum(cc_lr1_cc$LEdecomp)
sum(cc_lr2_2$LEdecomp); sum(cc_lr2_cc$LEdecomp)

#The sum is the same but the obtain different age-values.
cc_lr_sym <- lopez_ruzicka_sym(mx1 = mx1, mx2 = mx2, age = x)
cc_lr2_sym <- LEdecomp(mx1 = mx1, mx2 = mx2, age = x,
                       method = "lopez_ruzicka_sym")
sum(cc_lr_sym); sum(cc_lr2_sym$LEdecomp)

#LOPEZ-RUZICKA SYM - BY CAUSES
cc_lr2_sym_cau <- LEdecomp(mx1 = mx1_c, mx2 = mx2_c, age = x,
                           method = "lopez_ruzicka_sym")

sum(cc_lr2_sym$LEdecomp); sum(cc_lr2_sym_cau$LEdecomp)


plot(x, cc_lr1, type = "l", col = "blue")
lines(x, -cc_lr2, type = "l", col = "red")
lines(x, cc_lr_sym, type = "l", col = "green")

#Lopez-Ruzicka sensitivities
s_lr1_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_lopez_ruzicka")

s_lr1rri1_2$sens
sum(s_arri1_2$LEdecomp)

s_lr2_2 <- LEdecomp(mx2, mx1, age = x, method = "sen_lopez_ruzicka")

s_lr2_2$sens
sum(s_lr2_2$LEdecomp)
s_arri1_2$LEdecomp - s_lr2_2$LEdecomp

#Lopez-Ruzicka sensitivities BY CAUSES
s_lr1_2_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                        method = "sen_lopez_ruzicka")
s_lr2_2_cau <- LEdecomp(mx2_c, mx1_c, age = x,
                        method = "sen_lopez_ruzicka")

sum(s_lr2_2$LEdecomp); sum(s_lr1_2$LEdecomp)
sum(s_lr2_2_cau$LEdecomp); sum(s_lr1_2_cau$LEdecomp)


#Lopez-Ruzicka-symmetrical sensitivities
s_lr_sym_2 <- LEdecomp(mx1, mx2, age = x, method = "sen_lopez_ruzicka_sym")

s_lr_sym_2$sens
sum(s_lr_sym_2$LEdecomp)

#Lopez-Ruzicka-symmetrical sensitivities
s_lr_sym_2_cau <- LEdecomp(mx1_c, mx2_c, age = x,
                           method = "sen_lopez_ruzicka_sym")

sum(s_lr_sym_2$LEdecomp); sum(s_lr_sym_2_cau$LEdecomp)

plot(x, s_lr1_2$LEdecomp, type="l")
lines(x, -s_lr2_2$LEdecomp, col ="red")
lines(x, s_lr_sym_2$LEdecomp, col = "blue")

#Lopez-Ruzicka-Instantaneous sensitivities
s_lr_inst <- LEdecomp(mx1, mx2, age = x, method = "sen_lopez_ruzicka_inst")

s_lr_inst$sens
sum(s_lr_inst$LEdecomp)

#Lopez-Ruzicka-Instantaneous sensitivities BY CAUSES
s_lr_inst_cau <- LEdecomp(mx1_c, mx2_c, age = x, method = "sen_lopez_ruzicka_inst")

sum(s_lr_inst$LEdecomp); sum(s_lr_inst_cau$LEdecomp)

#Lopez-Ruzicka-Instantaneous2 sensitivities
s_lr_inst2 <- LEdecomp(mx1, mx2, age = x, method = "sen_lopez_ruzicka_inst")

s_lr_inst2$sens
sum(s_lr_inst2$LEdecomp);

#Lopez-Ruzicka-Instantaneous2 sensitivities BY CAUSES
s_lr_inst2_cau <- LEdecomp(mx1_c, mx2_c, age = x, method = "sen_lopez_ruzicka_inst")

sum(s_lr_inst2_cau$LEdecomp); sum(s_lr_inst2$LEdecomp)


plot(x, s_lr1_2$LEdecomp, type="l")
lines(x, -s_lr2_2$LEdecomp, col ="red")
lines(x, s_lr_sym_2$LEdecomp, col = "blue")
lines(x, s_lr_inst$LEdecomp, col = "green")
lines(x, s_lr_inst2$LEdecomp, col = "purple")

#Horiuchi Decomposition
cc_1 <- DemoDecomp::horiuchi(mx_to_e0, mx1, mx2, age = age, N = 20)
cc_LE_1 <- LEdecomp(mx1, mx2, age = age, Num_Intervals = 20, func = mx_to_e0,
                    method = "horiuchi")

cc_LE_cau <- LEdecomp(mx1_c, mx2_c, age = age, Num_Intervals = 20, func = mx_to_e0,
                    method = "horiuchi")

cc_1 - cc_LE_1$LEdecomp
sum(cc_1); sum(cc_LE_1$LEdecomp); sum(cc_LE_cau$LEdecomp)

cc_LE_2 <- LEdecomp(mx2, mx1, age = age, Num_Intervals = 20, func = mx_to_e0,
                    method = "horiuchi")

cc_LE_2

cc_LE_2$LEdecomp + cc_LE_1$LEdecomp

plot(x, cc_LE_1$LEdecomp, type = "l", col = "blue")
lines(x, -cc_LE_2$LEdecomp, type = "l", col = "red")

#STEPWISE - DEMODECOMP
dims <- c(101, 1)
Mxc2e0abrvec <- function(mxc, age, dims, sex = "t"){
  dim(mxc) <- dims
  mx <- rowSums(mxc)
  mx_to_e0(mx, age = age, sex = sex)
}
s_1 <- DemoDecomp::stepwise_replacement(func = Mxc2e0abrvec,
                                        pars1 = mx1, pars2 = mx2, dims = dims,
                                        age = age,
                                        symmetrical = TRUE, direction = "both")
s_2 <- LEdecomp(mx1, mx2, age = age, symmetrical = TRUE, func = Mxc2e0abrvec,
                direction = "both", method = "stepwise", dims = dims)
s_cau <- LEdecomp(mx1_c, mx2_c, age = age, symmetrical = TRUE, func = Mxc2e0abrvec,
                direction = "both", method = "stepwise", dims = dims)

sum(s_1); sum(s_2$LEdecomp)

#STEPWISE - DEMODECOMP

dle2 <- numDeriv::grad(mx_to_e0, (mx2+mx1)/2, age=x, sex = 't')
sum(dle2*(mx2 - mx1))
