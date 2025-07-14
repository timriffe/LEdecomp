#  a <- .001
#  b <- .07
#  x <- 0:100
#  mx1 <- a * exp(x * b)
#  mx2 <- a/2 * exp(x * b)
#  sex1 = "m"
#  sex2 = "m"
#  closeout = TRUE
lopez_ruzicka <- function(mx1,
                          mx2,
                          age = (1:length(mx1))-1,
                          nx = rep(1,legth(mx1)),
                          sex1 = 't',
                          sex2 = sex1,
                          closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1,
                  age = age,
                  nx = nx,
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2,
                  age = age,
                  nx = nx,
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  nx = nx,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  nx = nx,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1,
                   nx = nx)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2,
                   nx = nx)
  Tx1 <- rcumsum(Lx1)
  Tx2 <- rcumsum(Lx2)
  ex1 <- Tx1 / lx1
  ex2 <- Tx2 / lx2

  # from here on, everything uses lx and ex only;
  # we need to lag each of these, for easier code reading
  ex1_next <- shift(ex1, n = -1, fill = 0)
  ex2_next <- shift(ex2, n = -1, fill = 0)
  lx1_next <- shift(lx1, n = -1, fill = 0)
  lx2_next <- shift(lx2, n = -1, fill = 0)

  exclusively_effect <- (lx1/lx2)*(lx2*(ex2-ex1) - lx2_next*(ex2_next - ex1_next))

  interaction_effect <- (ex2_next - ex1_next)*(((lx1*lx2_next)/lx2) - lx1_next)

  decomp <- exclusively_effect + interaction_effect
  decomp

}

lopez_ruzicka_sym <- function(mx1,
                              mx2,
                              age = (1:length(mx1))-1,
                              nx = rep(1,legth(mx1)),
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  a1 <- lopez_ruzicka(mx1,
                      mx2,
                      age = age,
                      nx = nx,
                      sex1 = sex1,
                      sex2 = sex2,
                      closeout = closeout)
  a2 <- lopez_ruzicka(mx2,
                      mx1,
                      age = age,
                      nx = nx,
                      sex1 = sex2,
                      sex2 = sex1,
                      closeout = closeout)

  a_avg <- (a1 - a2) / 2
  a_avg

}
sen_lopez_ruzicka <- function(mx1,
                              mx2,
                              age = (1:length(mx1))-1,
                              nx = rep(1,legth(mx1)),
                              sex1 = 't',
                              sex2 = sex1,
                              closeout = TRUE){
  # TR: replaces redundant code from before
  decomp <- lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex1,
                          sex2 = sex2,
                          closeout = closeout)
  delta <- mx2 - mx1
  sen <- decomp/delta
  sen

}
sen_lopez_ruzicka_sym <- function(mx1, mx2,
                                  age = (1:length(mx1))-1,
                                  nx = rep(1,legth(mx1)),
                                  sex1 = 't',
                                  sex2 = sex1,
                                  closeout = TRUE){
  delta <- mx2 - mx1
  a_avg <- lopez_ruzicka_sym(mx1 = mx1,
                             mx2 = mx2,
                             age = age,
                             nx = nx,
                             sex1 = sex1,
                             sex2 = sex2,
                             closeout = closeout)
  a_avg / delta
}

sen_lopez_ruzicka_instantaneous <- function(mx,
                                            age = (1:length(mx1))-1,
                                            nx = rep(1,legth(mx1)),
                                            sex = 't',
                                            perturb = 1e-6,
                                            closeout = TRUE){
  mx1 <- mx * (1 / (1 - perturb))
  mx2 <- mx * (1 - perturb) / 1
  s1 <- sen_lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  s2 <- sen_lopez_ruzicka(mx1 = mx2,
                          mx2 = mx1,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  # TR: is this closeout actually needed here?
  # I have my doubts this was checked
  #To match the solution with the arriaga
  # if (closeout){
  #   s2[length(s2)] <- s2[length(s2)] * 2
  # }
  sen <- (s1 + s2) / 2

  sen
}
sen_lopez_ruzicka_instantaneous2 <- function(mx,
                                             age = (1:length(mx1))-1,
                                             nx = rep(1,legth(mx1)),
                                             sex = 't',
                                             perturb = 1e-6,
                                             closeout = TRUE){
  mx1 <- exp(log(mx) + perturb)
  mx2 <- exp(log(mx) - perturb)
  s1 <- sen_lopez_ruzicka(mx1 = mx1,
                          mx2 = mx2,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  s2 <- sen_lopez_ruzicka(mx1 = mx2,
                          mx2 = mx1,
                          age = age,
                          nx = nx,
                          sex1 = sex,
                          sex2 = sex,
                          closeout = closeout)
  # TR: is this closeout actually needed here?
  # I have my doubts this was checked
  #To match the solution with the arriaga
  # if (closeout){
  #   s2[length(s2)] <- s2[length(s2)] * 2
  # }
  sen <- (s1 + s2) / 2

  sen
}
