  
 #  a <- .001
 #  b <- .07
 #  x <- 0:100
 # mx1 <- a * exp(x * b)
 #  mx2 <- a/2 * exp(x * b)
 #  sex1 = "m"
 #  sex2 = "m"
 #  closeout = TRUE
  chandrasekaran <- function(mx1,
                           mx2, 
                           age, 
                           sex1, 
                           sex2, 
                           closeout = TRUE){
  ax1 <- mx_to_ax(mx = mx1, 
                  age = age, 
                  sex = sex1,
                  closeout = closeout)
  ax2 <- mx_to_ax(mx = mx2, 
                  age = age, 
                  sex = sex2,
                  closeout = closeout)
  qx1 <- mx_to_qx(mx = mx1,
                  ax = ax1,
                  closeout = closeout)
  qx2 <- mx_to_qx(mx = mx2,
                  ax = ax2,
                  closeout = closeout)
  lx1 <- qx_to_lx(qx1)
  lx2 <- qx_to_lx(qx2)
  dx1 <- lx_to_dx(lx1)
  dx2 <- lx_to_dx(lx2)
  Lx1 <- ald_to_Lx(ax = ax1,
                   lx = lx1,
                   dx = dx1)
  Lx2 <- ald_to_Lx(ax = ax2,
                   lx = lx2,
                   dx = dx2)
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
  
  # eq 1.1
  main_effect <- (lx1 / lx2) * 
                 (lx2 * (ex2 - ex1) - 
                    lx2_next * (ex2_next - ex1_next))
  
  # eq 1.2
  operative_effect <- (lx2 / lx1) * 
    (lx1 * (ex2 - ex1) - lx1_next * (ex2_next - ex1_next))
  
  # eq 1.3
  interaction_effect_deferred <- 
    lx2 * (ex2 - ex1) - lx2_next * (ex2_next - ex1_next)
  
  # eq 1.4
  interaction_effect_forwarded <-
    lx1 * (ex2 - ex1) - lx1_next * (ex2_next - ex1_next)
  
  # eq 1.5 (open age group already handled, since 
  # *_next values close with 0s)
  approachII <- ((ex2 - ex1) * (lx2 + lx1)) / 2 - ((ex2_next - ex1_next) * (lx2_next+lx1_next)) / 2

  # eq 1.7 (avg of operative and main effects)
  exclusive_effect <- ((ex2 - ex1) * (lx2 + lx1)) / 2 -
    ((ex2_next - ex1_next)*((lx1 * lx2_next)/lx2 + (lx2*lx1_next)/lx1)) / 2
  
  
  # eq 1.8 (toal interaction effect = difference between total and exclusive)
  approachIII <- (ex2_next - ex1_next) * ((((lx1 * lx2_next) / lx2 + (lx2 * lx1_next) / lx1)-(lx2_next + lx1_next)) / 2)
  
  # 1.9 and 1.10 are decompositions of 1.8 per Arriaga
  # eq 1.9
  direct_effect <- (lx2 + lx1) * ((ex2 - ex1) + (((lx1_next * ex1_next) / lx1) - ((lx2_next * ex2_next) / lx2 ))) * 1/2
    
  # eq 1.10
  indirect_effect <-  (ex2_next * (lx2_next - (lx2 * lx1_next) / lx1) -
    ex1_next * ( lx1_next - (lx1 * lx2_next) / lx2)) / 2
  
  sum(direct_effect)
  
  }
  
  
