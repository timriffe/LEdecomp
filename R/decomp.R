
#' LEdecomp
#' @title decompose life expectancy using a variety of methods
#' @description A variety of exact or asympototically exact life expectancy decomposition methods are implemented. These include the lifetable response experiment, using three potential sensitivity functions (lifetable-based `ltre_lt`, an instantaneous Arriaga-derived sensitivity `ltre_arriaga_instantaneous`, or numerical sensitivity `ltre_numerical`), the directional Arriaga method `arriaga`, a symmetrical Arriaga method `arriaga_symmetrical`. These all give similar results, except for directional `arriaga`, which is the most different.
#' @param mx1 numeric. age-structured mortality rates for population 1
#' @param mx2 numeric. age-structured mortality rates for population 2, must be same length as `mx1`
#' @param age integer. lower bound of each age group, must be same length as `mx1`
#' @param sex1 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param sex2 character. `"m"`,`"f"`, or `"t"`, affects a0 treatment.
#' @param method character. `"lifetable"`,`"arriaga_instantaneous"`,`"arriaga_instantaneous2"`,`"numerical"`,`"arriaga"`, `"arriaga_symmetrical"`
#' @param opt logical, default `TRUE`. For sensitivity-based decompositions, shall we optimize the rate averaging to eliminate the decomposition residual?
#' @param tol numeric, default `1e-10`, tolerance parameter for rate averaging optimization.
#' @param closeout logical. Do we handle closeout, or truncate at top age_
#' @param ... optional arguments passed to `numDeriv::grad()`
#' @export
# needs examples

 # a <- .001
 # b <- .07
 # x <- 0:100
 # mx1 <- a * exp(x * b)
 # mx2 <- a/2 * exp(x * b)
 # 
 # R <- seq(.1,.7,length=101)
 # mx1 <- cbind(mx1 * R, mx1 * (1-R))
 # mx2 <- cbind(mx2 * R, mx2 * (1-R))
 #LEdecomp(mx1,mx2,age=x,sex1='t',method = "arriaga_instantaneous")
 # 
 LEdecomp <- function(mx1, 
                      mx2, 
                      age, 
                      sex1 = 't', 
                      sex2 = sex1, 
                      method = c("lifetable","arriaga_instantaneous","arriaga_instantaneous2","numerical","arriaga", "arriaga_symmetrical"), 
                      closeout = TRUE,
                      opt = TRUE,
                      tol = 1e-10,
                      extras = FALSE,
                       ...){
  stopifnot(is.vector(mx1) | is.matrix(mx1))
  stopifnot(length(mx1) == length(mx2))
  
  
  # handle dimensions
  deez_dims <- dim(mx1)
  nages <- length(age)
  nmx <- length(mx1)
  if (nages < nmx){
    ncauses <- nmx / nages
    ncauses <- round(ncauses)
    stopifnot((ncauses * nages) == nmx)
    dim(mx1) <- c(nages,ncauses)
    dim(mx2) <- c(nages,ncauses)  
  }
  
  # add optimization option.
  
   method <- tolower(method)
   method <- match.arg( method,
                       choices = c("lifetable","arriaga_instantaneous","arriaga_instantaneous2","numerical","arriaga", "arriaga_symmetrical"))
   ccDONE <- FALSE
   # if sex1 != sex2, we re-call the function twice and take the average...
   if ((sex1 != sex2) & !method %in% c("arriaga","arriaga_symmetrical")){
     .sex1 <- sex1
     .sex2 <- sex2
     cc1 <- decomp_e0(mx1 = mx1, 
                      mx2 = mx2, 
                      age = age, 
                      sex1 = .sex1, 
                      sex2 = .sex1, 
                      method = method,
                      closeout = closeout,
                      N=N,
                      ...)
     cc2 <- decomp_e0(mx1 = mx1, 
                      mx2 = mx2, 
                      age = age, 
                      sex1 = .sex2, 
                      sex2 = .sex2, 
                      method = method,
                      closeout = closeout,
                      N=N,
                      ...)
     cc <- (cc1 + cc2) / 2
     dim(cc) <- deez_dims
     return(cc)
   }
   
   if (length(dim(mx1)) == 2){
     mxc1 <- mx1
     mx1  <- rowSums(mx1)
     
     mxc2 <- mx2
     mx2  <- rowSums(mx2)
     
     delta_xc <- mxc2 - mxc1
     causes   <- TRUE
   } else {
     causes <- FALSE
   }
   delta_x <- mx2 - mx1
   
    if (method == "arriaga"){
     cc <- arriaga(mx1 = mx1,
                   mx2 = mx2,
                   sex1 = sex1,
                   sex2 = sex2,
                   age = age,
                   closeout = closeout)
     if (causes){
       warning("Careful doing a cause-of-death decomposition with the classic Arriaga method, as results can be unstable. We recommend using a sensitivity-based approach for cause of death decompositions.")
       cc <- cc * (delta_xc / delta_x)
     }
     dim(cc) <- deez_dims
     return(cc)
   }
   if (method == "arriaga_symmetrical"){
     cc <- arriaga_sym(mx1 = mx1,
                   mx2 = mx2,
                   sex1 = sex1,
                   sex2 = sex2,
                   age = age,
                   closeout = closeout)
     if (causes){
       warning("Careful doing a cause-of-death decomposition with the classic Arriaga method, as results can be unstable. We recommend using a sensitivity-based approach for cause of death decompositions.")
       cc <- cc * (delta_xc / delta_x)
     }
     dim(cc) <- deez_dims
     return(cc)
   }
   
   # In this case we have a sensitivity-based decomposition, yay.
   # two blocks: either we optimize or we don't
   
   sen_fun <- switch(method,
                     arriaga_instantaneous = sen_arriaga_instantaneous,
                     arriaga_instantaneous2 = sen_arriaga_instantaneous2,
                     lifetable = sen_e0_mx_lt,
                     numerical = sen_num
   )
   if (!opt){
     mx <- (mx1 + mx2) / 2
     s <- sen_e0_mx(mx, age = age, sen_fun = sen_fun, sex = sex1)
  
   } else {
     # otherwise, we optimize the rate averaging
     s <- sen_min(mx1 = mx1,
                  mx2 = mx2,
                  age = age,
                  sex1 = sex1,
                  sex2 = sex2,
                  closeout = closeout,
                  sen_fun = sen_fun, 
                  tol = tol,
                  ...)
     }
   if (causes){
     cc <- s * delta_xc  
   } else {
     cc <- s * delta_x   
   }
   dim(cc) <- deez_dims
   return(cc)
 }

# a wrapper for sensitivity functions:
sen_e0_mx <- function(mx, age = 0:(length(mx)-1), 
                      sen_fun = sen_arriaga_instantaneous,
                      sex = 't',
                      closeout = TRUE){
    s = sen_fun(mx, age = age, sex = sex, closeout = closeout)
}












