load_all()
?plot.LEdecomp
data("US_data_CoD", package = "LEdecomp")
allc <- subset(US_data_CoD, Period == 2010 & cause == "All-causes") |>
  as.data.frame()

# Make Female vs Male all-cause schedules, Age 0:100
ac_w <- reshape(allc[, c("Gender","Age","mxt")],
                timevar = "Gender", idvar = "Age", direction = "wide")
names(ac_w) <- sub("^mxt\\.", "", names(ac_w))
ac_w <- ac_w[order(ac_w$Age), ]

dec_ac <- LEdecomp(
  mx1 = ac_w$Male,
  mx2 = ac_w$Female,
  age = 0:100,
  method = "sen_arriaga"
)

# Simple single-line plot

plot(dec_ac, main = "All-cause Arriaga, 2010 Female vs Male",what="sens")

## End(Not run)
## Example 2: Cause of death, one year, Female vs Male
cod <- subset(US_data_CoD, Period == 2010 & cause != "All-causes")
cod_w <- reshape(cod[, c("Gender","Age","cause","mxt")],
                 timevar = "Gender", idvar = c("cause","Age"),
                 direction = "wide")|>
  as.data.frame()
names(cod_w) <- sub("^mxt\\.", "", names(cod_w))
cod_w <- cod_w[order(cod_w$cause, cod_w$Age), ]

dec_cod <- LEdecomp(
  mx1 = cod_w$Male,
  mx2 = cod_w$Female,
  age = 0:100,
  n_causes = length(unique(cod_w$cause)),
  method = "sen_arriaga",
  cause_names = unique(cod$cause_id)
)

# Overlay of all causes

plot(dec_cod, layout = "overlay", main = "Arriaga CoD, 2010 Female vs Male", legend.pos = "top")

# Facet by cause (3 columns)
plot(dec_cod, layout = "facet", ncol = 3, main = "Arriaga by cause (faceted)")
plot(dec_cod, layout = "facet", ncol = 3, main = "Arriaga by cause (faceted)",what="sens")
age_ab <- c(0L, 1L, seq.int(5L, 110L, by = 5L))
n_age  <- length(age_ab)
# nx inferred from age grid (repeat last width for open age group)
nx_ab  <- c(diff(age_ab), tail(diff(age_ab), 1L))

# Synthetic Gompertz rates on abridged lower bounds
a <- 0.0001
b <- 0.07
mx1_all <- a * exp(age_ab * b)
mx2_all <- (a / 2) * exp(age_ab * b)

# Make 3 causes by random weights per age, normalized to 1
k <- 3
w1 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
w2 <- matrix(runif(n_age * k, 0.9, 1.1), nrow = n_age)
w1 <- w1 / rowSums(w1)
w2 <- w2 / rowSums(w2)

mx1_mat <- mx1_all * w1
mx2_mat <- mx2_all * w2

# Matrix baseline (explicit age and nx)
res_mat <- LEdecomp(
  mx1 = mx1_mat, mx2 = mx2_mat,
  age = age_ab, nx = nx_ab,
  sex1 = "t", sex2 = "t",
  method = "sen_arriaga", opt = TRUE,facet=TRUE
)
plot(res_mat, what = "sens")
res_mat

library(devtools)library(devtools)cause_id
library(usethis)
library(rhub)
document()

devtools::check()
check_win_devel()      #
check_win_release()    #
check_win_oldrelease() #
rhub_platforms()
rhub::rhub_check(platforms = c("linux","macos","macos-arm64","windows"))

library(spelling)
spell_check()

# doesn't work until on CRAN?
#library(revdepcheck)
#revdepcheck::revdep_check(cran = FALSE)

devtools::release()
