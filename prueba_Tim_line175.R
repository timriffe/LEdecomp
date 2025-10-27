library(dbplyr)
library(roxygen2)
library(devtools)

devtools::load_all()
load(file = "US_mortalityData.RData")
load(file = "US_mortalityData_CoD.RData")

#Load the US mortality data from HMD for Male and Female from 2000
US_mortdata
US_mortdata$mxt[1:100] - US_mortdata$Dxt[1:100]/US_mortdata$Ext[1:100]
US_mortdata_cod
US_mortdata_cod$mxt[1:100] -
  US_mortdata_cod$Dxt[1:100]/US_mortdata_cod$Ext[1:100]

print.LEdecompData <- function(x) {
  cat("Mortality Data\n")
  cat(attributes(x)$label, "including", attributes(x)$series ,"\n")
  cat("Periods", c(min(x$Period),":", max(x$Period)),"\n")
  cat("Complete Life Table with ages from", c(min(x$Age)," to ", max(x$Age)), "\n")
}

US_data <- structure(list(Age = US_mortdata$Age,
                          Gender = US_mortdata$Gender,
                          Period = US_mortdata$Period,
                          Ext  = US_mortdata$Ext,
                          Dxt  = US_mortdata$Dxt,
                          mxt = US_mortdata$mxt))

US_data <- as.data.frame(US_data)
attr(US_data, "label") <- "US total population"
attr(US_data, "series") <- "males and females"

class(US_data) <- c("LEdecompData", class(US_data))

class(US_data)
inherits(US_data, "data.frame")

US_data_CoD

US_data_CoD <- structure(list(Age = US_mortdata_cod$Age,
                              Gender = US_mortdata_cod$Gender,
                              Period = US_mortdata_cod$Period,
                              Ext  = US_mortdata_cod$Ext,
                              Dxt  = US_mortdata_cod$Dxt,
                              mxt = US_mortdata_cod$mxt,
                              cause = US_mortdata_cod$cause,
                              cause_id = US_mortdata_cod$cause_id))

US_data_CoD <- as.data.frame(US_data_CoD)
attr(US_data_CoD, "label") <- "US cause-of-death data"
attr(US_data_CoD, "series") <- "males and females"

class(US_data_CoD) <- c("LEdecompData", class(US_data_CoD))

class(US_data_CoD)
inherits(US_data_CoD, "data.frame")

US_data_CoD


library(dplyr)

save(US_data, file = "data/US_data.rda")
save(US_data_CoD, file = "data/US_data_CoD.rda")
usethis::use_data(US_data, US_data_CoD, internal = FALSE,
                  overwrite = TRUE)

#include in the file NAMESPACE:
library(Rdpack)
library(usethis)
S3method(print,US_data)

US_male2015 <- US_data[US_data$Period == 2015 &
                         US_data$Gender == 'Male',]
US_female2015 <- US_data[US_data$Period == 2015 &
                           US_data$Gender == 'Female',]

US_male2010 <- US_data[US_data$Period == 2010 &
                         US_data$Gender == 'Male',]
US_female2010 <- US_data[US_data$Period == 2010 &
                           US_data$Gender == 'Female',]


US_male2020 <- US_data[US_data$Period == 2020 &
                         US_data$Gender == 'Male',]
US_female2020 <- US_data[US_data$Period == 2020 &
                           US_data$Gender == 'Female',]

ar1 <- arriaga(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100))
sum(ar1)
ar2 <- arriaga(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100))
sum(ar2)

arriaga1 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
         method = "arriaga")
arriaga2 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                     method = "arriaga")
arriaga_sym20 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                        method = "arriaga_sym")
arriaga_sym15 <- LEdecomp(mx1 = US_male2015$mxt, mx2 = US_female2015$mxt, age = c(0:100),
                          method = "arriaga_sym")
arriaga_sym10 <- LEdecomp(mx1 = US_male2010$mxt, mx2 = US_female2010$mxt, age = c(0:100),
                          method = "arriaga_sym")

arriaga_sym10

plot(arriaga_sym10) +
  ggplot2::scale_y_continuous(limits = c(-0.0001, 0.115)) +
  ggplot2::labs(title = "(a)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)  # centrado
  )
plot(arriaga_sym15) +
  ggplot2::scale_y_continuous(limits = c(-0.0001, 0.115)) +
  ggplot2::labs(title = "(b)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)  # centrado
  )
plot(arriaga_sym20) +
  ggplot2::scale_y_continuous(limits = c(-0.0001, 0.115)) +
  ggplot2::labs(title = "(c)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5)  # centrado
  )

plot(c(50:90), arriaga_sym10$LEdecomp[50:90], type="l", ylim = c(0, 0.115))
plot(c(50:90), arriaga_sym15$LEdecomp[50:90], type="l", ylim = c(0, 0.115))
plot(c(50:90), arriaga_sym20$LEdecomp[50:90], type="l", ylim = c(0, 0.115))

sum(arriaga_sym10$LEdecomp[21:41])
sum(arriaga_sym15$LEdecomp[21:41])
sum(arriaga_sym20$LEdecomp[21:41])

sum(arriaga_sym10$LEdecomp[60:80])/sum(arriaga_sym10$LEdecomp)
sum(arriaga_sym15$LEdecomp[60:80])/sum(arriaga_sym15$LEdecomp)
sum(arriaga_sym20$LEdecomp[60:80])/sum(arriaga_sym20$LEdecomp)

#CAUSE-OF-DEATH ANALYSIS
#First, we select the data for period 2010, 2015, and 2020 for males and females
US_male2010_cod <- US_data_CoD[US_data_CoD$Period == 2010 &
                             US_data_CoD$Gender == 'Male',]
US_female2010_cod <- US_data_CoD[US_data_CoD$Period == 2010 &
                           US_data_CoD$Gender == 'Female',]
#Second, we keep only the data from cause-of-death -> cause =! "All-causes"
US_male2010_cod <- US_male2010_cod[US_male2010_cod$cause
                                   != 'All-causes',]
US_female2010_cod <- US_female2010_cod[US_female2010_cod$cause
                                       != 'All-causes',]
#Repeat with 2015
US_male2015_cod <- US_data_CoD[US_data_CoD$Period == 2015 &
                         US_data_CoD$Gender == 'Male',]
US_female2015_cod <- US_data_CoD[US_data_CoD$Period == 2015 &
                           US_data_CoD$Gender == 'Female',]
US_male2015_cod <- US_male2015_cod[US_male2015_cod$cause
                                   != 'All-causes',]
US_female2015_cod <- US_female2015_cod[US_female2015_cod$cause
                                       != 'All-causes',]

#Finally, with 2020
US_male2020_cod <- US_data_CoD[US_data_CoD$Period == 2020 &
                                 US_data_CoD$Gender == 'Male',]
US_female2020_cod <- US_data_CoD[US_data_CoD$Period == 2020 &
                                   US_data_CoD$Gender == 'Female',]
US_male2020_cod <- US_male2020_cod[US_male2020_cod$cause
                                   != 'All-causes',]
US_female2020_cod <- US_female2020_cod[US_female2020_cod$cause
                                       != 'All-causes',]

US_cod_2010<- LEdecomp(mx1 = matrix(US_male2010_cod$mxt, 101, 18),
                          mx2 = matrix(US_female2010_cod$mxt, 101, 18),
                          age = c(0:100), method = "arriaga")
US_cod_2010
arriaga_sym10

#CHECK TIM

cause_id <- unique(US_male2010_cod$cause_id)

max(US_cod_2010)

plot(US_cod_2010) +
  ggplot2::scale_y_continuous(limits = c(-0.01, 0.115)) +
  ggplot2::labs(title = "(a)") +
  ggplot2::scale_fill_discrete(
    labels = c(cause_id[1:18])) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),  # centrado
    legend.position = "bottom",
    legend.box = "horizontal",
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))


arriaga_sym10_cod <- LEdecomp(mx1 =  matrix(US_male2010_cod$mxt, 101, 18, byrow = T),
                              mx2 = matrix(US_female2010_cod$mxt, 101, 18, byrow = T),
                              age = c(0:100),
                              method = "arriaga_sym")
arriaga_sym10_cod
colSums(arriaga_sym10_cod$LEdecomp)
colSums(arriaga_sym15_cod$LEdecomp)
colSums(arriaga_sym20_cod$LEdecomp)

arriaga_sym15_cod <- LEdecomp(mx1 =  matrix(US_male2015_cod$mxt, 101, 18, byrow = T),
                              mx2 = matrix(US_female2015_cod$mxt, 101, 18, byrow = T),
                              age = c(0:100),
                              method = "arriaga_sym")
arriaga_sym15_cod
arriaga_sym15

arriaga_sym20_cod <- LEdecomp(mx1 =  matrix(US_male2020_cod$mxt, 101, 18, byrow = T),
                              mx2 = matrix(US_female2020_cod$mxt, 101, 18, byrow = T),
                              age = c(0:100),
                              method = "arriaga_sym")
arriaga_sym20_cod
arriaga_sym20

colSums(arriaga_sym10_cod$LEdecomp)
colSums(arriaga_sym15_cod$LEdecomp)
colSums(arriaga_sym20_cod$LEdecomp)

#Age 0
round(arriaga_sym10_cod$LEdecomp[1,], 4)
round(arriaga_sym15_cod$LEdecomp[1,], 4)
round(arriaga_sym20_cod$LEdecomp[1,], 4)

arriaga_sym10_cod$LEdecomp

max(rowSums(arriaga_sym10_cod$LEdecomp), rowSums(arriaga_sym15_cod$LEdecomp),
    rowSums(arriaga_sym20_cod$LEdecomp))
min(rowSums(arriaga_sym10_cod$LEdecomp), rowSums(arriaga_sym15_cod$LEdecomp),
    rowSums(arriaga_sym20_cod$LEdecomp))


as10 <- plot(arriaga_sym10_cod) +
  ggplot2::scale_y_continuous(limits = c(-0.005, 0.12)) +
  ggplot2::labs(title = "(a)") +
  ggplot2::scale_fill_discrete(
    labels = c(cause_id[1:18])) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),  # centrado
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8)
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))

as15 <- plot(arriaga_sym15_cod) +
  ggplot2::scale_y_continuous(limits = c(-0.005, 0.12)) +
  ggplot2::labs(title = "(b)") +
  ggplot2::scale_fill_discrete(
    labels = c(cause_id[1:18])) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),  # centrado
    legend.position = "bottom",
    legend.box = "horizontal",
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))

as20 <- plot(arriaga_sym20_cod) +
  ggplot2::scale_y_continuous(limits = c(-0.005, 0.12)) +
  ggplot2::labs(title = "(c)") +
  ggplot2::scale_fill_discrete(
    labels = c(cause_id[1:18])) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),  # centrado
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))

final_plot <- as10 | as15 | as20
final_plot

install.packages("patchwork")
library(patchwork)

plot(c(0:100), arriaga1$LEdecomp, col = "red", type = "l")
lines(c(0:100), ar1)
lines(c(0:100), -arriaga2$LEdecomp, col = "blue")
lines(c(0:100), -ar2, col = "green")

chan21 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                     method = "chandrasekaran_ii")
chan22 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                   method = "chandrasekaran_ii")
plot(c(0:100), chan21$LEdecomp, type = "l")
lines(c(0:100), -chan22$LEdecomp, col = "red")
lines(c(0:100), arriaga1$LEdecomp, col = "green")
lines(c(0:100), -arriaga2$LEdecomp, col = "green")
lines(c(0:100), arriaga_sym$LEdecomp, col = "blue")

chan31 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                   method = "chandrasekaran_iii")
chan32 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                   method = "chandrasekaran_iii")
plot(c(0:100), chan31$LEdecomp, type = "l")
lines(c(0:100), chan21$LEdecomp, col = "red")
lines(c(0:100), -chan32$LEdecomp, col = "black")
lines(c(0:100), -chan22$LEdecomp, col = "red")
lines(c(0:100), arriaga1$LEdecomp, col = "green")
lines(c(0:100), -arriaga2$LEdecomp, col = "green")
lines(c(0:100), arriaga_sym$LEdecomp, col = "blue")

lr1 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                   method = "lopez_ruzicka")
lr2 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                   method = "lopez_ruzicka")
plot(c(0:100), lr1$LEdecomp, type = "l")
lines(c(0:100), -lr2$LEdecomp, col = "red")
lines(c(0:100), chan31$LEdecomp, col = "black")
lines(c(0:100), chan21$LEdecomp, col = "blue")
lines(c(0:100), arriaga1$LEdecomp, col = "green")
lines(c(0:100), -arriaga2$LEdecomp, col = "green")
lines(c(0:100), arriaga_sym$LEdecomp, col = "blue")

step1 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                method = "stepwise")
step2 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                method = "stepwise")

plot(c(0:100), step1$LEdecomp, type = "l")
lines(c(0:100), -step2$LEdecomp, col = "red")
lines(c(0:100), arriaga_sym$LEdecomp, col = "blue")

hor1 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                  method = "horiuchi")
hor2 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                  method = "horiuchi")

plot(c(0:100), hor1$LEdecomp, type = "l")
lines(c(0:100), -hor2$LEdecomp, col = "red")
lines(c(0:100), arriaga_sym$LEdecomp, col = "blue")

step1 <- LEdecomp(mx1 = US_male2020$mxt, mx2 = US_female2020$mxt, age = c(0:100),
                 method = "stepwise")
step2 <- LEdecomp(mx2 = US_male2020$mxt, mx1 = US_female2020$mxt, age = c(0:100),
                 method = "stepwise")
sum(step1$LEdecomp)
sum(step2$LEdecomp)

plot(c(0:100), step1$LEdecomp, type = "l")
lines(c(0:100), -step2$LEdecomp, col = "red")
