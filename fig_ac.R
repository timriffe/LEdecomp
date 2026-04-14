# custom all-cause decomposition figure,
# the markup and some finer parts of plot creation were
# created with assistance of ChatGPT 5.2

library(devtools)
library(LEdecomp)
library(tidyverse)
library(colorspace)
library(ggrepel)

dec_data<-
  US_data |>
  ungroup() |>
  filter(year %% 5 == 0,
         year >= 2010) |>
  pivot_wider(names_from = sex, values_from = mx)

all_decs <-
  dec_data |>
  group_by(year) |>
  mutate(contribution =
           LEdecomp(mx1 = Male,
                    mx2 = Female,
                    age = age,
                    sex1 = "m",
                    sex2 = "f",
                    method = "sen_arriaga_sym_inst") |>
           pluck("LEdecomp"))

all_decs |>
  group_by(year) |>
  summarize(Delta = sum(contribution))
dec2010 <- dec_data |> filter(year == 2010)
dec2015 <- dec_data |> filter(year == 2015)
dec2020 <- dec_data |> filter(year == 2020)

dec2010 <- LEdecomp(mx1 = dec2010$Male,
                    mx2 = dec2010$Female,
                    age = 0:100,
                    sex1 = "m",
                    sex2 = "f",
                    method = "arriaga_sym")
dec2015 <- LEdecomp(mx1 = dec2015$Male,
                    mx2 = dec2015$Female,
                    age = 0:100,
                    sex1 = "m",
                    sex2 = "f",
                    method = "arriaga_sym")
dec2020 <- LEdecomp(mx1 = dec2020$Male,
                    mx2 = dec2020$Female,
                    age = 0:100,
                    sex1 = "m",
                    sex2 = "f",
                    method = "arriaga_sym")

g1 <- dec2010$LE2-dec2010$LE1
g2 <- dec2015$LE2-dec2015$LE1
g3 <- dec2020$LE2-dec2020$LE1
gaps<- c(g1,g2,g3) |> round(digits = 2)
all_decs |>
  ggplot(aes(x=age,y=contribution,color = as.factor(year))) +
  geom_line()+
  theme_minimal()


# begin plot setup
# labels at left edge (min age) per year
labs_left <- all_decs |>
  group_by(year)|>
  filter(age == min(age, na.rm = TRUE)) |>
  slice_tail(n = 1) |>
  ungroup() |>
  mutate(gap = gaps,
         label = paste0(year," (",gap,"-year gap)"))

# points at age 0
pts0 <- all_decs  |>
  group_by(year)|>
  filter(age == min(age, na.rm = TRUE)) |>
  ungroup()

# unit box coordinates
x1 <- 50; x2 <- 75
y1 <- 0;  y2 <- 0.04

# corner length as fractions of box width/height
dx <- 3.5        # years on x-axis
dy <- 0.004      # years on y-axis

corner_col <- "#6F879B"
clw = .45

# plot construction (one big chain)

p <- all_decs %>%
  ggplot(aes(x = age, y = contribution, color = year, group = year)) +
  geom_line(linewidth = 0.8) +

  # add points for age 0 since lines largely overlap
  geom_point(
    data = pts0,
    aes(x = age, y = contribution),
    size = 2.2,
    shape = 21,
    fill = "white",
    stroke = 0.9
  ) +

  # custom direct year - labels on the left
  geom_text_repel(
    data = labs_left,
    aes(label = label),
    direction = "y",
    nudge_x = 2,      # move to the right
    hjust = 0,        # left-justify
    seed = 1,
    box.padding = 0.25,
    point.padding = 0.2,
    segment.size = 0.25,
    size = 4
  ) +

  scale_color_binned_sequential(palette = "YlOrRd",
                                guide = "none", begin = .3, end = 1) +

  labs(x = "Age",
       y = "Contribution to sex gap (female - male), years") +

  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    panel.grid.major.x = element_line(linewidth = 0.35),
    panel.grid.major.y = element_line(linewidth = 0.35),
    panel.grid.minor.y = element_line(linewidth = 0.2),
    panel.grid.minor.x = element_line(linewidth = 0.25)
  ) +
  # use x and y grids that make it easy to identify areas that sum to
  # a year, or a half-year.
  scale_x_continuous(
    breaks = seq(0, 100, 25),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.12, 0.05))) +
  scale_y_continuous(
    breaks = seq(0, 0.12, 0.04),
    minor_breaks = seq(0, 0.12, 0.02)) +

  # identify clean 1-year box using corners
  # bottom-left
  annotate("segment", x = x1, y = y1, xend = x1 + dx, yend = y1,
           linewidth = clw, color = corner_col) +
  annotate("segment", x = x1, y = y1, xend = x1, yend = y1 + dy,
           linewidth = clw, color = corner_col) +

  # top-right
  annotate("segment", x = x2, y = y2, xend = x2 - dx, yend = y2,
           linewidth = clw, color = corner_col) +
  annotate("segment", x = x2, y = y2, xend = x2, yend = y2 - dy,
           linewidth = clw, color = corner_col) +

  # 1-year box label
  annotate("text",
           x = (x1 + x2) / 2,
           y = (y1 + y2) / 2,
           label = "1 year of gap",
           color = corner_col,
           size = 4.2,
           fontface = "plain")
p
ggsave("dec_ac.pdf",p,width=15,height=15,units="cm")
