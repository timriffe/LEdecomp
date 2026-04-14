# custom cause-specific decomposition figures,
# the markup setup and some finer parts of plot creation were
# created with assistance of ChatGPT 5.2

library(tidyverse)
library(colorspace)
library(grid)
library(LEdecomp)

# start w cause by age plot, lots of data prep and prelims to making this plot
all_dec_cod <-
  US_data_CoD |>
  select(-cause,-cause_id) |>
  filter(year %in% c(2010,2015,2020)) |>
  pivot_wider(names_from = sex, values_from = mxc) |>
  arrange(year, cause_short, age) |>
  group_by(year) |>
  mutate(contribution =
           LEdecomp(mx1 = Male,
                    mx2 = Female,
                    n_causes = 18,
                    age = age,
                    sex1 = "m",
                    sex2 = "f",
                    method = "sen_arriaga_sym_inst",
                    opt = TRUE) |>
           pluck("LEdecomp"))

# 0) Inputs
all_dec_cod$cause_short |> unique()

# 1) Recode causes (COVID rule + Other/ill-defined)
all_dec_cod <- all_dec_cod |>
  mutate(
    year = factor(year),
    cause_plot = case_when(
      str_detect(cause_short, "^Special codes") & year == "2020" ~ "COVID-19",
      str_detect(cause_short, "^Special codes")                  ~ "Other / ill-defined",
      cause_short %in% c("Symptoms & ill-defined conditions (R00-R99)",
                         "Maternal (O00-O99)",
                         "Skin (L00-L98)")                        ~ "Other / ill-defined",
      TRUE                                                       ~ cause_short
    )
  )

# 2) total per cause-year
cause_year_tot <-
  all_dec_cod |>
  group_by(year, cause_plot) |>
  summarize(total = sum(contribution, na.rm = TRUE), .groups = "drop")

# mean abs total across years
cause_score <-
  cause_year_tot |>
  group_by(cause_plot) |>
  summarize(score = mean(abs(total), na.rm = TRUE),
            .groups = "drop")

thr <- 0.02

keep_always <- c(
  "Congenital (Q00-Q99)" # because visible in infants
)

small <- cause_score|>
  filter(
    score < thr,
    !(cause_plot %in% c("COVID-19", "Other / ill-defined")),
    !(cause_plot %in% keep_always)
  )  |>
  pull(cause_plot)

# 4) collapse
all_dec_cod <- all_dec_cod |>
  mutate(
    cause_plot = if_else(cause_plot %in% small,
                         "Other / ill-defined",
                         cause_plot)
  )

# 3) stacking order (bottom -> top)
stack_levels <- c(
  "Other / ill-defined",
  "Perinatal (P00-P96)",
  "Congenital (Q00-Q99)",
  "Infectious (A00-B99)",
  "Blood & immune (D50-D89)",
  "Endocrine & metabolic (E00-E88)",
  "Digestive (K00-K92)",
  "Respiratory (J00-J98)",
  "Musculoskeletal (M00-M99)",
  "Genitourinary (N00-N98)",
  "Mental (F01-F99)",
  "Nervous system (G00-G98)",
  "Neoplasms (C00-D48)",
  "Circulatory (I00-I99)",
  "External causes (V01-Y89)",
  "COVID-19"
)

all_dec_cod <- all_dec_cod |>
  mutate(
    cause_plot = factor(cause_plot, levels = stack_levels)
  )



## 4) Palette

levs <- levels(droplevels(all_dec_cod$cause_plot))
pal <- setNames(rep(NA_character_, length(levs)), levs)

# anchors
if ("Other / ill-defined" %in% levs) pal["Other / ill-defined"] <- "#C9CED6"
if ("External causes (V01-Y89)" %in% levs) pal["External causes (V01-Y89)"] <- "#E07A2F"
if ("Circulatory (I00-I99)" %in% levs)      pal["Circulatory (I00-I99)"]      <- "#C23B30"
if ("Neoplasms (C00-D48)" %in% levs)        pal["Neoplasms (C00-D48)"]        <- "#6A51A3"
if ("COVID-19" %in% levs)                   pal["COVID-19"]                   <- "#4D4D4D"

assign_seq <- function(names_in, palette, ...) {
  names_in <- intersect(names_in, levs)
  if (length(names_in) == 0) return(invisible(NULL))
  cols <- colorspace::sequential_hcl(length(names_in), palette = palette, ...)
  pal[names_in] <<- cols
}

# families
assign_seq(c("Nervous system (G00-G98)", "Mental (F01-F99)",
             "Genitourinary (N00-N98)", "Musculoskeletal (M00-M99)"),
           palette = "Blue-Yellow", c1 = 40, c2 = 65, l1 = 85, l2 = 60)

assign_seq(c("Respiratory (J00-J98)", "Digestive (K00-K92)", "Endocrine & metabolic (E00-E88)"),
           palette = "Teal", c1 = 35, c2 = 60, l1 = 85, l2 = 60)

assign_seq(c("Perinatal (P00-P96)", "Congenital (Q00-Q99)"),
           palette = "OrYel", c1 = 25, c2 = 55, l1 = 90, l2 = 70)

assign_seq(c("Infectious (A00-B99)", "Blood & immune (D50-D89)"),
           palette = "Greens 3", c1 = 25, c2 = 55, l1 = 90, l2 = 65)

# Use two clearly separated hues (muted but distinct)- the abov emakes these too close
if ("Mental (F01-F99)" %in% levs) {
  pal["Mental (F01-F99)"] <- "#2F5597"   # deep muted blue
}
if ("Endocrine & metabolic (E00-E88)" %in% levs) {
  pal["Endocrine & metabolic (E00-E88)"] <- "#2A7F62"  # deep muted green-teal
}

# fallback if anything left
missing <- names(pal)[is.na(pal)]
if (length(missing) > 0) {
  pal[missing] <- colorspace::qualitative_hcl(length(missing), palette = "Dark 3", c = 55, l = 70)
}

## Legend order
direct_or_callout <- c(
  "External causes (V01-Y89)",
  "Circulatory (I00-I99)",
  "Neoplasms (C00-D48)",
  "COVID-19",
  "Perinatal (P00-P96)",
  "Congenital (Q00-Q99)",
  "Respiratory (J00-J98)",
  "Digestive (K00-K92)",
  "Endocrine & metabolic (E00-E88)",
  "Infectious (A00-B99)"
)

legend_levels_full <- c(
  "External causes (V01-Y89)",
  "Circulatory (I00-I99)",
  "Neoplasms (C00-D48)",
  "COVID-19",
  "Nervous system (G00-G98)",
  "Mental (F01-F99)",
  "Genitourinary (N00-N98)",
  "Musculoskeletal (M00-M99)",
  "Respiratory (J00-J98)",
  "Digestive (K00-K92)",
  "Endocrine & metabolic (E00-E88)",
  "Infectious (A00-B99)",
  "Blood & immune (D50-D89)",
  "Perinatal (P00-P96)",
  "Congenital (Q00-Q99)",
  "Other / ill-defined"
)

legend_levels_full <- legend_levels_full[legend_levels_full %in% levs]

legend_levels <-
  legend_levels_full[!(legend_levels_full %in% direct_or_callout)]


## 6) Labels (2020-only) ----

# In-area labels
labs_text_2020 <- tribble(
  ~year,  ~label,        ~x,  ~y,    ~col,    ~fontface,
  "2020", "External",    32,  0.030, "black", "bold",
  "2020", "COVID-19",    75,  0.008, "white", "bold",
  "2020", "Circulatory", 60,  0.046, "white", "bold",
  "2020", "Neoplasms",   72,  0.067, "white", "bold"
)

# Callouts (hand-picked coords, at given device dimension, super finnicky)
labs_call_2020 <- tribble(
  ~year,  ~label,        ~x_text, ~y_text, ~x_point, ~y_point, ~hjust,

  # early-life
  "2020", "Perinatal",      3,     0.045,      0,     0.050,     0,
  "2020", "Congenital",     3,     0.037,      0,     0.020,     0,

  # mid/top strata callouts (manual)
  "2020", "Endocrine",     85,     0.098,     78,     0.089,     0,
  "2020", "Respiratory",   85,     0.075,     83,     0.062,     0,

  # these two: right-align at x_text
  "2020", "Infectious",    50,     0.100,     56,     0.091,     1,
  "2020", "Digestive",     50,     0.092,     55,     0.077,     1
)
x_major <- seq(25, 100, 25)
y_major <- seq(0.04, 0.12, 0.04)
## 7) Plot ----
p <- ggplot(all_dec_cod, aes(x = age, y = contribution, fill = cause_plot)) +
  geom_col(
    width = 1,
    linewidth = 0.08,
    aes(color = after_scale(fill)),
    show.legend = FALSE
  ) +

  geom_vline(
    xintercept = x_major,
    linewidth = 0.3,
    color = "white",
    alpha = 0.65
  ) +

  geom_hline(
    yintercept = y_major,
    linewidth = 0.3,
    color = "white",
    alpha = 0.65
  ) +
  facet_grid(year ~ ., switch = "y") +
  scale_fill_manual(
    values = pal,
    breaks = legend_levels,
    labels = function(x) stringr::str_remove(x, " \\([A-Z0-9].*\\)$"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 25),
    minor_breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.12, by = 0.04),
    limits = c(-.01, .12),                      # keep lower at 0
    labels = function(z) sprintf("%.2f", z)     # show 0.10 etc clearly
  ) +
  labs(
    x = "Age",
    y = "Contribution to sex gap (female - male), years",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 13),
    # grid tuning (fainter; keep minor)
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey70"),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0,size=14),

    # legend below
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width  = unit(0.35, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +

  # Callout segments (2020 only)
  geom_segment(
    data = labs_call_2020,
    aes(x = x_text, y = y_text, xend = x_point, yend = y_point),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "grey20"
  ) +
  geom_text(
    data = labs_call_2020,
    aes(x = x_text, y = y_text, label = label, hjust = hjust),
    inherit.aes = FALSE,
    vjust = 0.5,
    size = 3.4,
    color = "grey20"
  ) +
  # Direct labels (2020 only)
  geom_text(
    data = labs_text_2020,
    aes(x = x, y = y, label = label, color = I(col), fontface = fontface),
    inherit.aes = FALSE,
    size = 4.0
  )

p
ggsave("dec_cod.pdf",p,width=15,height=24,units="cm")


#-------------------------------------------------------
# plot for cause totals
# ------------------------------------------------------
df_w <- all_dec_cod |>
  group_by(year, cause_short) |>
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") |>
  mutate(year = as.integer(as.character(year))) |>
  filter(year %in% c(2010, 2015, 2020)) |>
  pivot_wider(
    names_from  = year,
    values_from = contribution,
    names_prefix = "y"
  )

# Order: descending by 2020 (so External etc at top)
df_w <- df_w |>
  arrange(desc(y2020)) |>
  mutate(cause_short = gsub(pattern = "\\(",replacement = "\n(",x=cause_short)) |>
  mutate(cause = factor(cause_short, levels = unique(cause_short)))

# interval table (one row per cause)
df_int <- df_w |>
  transmute(
    cause,
    c10 = y2010, c15 = y2015, c20 = y2020
  ) |>
  mutate(
    y_base = as.integer(cause),
    y_10_15 = y_base + 0.12,
    y_15_20 = y_base - 0.12,

    x0_10_15 = c10, x1_10_15 = c15,
    x0_15_20 = c15, x1_15_20 = c20,

    # net direction (2010 -> 2020)
    dir_net = if_else(c20 >= c10, "Increase", "Decrease")
  )


cols_dir <- c("Increase" = "#5c0084", "Decrease" = "#008444")

# data frames for annotations
lab_external <- df_int |>
  filter(grepl("^External causes", as.character(cause))) |>
  slice(1) |>
  transmute(
    x = c20-.3, y = y_15_20 +.7,
    label = paste0(round(c20-c10,2),"-year increase"),
    col = cols_dir["Increase"]
  )

lab_neopl <- df_int |>
  filter(grepl("^Neoplasms", as.character(cause))) |>
  slice(1) |>
  transmute(
    x = c20+.3, y = y_15_20 + 0.8,
    label = paste0(round(c10-c20,2),"-year decrease"),
    col = cols_dir["Decrease"]
  )

lab_special <- df_int |>
  filter(grepl("^Special codes", as.character(cause))) |>
  slice(1) |>
  transmute(
    x = c20, y = y_15_20 - 0.25,
    label = "(mostly COVID-19)",
    col = cols_dir["Increase"]
  )

head_w <- 0.03   # x-units
head_h <- 0.2    # y-units

df_int2 <- df_int |>
  mutate(
    dir_last = if_else(dir_net == "Increase", "right", "left")
  )

tri_df <- df_int2 |>
  transmute(
    cause, dir_net, dir_last,
    x_butt = x1_15_20,  # 2020 contribution
    x_tip  = if_else(dir_last == "right", x1_15_20 + head_w, x1_15_20 - head_w),
    y_mid  = y_15_20,
    y_top  = y_15_20 + head_h/2,
    y_bot  = y_15_20 - head_h/2
  ) |>
  tidyr::uncount(3, .id = "v") |>
  mutate(
    x = dplyr::case_when(
      v == 1 ~ x_tip,
      v %in% c(2, 3) ~ x_butt
    ),
    y = dplyr::case_when(
      v == 1 ~ y_mid,
      v == 2 ~ y_top,
      v == 3 ~ y_bot
    )
  )

p_step <- ggplot(df_int2) +
  # 2010 -> 2015
  geom_segment(
    aes(x = x0_10_15, xend = x1_10_15, y = y_10_15, yend = y_10_15, color = dir_net),
    linewidth = 0.9,
    lineend = "round"
  ) +
  # vertical connector
  geom_segment(
    aes(x = c15, xend = c15, y = y_10_15, yend = y_15_20, color = dir_net),
    linewidth = 0.9,
    lineend = "round"
  ) +
  # 2015 -> 2020
  geom_segment(
    aes(x = x0_15_20, xend = x1_15_20, y = y_15_20, yend = y_15_20, color = dir_net),
    linewidth = 0.9,
    lineend = "round"
  ) +
  # arrowhead polygon
  geom_polygon(
    data = tri_df,
    aes(x = x, y = y, group = cause, fill = dir_net),
    inherit.aes = FALSE,
    linewidth = 0
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
  scale_color_manual(values = cols_dir, guide = "none") +
  scale_fill_manual(values = cols_dir, guide = "none") +
  scale_y_continuous(
    breaks = df_int2$y_base,
    labels = levels(df_int2$cause),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  #
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = "Cause contribution to sex gap (female - male), years",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 12)
  ) +
  # annotations
  geom_text(
    data = lab_external,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 4,
    color = lab_external$col
  ) +
  geom_text(
    data = lab_neopl,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 4,
    color = lab_neopl$col
  ) +
  geom_text(
    data = lab_special,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 4,
    color = lab_special$col
  )

p_step

ggsave(
  p_step,
  filename = "dec_cod_tot.pdf",
  width = 25, height = 20, units = "cm",
  device = cairo_pdf
)


all_dec_cod$cause_short |> unique()

all_dec_cod |>
  filter(cause_short == "External causes (V01-Y89)" ) |>
  ggplot(aes(x=age, y = Female, color = year, group=year)) +
  geom_line()
