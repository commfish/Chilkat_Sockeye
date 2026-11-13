# run the script 1_coho_apportionment_ratio.R first
# figure of daily coho counts at the Chilkat weir, highlighting 2025
# reads coho_daily_wide_count.csv (daily counts after apportionment with the imputed ratios)
# and plots 2025 against 2018-2024, shading the 2025 days filled by imputation
# the file coho_daily_wide_count.csv needs to be created based on the data file located
# here: C:\chilkat_sockeye\coho_apportionment\data\raw_files

#  load ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
source('coho_apportionment/code/missfill_functions.r') # for tickr()

# settings ----
gap_start <- 236   # first imputed day in 2025
gap_end   <- 250   # last imputed day in 2025
focus_year <- 2025

theme_sleek <- function(base_size = 12, base_family = "serif") {
  half_line <- base_size / 2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.key.size = unit(2, "lines"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}
theme_set(theme_sleek())

# data ----
wide <- read.csv('coho_apportionment/output/coho_daily_wide_count.csv')

long <- wide %>%
  pivot_longer(
    cols = -JD,
    names_to = "year",
    values_to = "total_count",
    names_transform = list(year = readr::parse_number)) %>%
  rename(date = JD) %>%
  select(year, date, total_count) %>%
  arrange(year, date)

write.csv(long, "coho_apportionment/output/coho_daily_long_count.csv", row.names = FALSE)

# figure ----
other_label <- paste0(min(long$year), "-", focus_year - 1)   # e.g. "2018-2024"
plot_data <- long %>%
  filter(!is.na(total_count)) %>%                            # drop out-of-season days
  mutate(group_label = ifelse(year == focus_year, as.character(focus_year), other_label)) %>%
  arrange(year == focus_year)                                # draw focus year last (on top)

x_range <- c(160, 300)
axisf <- tickr(data.frame(date = x_range[1]:x_range[2]), date, 20)

fig <- ggplot(plot_data, aes(x = date, y = total_count, color = group_label, group = year)) +
  annotate("rect", xmin = gap_start, xmax = gap_end, ymin = -Inf, ymax = Inf,
           fill = "grey50", alpha = 0.15) +
  geom_vline(xintercept = c(gap_start, gap_end), linetype = "dashed", color = "black") +
  geom_line() +
  scale_x_continuous(breaks = axisf$breaks, labels = axisf$labels) +
  coord_cartesian(xlim = x_range) +
  scale_color_manual(values = setNames(c("red", "grey60"),
                                       c(as.character(focus_year), other_label))) +
  labs(x = "Julian day", y = "Coho count", color = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.8))

ggsave("coho_apportionment/output/coho_missfill.png", plot = fig,
       dpi = 500, height = 4, width = 7, units = "in")
