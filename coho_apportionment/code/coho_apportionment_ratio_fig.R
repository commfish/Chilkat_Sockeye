# using the proportions, create the coho_daily_wide_count.csv file

# load ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# create a figure 
tickryr <- data.frame(date = 160:300)
axisf <- tickr(tickryr, date, 20)

theme_sleek <- function(base_size = 12, base_family = "Times") {
  half_line <- base_size/2
  theme_light(base_size = 12, base_family = "Times") +
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

read.csv('coho_apportionment/data/coho_daily_wide_count.csv')-> wide

long <- wide %>%
  pivot_longer(
    cols            = -JD,
    names_to        = "year",
    values_to       = "total_count",
    names_transform = list(year = readr::parse_number)
  ) %>%
  rename(date = JD) %>%
  select(year, date, total_count) %>%
  arrange(year, date)%>%
  as.data.frame

write.csv(long, "coho_apportionment/output/global_imputed_fig.csv")

theme_set(theme_sleek())
read.csv("coho_apportionment/output/global_imputed_fig.csv") %>%
  ggplot(aes(x = date, y = total_count,
             color = ifelse(year == 2025, "2025", "2018-2024"),
             group = year)) +
  geom_line() +
  scale_x_continuous(limits = c(min(tickryr$date), max(tickryr$date)),
                     breaks = axisf$breaks, labels = axisf$labels) +
  scale_color_manual(values = c("2025" = "red", "2021-2024" = "grey")) +
  labs(x = "Julian day", y = "Coho count", color = NULL) +
  theme_classic() +
  theme(legend.position = c(0.2, 0.8))+
  geom_vline(xintercept = 236, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 250, linetype = "dashed", color = "black")-> fig

ggsave(paste0("coho_apportionment/output/coho_missfill.png"), plot = fig, dpi = 500, height = 4, width = 7, units = "in")




