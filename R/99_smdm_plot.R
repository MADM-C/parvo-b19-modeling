library(ggplot2)
library(viridis)
df_smdm <- data.frame(
  strategy = c(rep("Screening", 3), rep("Vaccination", 3), rep("Screening + Vaccination", 3)),
  outcome = rep(c("B19 Fetal Deaths", "B19 Stillbirths", "Intrauterine Transfusions"), 3),
  value = c(482, 482, -403, 5875, 651, 1288, 6140, 916, 1066)
)

df_smdm$strategy <- factor(df_smdm$strategy,
                           levels = c("Screening", "Vaccination", "Screening + Vaccination"))

ggplot(data = df_smdm, aes(fill = strategy, y = value, x = outcome)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values=c("#1AE4B6FF", "#FABA39FF", "#7A0403FF")) +
  ylab("Events Averted") +
  xlab("Outcome") +
  labs(fill = "Strategy")

ggsave(path = "results/Phase 3 Plots", filename = "base_case_bar.png",
       width = 9, height = 5)
