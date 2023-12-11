library(readxl)
library(reshape2)
library(ggplot2)

# Read the Excel file
AcidityResults <- read_excel("AcidityResults.xlsx", col_names = TRUE)

well_watered <- rowMeans(AcidityResults[, c("C1", "C2", "C3")])
salt_treated <- rowMeans(AcidityResults[, c("S1", "S2", "S3")])

# Create a new data frame with the original values and the corresponding day
averages <- data.frame(Day = factor(AcidityResults$Day, levels = AcidityResults$Day),
                       Well_Watered = c(AcidityResults$C1, AcidityResults$C2, AcidityResults$C3),
                       Salt_Treated = c(AcidityResults$S1, AcidityResults$S2, AcidityResults$S3))

# Reshape the data frame to a longer format
averages_long <- reshape2::melt(averages, id.vars = "Day", variable.name = "Group", value.name = "Value")

# Create the boxplot using ggplot
plot <- ggplot(averages_long, aes(x = Day, y = Value, fill = Group)) +
  geom_boxplot() +
  labs(title = "Titratable acidity",
       x = "Days of salt treatment", y = "μmol[H+].g-1") +
  scale_fill_manual(values = c("Well_Watered" = "blue", "Salt_Treated" = "red")) +
  theme_bw()

ggsave("images/boxplot_acidity.png", plot, dpi = 300, width=9, height=8)
