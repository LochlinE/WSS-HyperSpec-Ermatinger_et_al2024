
library(readxl)
library(dplyr)
library(ggplot2)
library(ggridges)

infestationData <- read_xlsx("C:\\Users\\lochl\\OneDrive - Montana State University\\Documents\\SawflyResearch\\Analysis\\ProximalSensing2023\\ProximalInfestationData2023.xlsx")

modeledData <- infestationData %>%
  filter(rep %in% c("G","H","I","J")) %>%
  filter(sampleName != "PII05") %>%
  mutate(plantingGroup = case_match(rep,
                                    "G" ~ "A",
                                    "H" ~ "B",
                                    "I" ~ "C",
                                    "J" ~ "D"))
combinedData <- modeledData %>%
  mutate(plantingGroup = "Combined Groups")
combinedModeledData <- bind_rows(modeledData, combinedData)

# Set the order of the groups
combinedModeledData$plantingGroup <- factor(combinedModeledData$plantingGroup, 
                                            levels = c("A", "B", "C", "D", "Combined Groups"))
# Define custom colors
custom_colors <- c("A" = "#74c476", "B" = "#74c476", "C" = "#74c476", "D" = "#74c476", "Combined Groups" = "#006d2c")

# Create the ridge density plot
FigureA3 <- ggplot(combinedModeledData, aes(x = percentSigInfested, y = plantingGroup, fill = plantingGroup)) +
  geom_density_ridges(scale = 1, alpha = .75) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = "Times New Roman", size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),) +
  labs( x = "Proportion of Adequately Infested WSS Stems", y = "") +
  xlim(0, NA) 