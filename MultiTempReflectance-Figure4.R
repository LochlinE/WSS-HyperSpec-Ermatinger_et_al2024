# Code used to create the plot shown in Figure 4

####################### Libraries & Data
library(readxl)
library(dplyr)
library(tidyr)
library(plotly)
library(extrafont)

refs <- read.csv2(file = "Data/ReflectanceData23.csv", header = T, sep = "," )
names(refs)[2:2152] <- c(350:2500)
refs[2:2152] <- lapply(refs[2:2152],as.numeric)

# filter out wavelengths below 520 nm & above 1870 nm.
refsLongFull <- refs %>% select(c(`520`:`1870`,DAI,sampleName,plantName)) %>%
  pivot_longer(-c(DAI,sampleName,plantName),
               names_to = "wavelength",
               values_to = "reflectance")
refsLong$wavelength <- as.numeric(refsLong$wavelength)

specMu <- refsLongFull %>%
  group_by(DAI,wavelength) %>%
  summarise(mean_reflectance = mean(reflectance)) %>%
  ungroup() %>%
  setNames(c("time","wavelength","mean_reflectance"))

# Convert mean_df to a wide format
mean_wide_df <- specMu %>%
  pivot_wider(names_from = wavelength, values_from = mean_reflectance)

# Create a grid of unique time and wavelength values
time <- sort(unique(specMu$time))
wavelength <- sort(unique(specMu$wavelength))
Z <- matrix(NA, nrow = length(time), ncol = length(wavelength))
for (i in 1:length(time)) {
  for (j in 1:length(wavelength)) {
    Z[i, j] <- specMu %>% 
      filter(time == time[i], wavelength == wavelength[j]) %>% 
      pull(mean_reflectance)
  }
}
refsLongFull$wavelength <- as.numeric(refsLongFull$wavelength)
mean_df <- refsLongFull %>%
  group_by(DAI, wavelength) %>%
  summarise(mean_reflectance = mean(reflectance)) %>%
  spread(key = wavelength, value = mean_reflectance)

# Extract time and wavelength values
time <- unique(mean_df$DAI)
wavelength <- names(mean_df)[-1]

# Create Z matrix 
Z <- as.matrix(mean_df[, -1])
surfacecolor <- t(matrix(rep(time, length(wavelength)), nrow = length(wavelength), byrow = TRUE))

plot_ly(z = Z,y = time, x = wavelength, type = "surface", 
        surfacecolor = surfacecolor,
        colors = c("#33CC33","#33CC33","greenyellow","greenyellow","yellow")) %>%
  layout(scene = list(
    xaxis = list(title = "Wavelength nm", range = c(520,1870)),
    yaxis = list(title = "Days After Infestation"),
    zaxis = list(title = "Mean Reflectance", range = c(0,.6))),
    font = list(family = 'Times New Roman'),
    title = NULL)

