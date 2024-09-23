# This script contains the Sparse Multiway Partial Least Squares Regression (sNPLS) analysis
# and accompanying figures 5-8. Because the figures rely on data from the model, it is most convenient 
# to store the code for these figures immediately following the creation of the model in a single script.

####################### Script Index
# Scale and center reflectance data - lines: 38 - 62
# Calibration and Validation data partition - lines: 65 - 159
# Fit and predict sNPLS - lines 172 - 182
# Figures 5 - lines: 216 - 267
# Figure 6 - lines: 269 - 323
# beta coefficients - lines: 326 - 357
# Figure 7 - lines: 359 - 388
# Figure 8 - lines: 390 - 696

####################### Libraries & Data
library(readxl)
library(dplyr)
library(tidyr)
library(sNPLS)
library(ggplot2)
library(reshape2)
library(extrafont)

infestationData <- read_xlsx("Data/ProximalInfestationData2023.xlsx")
refs <- read.csv2(file = "Data/ReflectanceData23.csv", header = T, sep = "," )
names(refs)[2:2152] <- c(350:2500)
refs[2:2152] <- lapply(refs[2:2152],as.numeric)

# filter out wavelengths below 520 nm & above 1870 nm.
refsLong <- refs %>% 
  dplyr::select(c(`520`:`1870`,DAI,sampleName,plantName)) %>%
  pivot_longer(-c(DAI,sampleName,plantName),
               names_to = "wavelength",
               values_to = "reflectance")
refsLong$wavelength <- as.numeric(refsLong$wavelength) # convert to numeric mostly for plotting

# make a wide version
refsWide <- refsLong %>%
  pivot_wider(names_from = wavelength,
              values_from = reflectance,
              id_cols = c(sampleName,DAI))

####################### Scale and Center Reflectance data by DAI
df <- refsLong %>% 
  inner_join(infestationData, by = "plantName") %>%
  mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
  filter(plantName != "PII05") %>% # missed DAI 0
  filter(rep != "K") %>% # incomplete planting group
  dplyr::select(DAI,plantName,wavelength,reflectance,sigInf.Prop)

refsLongUnfold <- df %>% unite("DAI_Wavelength", c(DAI,wavelength), sep = "_") %>%
  pivot_wider(names_from = DAI_Wavelength,
              values_from = "reflectance")
refsLongUnfold.Scaled <- scale(refsLongUnfold[,3:length(refsLongUnfold)], center = T, scale = T)
refsWideScaled <- as.data.frame(cbind(refsLongUnfold[,1:2],refsLongUnfold.Scaled))

refsLongScaled <- refsWideScaled %>%
  pivot_longer(
    cols = -c(plantName,sigInf.Prop),
    names_to = "DAI_wavelength",
    values_to = "reflectance")

df.long.scaled <- refsLongScaled %>%
  separate(
    DAI_wavelength,
    into = c("DAI", "wavelength"),
    sep = "_")

####################### Partition Calibration Data and Populate Tensor

infNamesARD <- infestationData %>%
  filter(rep != "K") %>%
  filter(plantName != "PII05") %>%
  mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
  select(plantName,sigInf.Prop) %>% 
  arrange(desc(sigInf.Prop))

TestNames <- infNamesARD$plantName[seq(from = 1,
                                       to = length(infNamesARD$plantName),
                                       by = 7)]

df <- df.long.scaled %>% 
  inner_join(infestationData, by = "plantName") %>%
  mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
  filter(rep != "K") %>% # incomplete planting group
  filter(plantName != "PII05") %>% # missed DAI 0
  filter(!plantName %in% TestNames ) %>%
  select(DAI,plantName,wavelength,reflectance,sigInf.Prop)

# Sort the unique values for indexing
plant_names <- sort(unique(df$plantName))
dai <- sort(unique(df$DAI))
wavelengths <- sort(unique(df$wavelength))

# Populate tensor
n_plants <- length(plant_names)
n_dai <- length(dai)
n_wavelengths <- length(wavelengths)
data_array <- array(NA, dim = c(n_plants, n_dai, n_wavelengths))
for (i in 1:n_plants) {
  for (j in 1:n_dai) {
    plant_df <- df[df$plantName == plant_names[i] & df$DAI == dai[j], ]
    for (k in 1:n_wavelengths) {
      reflectance_val <- plant_df$reflectance[plant_df$wavelength == wavelengths[k]]
      if (length(reflectance_val) > 0) {
        data_array[i, j, k] <- reflectance_val
      }
    }
  }
}

X.Train <- data_array

# dependent (y) data array [I,1]
Y.NPLS.train <- matrix(unlist(infestationData %>%
                                filter(rep != "K") %>%
                                filter(plantName != "PII05") %>%
                                filter(!plantName %in% TestNames) %>%
                                mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
                                select(sigInf.Prop)), ncol = 1)

# Transform to logit
Y.NPLS.Train.Logit <- log((Y.NPLS.train+0.01) / (1 - (Y.NPLS.train+0.01)))

####################### Partition Validation Data and Populate Tensor

# recreate df to grab names not found in "testnames"
df <- df.long.scaled %>% 
  inner_join(infestationData, by = "plantName") %>%
  mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
  filter(rep != "K") %>% # incomplete planting group
  filter(plantName != "PII05") %>% # missed DAI 0
  filter(plantName %in% TestNames ) %>%
  select(DAI,plantName,wavelength,reflectance,sigInf.Prop)

# Sort the unique values for indexing
plant_names <- sort(unique(df$plantName))
dai <- sort(unique(df$DAI))
wavelengths <- sort(unique(df$wavelength))

# Populate tensor
n_plants <- length(plant_names)
n_dai <- length(dai)
n_wavelengths <- length(wavelengths)
data_array <- array(NA, dim = c(n_plants, n_dai, n_wavelengths))
for (i in 1:n_plants) {
  for (j in 1:n_dai) {
    plant_df <- df[df$plantName == plant_names[i] & df$DAI == dai[j], ]
    for (k in 1:n_wavelengths) {
      reflectance_val <- plant_df$reflectance[plant_df$wavelength == wavelengths[k]]
      if (length(reflectance_val) > 0) {
        data_array[i, j, k] <- reflectance_val
      }
    }
  }
}

X.Test <- data_array

# dependent (y) data array [I,1]
Y.NPLS.test <- matrix(unlist(infestationData %>%
                               filter(rep != "K") %>%
                               filter(plantName != "PII05") %>%
                               filter(plantName %in% TestNames) %>%
                               mutate(sigInf.Prop = totalSignificantInfestations/totalStems) %>%
                               select(sigInf.Prop)), ncol = 1)
# Transform to logit
Y.NPLS.Test.Logit <- log((Y.NPLS.test+0.01) / (1 - (Y.NPLS.test+0.01)))

####################### Sparse Multiway Partial Least Squares Regression Model
train.NPLS6 <- sNPLS(XN = X.Train,
                     Y = Y.NPLS.Train.Logit,
                     scale.X = F, # Scaling & center was done slab wise (based on DAI), thus we dont use the internal functions
                     center.X = F,
                     scale.Y = T,
                     center.Y = T,
                     keepJ = rep(1:ncol(X.Train),13),
                     keepK = rep(1:dim(X.Train)[3],13),
                     ncomp = 13)

preds6.2 <- predict(train.NPLS6, X.Test, rescale = T) # this will take awhile

# Define the inverse logit function
inv_logit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

#  Apply the inverse logit function
modelOG <- inv_logit(Y.NPLS.Train.Logit)
model.fitOG <- inv_logit(train.NPLS6$Yadj)
testOG <- inv_logit(Y.NPLS.Test.Logit)
predsOG <- inv_logit(preds6.2)

# Subtract 0.01 to get back the original Y.NPLS.test values
modelOG <- modelOG - 0.01
model.fitOG <- model.fitOG - 0.01
testOG <- testOG - 0.01
predsOG <- predsOG - 0.01

calculate_metrics <- function(Y_true, Y_pred) { 
  residuals <- Y_true - Y_pred
  rss <- sum(residuals^2)
  tss <- sum((Y_true - mean(Y_true))^2)
  r_squared <- 1 - rss/tss
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  return(list(R2 = r_squared, RMSE = rmse, MAE = mae))
}

####################### Figure 5

# The data contained in 'ModelCalibrationPlot' were calculated by rerunning the code above with each number of latent variables to determine the optimal number for the final model
# RMSE and the coefficient of determination (R2) was calculated using the defined function on lines 201 - 209 'calculate_metrics'
ModelCalibrationPlot <- 
  data.frame(cbind(
    components <- 1:20,
    calRMSE <- c(.21,.21,.2,.19,.19,
                 .18,.16,.15,.13,.13,
                 .13,.12,.13,.12,.12,
                 .11,.11,.10,.10,.10),
    valRMSE <- c(.26,.26,.26,.24,.25,
                 .26,.22,.23,.21,.17,
                 .17,.14,.14,.15,.22,
                 .22,.19,.24,.23,.28))) %>%
  setNames(c("Component","Calibration RMSE","Validation RMSE"))

# Melt the data frame for plotting
ModelCalibrationPlotMelt <- melt(ModelCalibrationPlot, id.vars = "Component")

# Plot the data
ggplot(ModelCalibrationPlotMelt, aes(x = Component, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  labs(title = "RMSE vs. Number of Components",
       x = "Number of Components",
       y = "RMSE",
       color = "Series") +
  theme_bw()

Figure5 <- ggplot(ModelCalibrationPlotMelt, aes(x = Component, y = value, color = variable)) +
  geom_line(size = .75) +
  geom_point(size = 1.75) +
  geom_vline(xintercept = 13, linetype = "dashed", size = .75) +
  scale_color_manual(values = c("Calibration RMSE" = "#33CC33", "Validation RMSE" = "purple")) +
  labs(x = "Number of Components",
       y = "RMSE",
       color = "Series") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.position = c(0.75, 0.9),
    legend.background = element_rect(fill = NA, color = NA))

####################### Figure 6

Train <- as.data.frame(cbind(modelOG,model.fitOG)) %>% 
  setNames(c("Observed","Predicted")) %>%
  mutate(Class = "calibration")
Test <- as.data.frame(cbind(testOG,predsOG)) %>% 
  setNames(c("Observed","Predicted")) %>%
  mutate(Class = "validation")
Cumulative <- as.data.frame(rbind(Train,Test))

# Function to calculate R-squared
calc_r_squared <- function(observed, predicted) {
  ss_total <- sum((observed - mean(observed))^2)
  ss_residual <- sum((observed - predicted)^2)
  r_squared <- 1 - (ss_residual / ss_total)
  return(r_squared)
}

# Calculate RMSE, Bias, R-squared, and SE
rmse_train <- sqrt(mean((Train$Predicted - Train$Observed)^2))
bias_train <- mean(Train$Predicted - Train$Observed)
r_squared_train <- calc_r_squared(Train$Observed, Train$Predicted)
se_train <- sd(Train$Predicted - Train$Observed) / sqrt(nrow(Train))

rmse_test <- sqrt(mean((Test$Predicted - Test$Observed)^2))
bias_test <- mean(Test$Predicted - Test$Observed)
r_squared_test <- calc_r_squared(Test$Observed, Test$Predicted)
se_test <- sd(Test$Predicted - Test$Observed) / sqrt(nrow(Test))

Cumulative <- rbind(
  Train %>% mutate(Class = "calibration"),
  Test %>% mutate(Class = "validation"))

# Plotting with annotations
Figure6 <- ggplot(Cumulative, aes(x = Predicted, y = Observed, color = Class)) +
  stat_smooth(data = subset(Cumulative, Class == "calibration"), method = "lm", se = TRUE, aes(group = 1), color = "#33CC33", linetype = "dashed", size = 3) +
  geom_point(shape = 19, size = 5) +
  geom_abline() +
  xlim(0, 1) + ylim(0, 1) +
  scale_color_manual(values = c("validation" = "purple", "calibration" = "#33CC33")) +
  labs(x = "Predicted Proportion of Adequate WSS Infestation",
       y = "") +
  theme_bw() + 
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.90, 0.15),
    legend.text = element_text(size = 20),
    legend.background = element_rect(fill = NA, color = NA),
    axis.title.y = element_text(margin = margin(t = 40, r = 0, b = 0, l = 0))) +
  annotate("text", x = 0.0, y = 0.90, label = paste0("RMSE = ", round(rmse_train, 3)*100, "%, R² = ", round(r_squared_train, 3), ", SE = ", round(se_train, 3), ", bias = ", round(bias_train, 3)), color = "#33CC33", size = 8, hjust = 0, family = "Times New Roman") +
  annotate("text", x = 0.0, y = 0.95, label = paste0("RMSE = ", round(rmse_test, 3)*100, "%, R² = ", round(r_squared_test, 3), ", SE = ", round(se_test, 3), ", bias = ", round(bias_test, 3)), color = "purple", size = 8, hjust = 0, family = "Times New Roman")

####################### Figure 7

# Grab beta coefficents from model
N.PLS6.Coefs <- coef(train.NPLS6, as.matrix = T) # this will take awhile

mat.T <- as.data.frame(as.matrix(rbind(
  unname(N.PLS6.Coefs[1,,]),
  unname(N.PLS6.Coefs[2,,]),
  unname(N.PLS6.Coefs[3,,]),
  unname(N.PLS6.Coefs[4,,]),
  unname(N.PLS6.Coefs[5,,]),
  unname(N.PLS6.Coefs[6,,]),
  unname(N.PLS6.Coefs[7,,]),
  unname(N.PLS6.Coefs[8,,]),
  unname(N.PLS6.Coefs[9,,])))) %>%
  setNames(520:1870) %>%
  mutate(DAI = c(0,7,14,21,28,35,42,49,56))

mat.T.Long <- mat.T %>%
  pivot_longer(-DAI,names_to = "wavelength", values_to = "beta")

# reshaped for a surface
mat.T3D <- list(as.matrix(mat.T[,1:1351]), #z
                as.numeric(colnames(mat.T[,1:1351])), #x
                mat.T$DAI) %>% #y
  setNames(c("z","x","y"))

# Convert matrix to long format for ggplot2
mat_long <- melt(mat.T3D$z) %>%
  setNames(c("DAI","wavelength","beta"))

# Create a named vector using setNames
week_mapping <- setNames(c(0, 7, 14, 21, 28, 35, 42, 49, 56), 1:9)
mat_long <- mat_long %>%
  mutate(DAI = recode(DAI, !!!week_mapping))

Figure7 <- ggplot(mat_long, aes(x = wavelength, y = DAI, fill = sign(beta))) +
  geom_tile() +
  scale_fill_gradient2(low = "#663300", mid = "white", high = "#009999", midpoint = 0) +
  scale_x_continuous(
    breaks = seq(500, 1800, by = 100), 
    labels = c("500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800"),
    sec.axis = sec_axis(
      trans = ~ .,  # Keep the same scale as the primary X-axis
      breaks = c(575, 940, 1300, 1385, 1550, 1650),  # Specific tick positions for the top axis
      labels = c("575", "940", "1300", "1385", "1550", "1650"),  # Labels corresponding to the breaks
      name = "Relevant Cluster Centers, Wavelength nm"  # Title is set to NULL to avoid redundancy
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 56, by = 7), 
    labels = c("0", "7", "14", "21", "28", "35", "42", "49", "56")
  ) +
  labs(x = "Wavelength nm", y = "Days after Infestation") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 18),
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.title.x.top = element_text(size = 20),
    axis.text.x.top = element_text(size = 18),  # Adjust the font size of the top X-axis text if needed
    axis.text.y = element_text(size = 18),
    panel.grid = element_blank(),
    legend.position = "none")

####################### Figure 8

### 0 DAI
specMuMod <- refsLong %>%
  group_by(DAI,wavelength) %>%
  summarise(mean_reflectance = mean(reflectance)) %>%
  ungroup() %>%
  setNames(c("DAI","wavelength","mean_reflectance"))
specMuMod$wavelength <- as.numeric(specMuMod$wavelength)
spec0 <- specMuMod %>%
  filter(DAI == 0)

mat.T.Long0 <- mat.T.Long %>%
  filter(DAI == "0")

p0 <- ggplot(mat.T.Long0, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec0$mean_reflectance * 4, color = "Reflectance"), size = 1) +
  scale_y_continuous(
    name = "Beta Coefficients",
    limits = c(-0.5, 1.8),
    sec.axis = sec_axis(~ . / 4,
                        name = "Reflectance",
                        breaks = seq(0, .4, by = 0.1)  
    )
  ) +
  geom_line(aes(y = mat.T.Long0$beta, color = "Beta Coefficients"), size = 1) +
  scale_x_continuous(
    name = NULL,  
    sec.axis = sec_axis(~ ., 
                        name = "Wavelength nm",  
                        breaks = seq(500, 1850, by = 250)  
    )
  ) +
  scale_color_manual(values = c("Reflectance" = "#33CC33", "Beta Coefficients" = "black")) +
  labs(x = NULL) + 
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.left = element_text(size = 20),
    axis.text.y.left = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.top = element_text(size = 20),  
    axis.text.x.top = element_text(size = 16),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = 1700, y = 1.7, label = "0 DAI", size = 8, family = "Times New Roman")

### 7 DAI
spec7 <- specMuMod %>%
  filter(DAI == 7)

mat.T.Long7 <- mat.T.Long %>%
  filter(DAI == "7")

p7 <- ggplot(mat.T.Long7, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec7$mean_reflectance * 4, color = "Reflectance"), size = 1) +
  scale_y_continuous(
    name = "Beta Coefficients",
    limits = c(-0.5, 1.8),
    sec.axis = sec_axis(~ . / 4,
                        name = "Reflectance",
                        breaks = seq(0, .4, by = 0.1)  
    )
  ) +
  geom_line(aes(y = mat.T.Long7$beta, color = "Beta Coefficients"), size = 1) +
  scale_x_continuous(
    name = NULL,  
    sec.axis = sec_axis(~ ., 
                        name = "Wavelength nm", 
                        breaks = seq(500, 1850, by = 250)  
    )
  ) +
  scale_color_manual(values = c("Reflectance" = "#33CC33", "Beta Coefficients" = "black")) +
  labs(x = NULL) +  
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.top = element_text(size = 20),  
    axis.text.x.top = element_text(size = 16), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = 1700, y = 1.7, label = "7 DAI", size = 8, family = "Times New Roman")

### 14 DAI
spec14 <- specMuMod %>%
  filter(DAI == 14)

mat.T.Long14 <- mat.T.Long %>%
  filter(DAI == "14")

p14 <- ggplot(mat.T.Long14, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec14$mean_reflectance * 4, color = "Reflectance"), size = 1) +
  scale_y_continuous(
    name = "Beta Coefficients",
    limits = c(-0.5, 1.8),
    sec.axis = sec_axis(~ . / 4,
                        name = "Reflectance",
                        breaks = seq(0, .4, by = 0.1)
    )
  ) +
  geom_line(aes(y = mat.T.Long14$beta, color = "Beta Coefficients"), size = 1) +
  scale_x_continuous(
    name = NULL,  
    sec.axis = sec_axis(~ .,  
                        name = "Wavelength nm",
                        breaks = seq(500, 1850, by = 250) 
    )
  ) +
  scale_color_manual(values = c("Reflectance" = "#33CC33", "Beta Coefficients" = "black")) +
  labs(x = NULL) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title.y.right = element_text(size = 20),
    axis.text.y.right = element_text(size = 16),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.top = element_text(size = 20),  
    axis.text.x.top = element_text(size = 16),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = 1700, y = 1.7, label = "14 DAI", size = 8, family = "Times New Roman")

### 21 DAI
spec21 <- specMuMod %>%
  filter(DAI == 21)

mat.T.Long21 <- mat.T.Long %>%
  filter(DAI == "21")

p21 <- ggplot(mat.T.Long21, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec21$mean_reflectance*4, color = "Mean Reflectance"), size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Mean Reflectance")) + 
  geom_line(aes(y = mat.T.Long21$beta, color = "Beta Coefficients"), size = 1) +  
  scale_color_manual(values = c("Mean Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_text(size = 20),
        axis.text.y.left = element_text(size = 16),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  annotate("text", x = 1700, y = 1.7, label = "21 DAI", size = 8, family = "Times New Roman")

### 28 DAI
spec28 <- specMuMod %>%
  filter(DAI == 28)

mat.T.Long28 <- mat.T.Long %>%
  filter(DAI == "28")

p28 <- ggplot(mat.T.Long28, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec28$mean_reflectance*4, color = "Mean Reflectance"),size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Mean Reflectance")) + 
  geom_line(aes(y = mat.T.Long28$beta, color = "Beta Coefficients"), size = 1) +  
  scale_color_manual(values = c("Mean Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  annotate("text", x = 1700, y = 1.7, label = "28 DAI", size = 8, family = "Times New Roman")

### 35 DAI
spec35 <- specMuMod %>%
  filter(DAI == 35)

mat.T.Long35 <- mat.T.Long %>%
  filter(DAI == "35")

p35 <- ggplot(mat.T.Long35, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec35$mean_reflectance*4, color = "Reflectance"), size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Reflectance",breaks = seq(0, .4, by = 0.1))) + 
  geom_line(aes(y = mat.T.Long35$beta, color = "Beta Coefficients"), size = 1) +  
  scale_color_manual(values = c("Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        axis.title.y.right = element_text(size = 20),
        axis.text.y.right = element_text(size = 16),
        axis.title.y.left = element_blank(), 
        axis.text.y.left = element_blank(),
        axis.title.x = element_blank(),  
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  annotate("text", x = 1700, y = 1.7, label = "35 DAI", size = 8, family = "Times New Roman")

### 42 DAI
spec42 <- specMuMod %>%
  filter(DAI == 42)

mat.T.Long42 <- mat.T.Long %>%
  filter(DAI == "42")

p42 <- ggplot(mat.T.Long42, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec42$mean_reflectance*4, color = "Mean Reflectance"), size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Mean Reflectance")) + 
  scale_x_continuous(name = "Wavelength nm", breaks = seq(500, 1850, by = 250)) +
  geom_line(aes(y = mat.T.Long42$beta, color = "Beta Coefficients"), size = 1) +  
  scale_color_manual(values = c("Mean Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        axis.title.y.right = element_blank(),  
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_text(size = 20),
        axis.text.y.left = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  annotate("text", x = 1700, y = 1.7, label = "42 DAI", size = 8, family = "Times New Roman")

### 49 DAI
spec49 <- specMuMod %>%
  filter(DAI == 49)

mat.T.Long49 <- mat.T.Long %>%
  filter(DAI == "49")

p49 <- ggplot(mat.T.Long49, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec49$mean_reflectance*4, color = "Mean Reflectance"),size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Mean Reflectance")) + 
  scale_x_continuous(name = "Wavelength nm", breaks = seq(500, 1850, by = 250)) +
  geom_line(aes(y = mat.T.Long49$beta, color = "Beta Coefficients"), size = 1) +  
  scale_color_manual(values = c("Mean Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),  
        axis.text.y.left = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  annotate("text", x = 1700, y = 1.7, label = "49 DAI", size = 8, family = "Times New Roman")

### 56 DAI
spec56 <- specMuMod %>%
  filter(DAI == 56)

mat.T.Long56 <- mat.T.Long %>%
  filter(DAI == "56")

p56 <- ggplot(mat.T.Long56, aes(x = as.numeric(wavelength))) +
  geom_line(aes(y = spec56$mean_reflectance*4, color = "Reflectance"),size = 1) +
  scale_y_continuous(name = "Beta Coefficients",limits = c(-.5,1.8), sec.axis = sec_axis(~./4, name = "Reflectance",breaks = seq(0, .4, by = 0.1))) + 
  scale_x_continuous(name = "Wavelength nm", breaks = seq(500, 1850, by = 250)) +
  geom_line(aes(y = mat.T.Long56$beta, color = "Beta Coefficients"),size = 1) +  
  scale_color_manual(values = c("Reflectance" = "#33CC33", "Beta Coefficients" = "black")) + 
  labs(x = "Wavelength nm") + 
  theme_bw() + 
  theme(text = element_text(family = "Times New Roman"),
        # axis.title.y.right = element_blank(),  
        # axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank(),  
        axis.text.y.left = element_blank(),
        axis.title.y.right = element_text(size = 20),
        axis.text.y.right = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + 
  annotate("text", x = 1700, y = 1.7, label = "56 DAI", size = 8, family = "Times New Roman")

### Patch work to create panel plot 
library(patchwork)

Figure8 <- (p0 | p7 | p14) /
  (p21 | p28 | p35) /
  (p42 | p49 | p56) + theme(text = element_text(family = "Times New Roman"))

