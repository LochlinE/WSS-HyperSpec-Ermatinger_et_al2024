# Use of Morales et al 2021's Interband Redundancy Analysis (IBRA) to identify and remove and noisy spectra due to areas of low emission
# from the halogen lights used for sampling the hyperspectral reflectance



####################### Libraries & Data

library(readxl)
library(tidyverse)

infestationData <- read_xlsx("C:\\Users\\lochl\\OneDrive - Montana State University\\Documents\\SawflyResearch\\Analysis\\ProximalSensing2023\\ProximalInfestationData2023.xlsx")

refs <- read.csv2(file = "C:\\Users\\lochl\\OneDrive - Montana State University\\Documents\\SawflyResearch\\Analysis\\ProximalSensing2023\\IBRA\\ReflectanceData23.csv", header = T, sep = "," )
names(refs)[2:2152] <- c(350:2500)
refs[2:2152] <- lapply(refs[2:2152],as.numeric)

####################### define IBRA 
# translated to the R language as described in https://doi.org/10.3390/rs13183649 (repo: https://github.com/NISL-MSU/HSI-BandSelection)

IBRA <- function(spectral_matrix, threshold) { 
  if (!is.numeric(spectral_matrix)) {
    stop("Input matrix must be numeric.")
  }
  
  vifPair <- function(x, y) {
    model <- lm(x ~ y)
    rsquared <- summary(model)$r.squared
    vif <- 1 / (1 - rsquared)
    return(vif)
  }
  
  n <- ncol(spectral_matrix)
  vif_matrix <- matrix(NA, nrow = n, ncol = n)
  distances_left <- rep(0, n)
  distances_right <- rep(0, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        model <- lm(spectral_matrix[, i] ~ spectral_matrix[, j])
        rsq <- summary(model)$r.squared
        vif <- 1 / (1 - rsq)
        vif_matrix[i, j] <- vif
      } else {
        vif_matrix[i, j] <- NA 
      }
    }
  }
  
  for (band in 1:ncol(vif_matrix)) {
    d_left <- 1
    vifVal_left <- Inf
    while (vifVal_left > threshold && (band - d_left) > 0) {
      if (vif_matrix[band, band - d_left] == 0) {
        vif_matrix[band, band - d_left] <- vifPair(vif_matrix[, band], vif_matrix[, band - d_left])
        vif_matrix[band - d_left, band] <- vif_matrix[band, band - d_left]
      }
      vifVal_left <- vif_matrix[band, band - d_left]
      d_left <- d_left + 1
    }
    distances_left[band] <- d_left - 1
    
    d_right <- 1  
    vifVal_right <- Inf
    while (vifVal_right > threshold && (band + d_right) <= ncol(vif_matrix)) {
      if (vif_matrix[band, band + d_right] == 0) {
        vif_matrix[band, band + d_right] <- vifPair(vif_matrix[, band], vif_matrix[, band + d_right])
        vif_matrix[band + d_right, band] <- vif_matrix[band, band + d_right]
      }
      vifVal_right <- vif_matrix[band, band + d_right]
      d_right <- d_right + 1
    }
    distances_right[band] <- d_right - 1
  }
  
  return(list(vif_matrix = vif_matrix, clusters = abs(distances_left - distances_right)))
}

####################### Identify Noisy Spectra
refsT <- refs %>% 
  filter(DAI == 14) %>%
  dplyr::select(`350`:`2500`) %>%
  as.matrix(.)

IBRAStatic14 <- IBRA(refsT, threshold = 8) # this step will take awhile

refsTMu <- refsT %>%
  as.data.frame(.) %>%
  summarize(across(everything(),list(mean))) %>%
  setNames(350:2500)

IBRAStatic14$reflectance <- unlist(refsTMu)
IBRAStatic14$wavelength <- 350:2500
mat <- as.data.frame(IBRAStatic14[-1])

par(mfrow = c(1,2))
plot(mat$wavelength,mat$clusters, type = "l", xlab = "Wavelength nm", ylab = "VIF", main = "IBRA")
abline(v = 520, col = "red")
abline(v = 1870, col = "red")
plot(mat$wavelength,mat$reflectance, type = "l", col = "green", ylim = c(0.04,.4),xlab = "Wavelength nm", ylab = "reflectance", main = "Mean Reflectance")
abline(v = 520, col = "red")
abline(v = 1870, col = "red")
