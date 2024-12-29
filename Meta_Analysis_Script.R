## Created by : Mila Desi Anasanti
## email: mila.mld@nusamandiri.ac.id

# Set Working Directory
setwd("C:/Users/user1/Downloads/nasir")

# Load Required Libraries
library(meta)
library(metafor)
library(knitr)

# Load Dataset
data <- read.csv("cleaned_subgroup_analysis_data_for_manuscript.csv")  # Replace with your file path

# Inspect Dataset
str(data)
head(data)

# Data Cleaning: Remove Missing Core Meta-Analysis Variables
data <- data[!is.na(data$SMD) & !is.na(data$SE), ]

# Meta-Analysis -------------------------------------------------------------
m1 <- metagen(
  TE = SMD,                # Standardized Mean Difference
  seTE = SE,               # Standard Error
  studlab = Study,         # Study Label
  data = data,             # Dataset
  comb.fixed = TRUE,       # Fixed-effect model
  comb.random = TRUE,      # Random-effects model
  method.tau = "REML",     # Restricted Maximum Likelihood
  hakn = TRUE,             # Hartung-Knapp adjustment
  prediction = TRUE,       # Prediction interval
  sm = "SMD"               # Effect size type
)

# Summary of Meta-Analysis
summary(m1)

# Heterogeneity Statistics
cat("Heterogeneity:\n")
cat("Q Statistic:", m1$Q, "\n")
cat("I² Statistic:", m1$I2, "%\n")
cat("Tau²:", m1$tau2, "\n")

# Forest Plot
forest(m1, 
       common = FALSE,             # Use random-effects model only
       digital.sd = 2,             # Display numbers with 2 decimal places
       sortvar = TE,               # Sort studies by effect size (SMD)
       fontsize = 10,              # Font size for clarity
       col.diamond = "blue",       # Blue color for the overall effect diamond
       col.diamond.lines = "blue") # Blue border for the diamond

# Funnel Plot and Egger's Test ----------------------------------------------
funnel(m1, main = "Funnel Plot for Publication Bias")
egger_test <- metabias(m1, method = "linreg", k.min = 10)
print(egger_test)

# Subgroup Analysis ---------------------------------------------------------
# Function to Generate Subgroup Forest Plots
generate_subgroup_forest <- function(data, subgroup_var, plot_title, file_name) {
  if (!subgroup_var %in% names(data)) {
    stop(paste("Subgroup variable", subgroup_var, "not found in the dataset."))
  }
  
  data_filtered <- data[!is.na(data[[subgroup_var]]) & data[[subgroup_var]] != "", ]
  
  if (nrow(data_filtered) < 2) {
    warning(paste("Skipping subgroup analysis for", plot_title, 
                  "due to insufficient valid data."))
    return()
  }
  
  m_filtered <- metagen(
    TE = SMD,
    seTE = SE,
    studlab = Study,
    data = data_filtered,
    comb.fixed = TRUE,
    comb.random = TRUE,
    method.tau = "REML",
    hakn = TRUE,
    prediction = TRUE,
    sm = "SMD",
    byvar = data_filtered[[subgroup_var]]  # Define the subgroup
  )
  
  png(filename = paste0(file_name, ".png"), width = 1200, height = 800, res = 150)
  forest(m_filtered, 
         main = plot_title, 
         col.diamond = "blue", 
         col.diamond.lines = "blue")
  dev.off()
  
  cat(paste("Saved forest plot for", plot_title, "as", file_name, ".png\n"))
}

# Perform Subgroup Analysis
generate_subgroup_forest(data, "Intervention.Category", "Intervention Category", "subgroup_intervention")
generate_subgroup_forest(data, "Age.Category", "Age Category", "subgroup_age")
generate_subgroup_forest(data, "BMI.Category", "BMI Category", "subgroup_bmi")
generate_subgroup_forest(data, "Duration.Category", "Duration Category", "subgroup_duration")
generate_subgroup_forest(data, "Adherence.Level", "Adherence Level", "subgroup_adherence")

# Univariable Meta-Regression ------------------------------------------------
meta_reg_intervention <- metareg(m1, ~ Intervention.Category)
meta_reg_age <- metareg(m1, ~ Age.Category)
meta_reg_bmi <- metareg(m1, ~ BMI.Category)
meta_reg_adherence <- metareg(m1, ~ Adherence.Level)
meta_reg_duration <- metareg(m1, ~ Duration.Category)

# Print Univariable Results
print(meta_reg_intervention)
print(meta_reg_age)
print(meta_reg_bmi)
print(meta_reg_adherence)
print(meta_reg_duration)

# Multivariable Meta-Regression
formula_str <- "Intervention.Category + Age.Category + BMI.Category + Adherence.Level + Duration.Category"
meta_reg_multi <- metareg(m1, as.formula(paste("~", formula_str)))
summary(meta_reg_multi)


# Data Cleaning: Remove Missing Core Meta-Analysis Variables
data <- data[!is.na(data$SMD) & !is.na(data$SE), ]

# Meta-Analysis -------------------------------------------------------------
# re-run meta-analysis without missing data
m1 <- metagen(
  TE = SMD,                # Standardized Mean Difference
  seTE = SE,               # Standard Error
  studlab = Study,         # Study Label
  data = data,             # Dataset
  comb.fixed = TRUE,       # Fixed-effect model
  comb.random = TRUE,      # Random-effects model
  method.tau = "REML",     # Restricted Maximum Likelihood
  hakn = TRUE,             # Hartung-Knapp adjustment
  prediction = TRUE,       # Prediction interval
  sm = "SMD"               # Effect size type
)

# Summary of Meta-Analysis
summary(m1)

# Prepare data by removing rows with missing values in critical variables
data <- na.omit(data)

# Define the main meta-analysis object
m1 <- rma(yi = SMD, sei = SE, data = data, method = "REML")

# Interaction Meta-Regression
meta_reg_age_adherence <- rma(yi = SMD, sei = SE, mods = ~ Age.Category * Adherence.Level, data = data)
meta_reg_bmi_adherence <- rma(yi = SMD, sei = SE, mods = ~ BMI.Category * Adherence.Level, data = data)
meta_reg_adherence_duration <- rma(yi = SMD, sei = SE, mods = ~ Adherence.Level * Duration.Category, data = data)
meta_reg_age_bmi <- rma(yi = SMD, sei = SE, mods = ~ Age.Category * BMI.Category, data = data)
meta_reg_age_duration <- rma(yi = SMD, sei = SE, mods = ~ Age.Category * Duration.Category, data = data)
meta_reg_bmi_duration <- rma(yi = SMD, sei = SE, mods = ~ BMI.Category * Duration.Category, data = data)

# Print summaries
summary(meta_reg_age_adherence)
summary(meta_reg_bmi_adherence)
summary(meta_reg_adherence_duration)
summary(meta_reg_age_bmi)
summary(meta_reg_age_duration)
summary(meta_reg_bmi_duration)


###
#install.packages("ggplot2")
#install.packages("patchwork")

# Load necessary libraries
library(ggplot2)
library(patchwork)

create_bubble_plot <- function(meta_reg, title, xlab, ylab) {
  # Extract moderators, coefficients, and standard errors
  data <- data.frame(
    moderator = meta_reg$X[, 2], # First interaction term (assumes interaction is ~ Var1 * Var2)
    coefficient = meta_reg$yi,  # Observed outcomes
    size = 1 / meta_reg$vi      # Weights (inverse variance)
  )
  
  # Create the plot
  ggplot(data, aes(x = moderator, y = coefficient, size = size)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 10),
          legend.position = "none") +
    scale_size_continuous(range = c(2, 8))  # Adjust bubble sizes
}

# Create bubble plots for each interaction
plot1 <- create_bubble_plot(meta_reg_age_adherence,
                            "Age × Adherence Interaction",
                            "Age Category and Adherence Level",
                            "Standardized Mean Difference")
plot2 <- create_bubble_plot(meta_reg_bmi_adherence,
                            "BMI × Adherence Interaction",
                            "BMI Category and Adherence Level",
                            "Standardized Mean Difference")
plot3 <- create_bubble_plot(meta_reg_adherence_duration,
                            "Adherence × Duration Interaction",
                            "Adherence Level and Duration",
                            "Standardized Mean Difference")
plot4 <- create_bubble_plot(meta_reg_age_bmi,
                            "Age × BMI Interaction",
                            "Age Category and BMI Category",
                            "Standardized Mean Difference")
plot5 <- create_bubble_plot(meta_reg_age_duration,
                            "Age × Duration Interaction",
                            "Age Category and Duration",
                            "Standardized Mean Difference")
plot6 <- create_bubble_plot(meta_reg_bmi_duration,
                            "BMI × Duration Interaction",
                            "BMI Category and Duration",
                            "Standardized Mean Difference")

# Combine plots into one figure
combined_plot <- (plot1 + plot2 + plot3) / (plot4 + plot5 + plot6)

# Save the combined plot
ggsave("combined_bubble_plots.png", combined_plot, width = 12, height = 8, dpi = 300)

#### end of script ###

