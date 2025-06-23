# Survival Analysis Simulation Exercise
# Exponential Survival Data with Censoring

# Load required libraries
library(survival)
library(survminer)

# Set seed for reproducibility
set.seed(123)

# Step 1: Define Simulation Parameters
n <- 500                    # Total patients
prop_treatment <- 0.5       # Proportion in treatment arm
lambda_control <- 0.1       # Hazard rate for control
lambda_treatment <- 0.07    # Hazard rate for treatment
lambda_censor <- 0.05       # Censoring rate

# Randomly assign patients to groups
n_treatment <- sum(rbinom(n, 1, prop_treatment))
n_control <- n - n_treatment

cat("Patients in Control:", n_control, "\n")
cat("Patients in Treatment:", n_treatment, "\n\n")

# Step 2: Generate Survival and Censoring Times

# Control group
true_survival_control <- rexp(n_control, rate = lambda_control)
censoring_time_control <- rexp(n_control, rate = lambda_censor)
observed_time_control <- pmin(true_survival_control, censoring_time_control)
status_control <- as.integer(true_survival_control <= censoring_time_control)

# Treatment group
true_survival_treatment <- rexp(n_treatment, rate = lambda_treatment)
censoring_time_treatment <- rexp(n_treatment, rate = lambda_censor)
observed_time_treatment <- pmin(true_survival_treatment, censoring_time_treatment)
status_treatment <- as.integer(true_survival_treatment <= censoring_time_treatment)

# Step 3: Combine the Data
control_data <- data.frame(
  id = 1:n_control,
  group = "Control",
  time = observed_time_control,
  status = status_control,
  true_survival = true_survival_control,
  censor_time = censoring_time_control
)

treatment_data <- data.frame(
  id = (n_control + 1):n,
  group = "Treatment",
  time = observed_time_treatment,
  status = status_treatment,
  true_survival = true_survival_treatment,
  censor_time = censoring_time_treatment
)

# Combine datasets
surv_data <- rbind(control_data, treatment_data)
surv_data$group <- factor(surv_data$group, levels = c("Control", "Treatment"))

# Display summary statistics
cat("Dataset Summary:\n")
cat("Total observations:", nrow(surv_data), "\n")
cat("Events:", sum(surv_data$status), "\n")
cat("Censored:", sum(1 - surv_data$status), "\n")
cat("Censoring percentage:", round(mean(1 - surv_data$status) * 100, 1), "%\n\n")

# Summary by group
cat("Summary by Group:\n")
aggregate(status ~ group, data = surv_data, 
          FUN = function(x) c(n = length(x), 
                             events = sum(x), 
                             censored = sum(1-x),
                             event_rate = round(mean(x)*100, 1)))

# Step 4: Kaplan-Meier Curves

# Create survival object
surv_obj <- Surv(surv_data$time, surv_data$status)

# Fit Kaplan-Meier curves
km_fit <- survfit(surv_obj ~ group, data = surv_data)

# Print summary
print(km_fit)

# Base R plot
par(mfrow = c(1, 2))
plot(km_fit, 
     col = c("red", "blue"), 
     lty = 1:2,
     xlab = "Time", 
     ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curves (Base R)")
legend("topright", 
       legend = c("Control", "Treatment"), 
       col = c("red", "blue"), 
       lty = 1:2)

# ggsurvplot (if survminer is available)
gg_plot <- ggsurvplot(km_fit, 
                      data = surv_data,
                      pval = TRUE,
                      conf.int = TRUE,
                      risk.table = TRUE,
                      risk.table.col = "strata",
                      legend.labs = c("Control", "Treatment"),
                      palette = c("red", "blue"),
                      title = "Kaplan-Meier Survival Curves",
                      xlab = "Time",
                      ylab = "Survival Probability")
print(gg_plot)

# Step 5: Cox Proportional Hazards Model

# Fit Cox model
cox_model <- coxph(Surv(time, status) ~ group, data = surv_data)

# Display results
cat("\n\nCox Proportional Hazards Model Results:\n")
summary(cox_model)

# Extract hazard ratio and confidence interval
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))

cat("\n\nHazard Ratio (Treatment vs Control):", round(hr, 3), "\n")
cat("95% CI: [", round(ci[1], 3), ",", round(ci[2], 3), "]\n")
cat("p-value:", round(summary(cox_model)$coefficients[5], 4), "\n")

# Bonus Challenge: Add Age Covariate

# Generate age covariate (normally distributed, mean 60, sd 10)
surv_data$age <- rnorm(n, mean = 60, sd = 10)

# Center age for better interpretation
surv_data$age_centered <- surv_data$age - mean(surv_data$age)

# Fit Cox model with age
cox_model_age <- coxph(Surv(time, status) ~ group + age_centered, data = surv_data)

cat("\n\nCox Model with Age Covariate:\n")
summary(cox_model_age)

# Compare models
cat("\n\nModel Comparison:\n")
cat("Without age - HR:", round(exp(coef(cox_model)), 3), "\n")
cat("With age - HR:", round(exp(coef(cox_model_age))[1], 3), "\n")

# Visualize age distribution by group
par(mfrow = c(1, 2))
boxplot(age ~ group, data = surv_data, 
        main = "Age Distribution by Group",
        xlab = "Group", ylab = "Age")

# Check proportional hazards assumption
ph_test <- cox.zph(cox_model)
plot(ph_test, main = "Schoenfeld Residuals")

# Additional Analysis: Median Survival Times
median_surv <- survfit(Surv(time, status) ~ group, data = surv_data)
cat("\n\nMedian Survival Times:\n")
print(median_surv)

# Create a summary table
summary_table <- data.frame(
  Group = c("Control", "Treatment"),
  N = c(n_control, n_treatment),
  Events = aggregate(status ~ group, data = surv_data, sum)$status,
  Median_Survival = summary(median_surv)$table[, "median"],
  True_Median = c(log(2)/lambda_control, log(2)/lambda_treatment)
)

cat("\n\nSummary Table:\n")
print(summary_table)

# Verify theoretical vs observed
cat("\n\nTheoretical vs Observed Median Survival:\n")
cat("Control - Theoretical:", round(log(2)/lambda_control, 2), 
    "Observed:", round(summary_table$Median_Survival[1], 2), "\n")
cat("Treatment - Theoretical:", round(log(2)/lambda_treatment, 2), 
    "Observed:", round(summary_table$Median_Survival[2], 2), "\n")