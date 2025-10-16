# These scripts are adapted from the original code used in the study and have been modified for public release.
# They may depend on specific aspects of the developer’s computing environment; users should verify compatibility before running.
# For inquiries regarding the code or its usage, please contact Yoshito Saito (yoshitos@student.unimelb.edu.au).

# Group Comparison
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggseg)
library(ggrepel)
library(ppcor)

data_start= # column number where data you are going to check group difference starts
data_end= # column number where data you are going to check group difference ends
n=78 # n of HC and EOP
group1_name=0 # HC
group2_name=1 # EOP

# results data
results <- data.frame(
    column_name = character(),
    p_value = numeric(),
    z_statistic = numeric(),
    mean1 = numeric(),
    sd1 = numeric(),
    mean1 = numeric(),
    sd1 = numeric(),
    effect_size_d = numeric(),
    group_difference = numeric(),
    Higher_mean_group = character()
)

data <- read_excel("SWM_results_with_demographics.xlsx")
data <- data[data$dx_cat == group1_name | data$dx_cat == group2_name,]

get_higher_mean_group <- function(data_col, grouping_col) {
    means <- tapply(data_col, grouping_col, mean, na.rm = TRUE)
    return(names(means)[which.max(means)])
}

# Function to calculate Cohen's d
calculate_cohens_d <- function(group1, group2) {
    n1 <- length(group1)
    n2 <- length(group2)
    mean1 <- mean(group1, na.rm = TRUE)
    mean2 <- mean(group2, na.rm = TRUE)
    sd1 <- sd(group1, na.rm = TRUE)
    sd2 <- sd(group2, na.rm = TRUE)
    pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    d <- (mean2 - mean1) / pooled_sd
    group_difference <- ((mean2 - mean1)/ mean1) * 100
    return(list(mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2, cohens_d = d, group_difference = group_difference))
}

for(i in data_start:data_end) {
    target_col_name <- names(data)[i]
    group1_data <- data[data$dx_cat == group1_name, target_col_name, drop = TRUE]
    group2_data <- data[data$dx_cat == group2_name, target_col_name, drop = TRUE]
    n1 <- length(na.omit(group1_data))
    n2 <- length(na.omit(group2_data))
    
    if (n1 <= (n-1) | n2 <= (n-1)) {
        # If there are missing values, the analysis will not be performed
        results <- rbind(results, data.frame(
            column_name = target_col_name,
            p_value = NA,
            z_statistic = NA,
            mean1 = NA,
            sd1 = NA,
            mean2 = NA,
            sd2 = NA,
            effect_size_d = NA, 
            group_difference = NA,
            Higher_mean_group = NA
        ))
    } else {
        if (grepl("_FA$", target_col_name)) {
            model <- lm(formula = data[[target_col_name]] ~ dx_cat + interview_age + sex + Abs_mov + Rel_mov, data = data)
            summary_model <- summary(model)
        } else if (grepl("_fd$", target_col_name)) {
            model <- lm(formula = data[[target_col_name]] ~ dx_cat + interview_age + sex + Abs_mov + Rel_mov, data = data, na.action = na.omit)
            summary_model <- summary(model) # + eTIV
        } else {
            # Linear model with dx_cat, interview_age, and sex as predictors
            model <- lm(formula = data[[target_col_name]] ~ dx_cat + interview_age + sex, data = data)
            summary_model <- summary(model)
        }
        
        # Extract the coefficient, statistic, and p-value for dx_cat
        dx_cat_coef <- summary_model$coefficients['dx_cat', ]
        dx_cat_p_value <- dx_cat_coef["Pr(>|t|)"]
        dx_cat_statistic <- dx_cat_coef["t value"]
        
        # Calculate Cohen's d for effect size between group 1 and group 3
        effect_size_results <- calculate_cohens_d(group1_data, group2_data)
        mean1 <- effect_size_results$mean1
        mean2 <- effect_size_results$mean2
        sd1 <- effect_size_results$sd1
        sd2 <- effect_size_results$sd2
        cohens_d <- effect_size_results$cohens_d
        group_difference <- effect_size_results$group_difference
        
        # Determine the group with the higher mean
        higher_mean_group <- get_higher_mean_group(data[[target_col_name]], data$dx_cat)
        
        results <- rbind(results, data.frame(
                column_name = target_col_name,
                p_value = dx_cat_p_value,
                z_statistic = dx_cat_statistic,
                mean1 = mean1,
                sd1 = sd1,
                mean2 = mean2,
                sd2 = sd2,
                effect_size_d = cohens_d, 
                group_difference = group_difference,
                Higher_mean_group = higher_mean_group
        ))
    }
}

write.csv(results, "groupcomparison_results.csv", row.names = FALSE)

#BH correction
suffixes <- c("_GMWM_fd", "_L1_fd", "_GMWM_FA", "_L1_FA", "_area", "_thickness")
p_columns <- c("p_value")

results <- read.csv("groupcomparison_results.csv", row.names = NULL) 
for(col in p_columns) {
    results[paste0(col, "_BH_corrected")] <- NA
}
for(suffix in suffixes) {
  target_rows <- grep(paste0(suffix, "$"), results$column_name)
    if(length(target_rows) > 0) {
    for(col in p_columns) {
      original_p_values <- results[target_rows, col]
      corrected_p_values <- p.adjust(original_p_values, method = "BH")
      results[target_rows, paste0(col, "_BH_corrected")] <- corrected_p_values
    }
  }
}

write.csv(results, "groupcomparison_results_BHcorrect.csv", row.names = FALSE)


# Correlation between thickness and FA/FD based on effect size d
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggseg)
library(ggrepel)

data <- read.csv("groupcomparison_results_BHcorrect.csv")
region_mapping <- read.csv("region_mapping_des.csv") # It contains two columns, “label” and “region,” which list the name of each ROI and its corresponding brain region classification (e.g., Frontal, Temporal, etc.).
data <- data %>%
  mutate(regionname = gsub("_area$|_thickness$|_GMWM_fd$|_L1_fd$|_GMWM_FA$|_L1_FA$", "", column_name))
data <- data %>%
  left_join(region_mapping, by = c("regionname" = "label"))

thickness_data <- data %>% filter(grepl("_thickness$", column_name))
GMWM_fd_data <- data %>% filter(grepl("_GMWM_fd$", column_name))
L1_fd_data <- data %>% filter(grepl("_L1_fd$", column_name))
GMWM_FA_data <- data %>% filter(grepl("_GMWM_FA$", column_name))
L1_FA_data <- data %>% filter(grepl("_L1_FA$", column_name))

calculate_correlation <- function(data1, data2, region) {
  subset1 <- data1 %>% filter(region == !!region)
  subset2 <- data2 %>% filter(region == !!region)
  correlation <- cor.test(subset1$effect_size_d, subset2$effect_size_d, method = "pearson")
  return(list(r_value = correlation$estimate, p_value = correlation$p.value))
}

comparison_datasets <- list(
  "GMWM_fd" = GMWM_fd_data,
  "L1_fd" = L1_fd_data,
  "GMWM_FA" = GMWM_FA_data,
  "L1_FA" = L1_FA_data
)

unique_regions <- unique(thickness_data$region)
results <- list()

for (dataset_name in names(comparison_datasets)) {
  subset2_data <- comparison_datasets[[dataset_name]]
  
  for (region in unique_regions) {
    subset1 <- thickness_data %>% filter(region == !!region)
    subset2 <- subset2_data %>% filter(region == !!region)
    if (nrow(subset1) > 2 & nrow(subset2) > 2) {
      result <- calculate_correlation(subset1, subset2, region)
      results[[paste0(region, "_", dataset_name)]] <- result
      cat("Region:", region, " | Dataset:", dataset_name, "\n")
      cat("r =", result$r_value, ", p =", result$p_value, "\n\n")
    } else {
      results[[paste0(region, "_", dataset_name)]] <- list(r_value = NA, p_value = NA)
      cat("Region:", region, " | Dataset:", dataset_name, "\n")
      cat("Not enough data to calculate correlation.\n\n")
    }
  }
}  

plot_scatter <- function(x_data, y_data, y_label) {
    correlation <- cor.test(x_data$effect_size_d, y_data$effect_size_d, method = "pearson")
    p_value <- correlation$p.value
    r_value <- correlation$estimate
    r_squared <- r_value^2 
    
    if (p_value > 0.00049){
        ggplot() +
        geom_point(data = x_data, aes(x = effect_size_d, y = y_data$effect_size_d, color = region), size = 3) +
        geom_smooth(aes(x = x_data$effect_size_d, y = y_data$effect_size_d), method = "lm", se = TRUE, color = "black", size = 1.5) +
        scale_color_manual(values = c(
            "Frontal" = "red2",
            "Temporal" = "blue1",
            "Parietal" = "green4",
            "Occipital" = "purple2",
            "Cingulate" = "orange1",
            "Insula" = "brown"
        )) +
        labs(x = "Effect size: CT", y = paste("Effect size: ", y_label, sep = "")) +
        ggtitle(paste("p =", round(p_value, 3), ", r =", round(r_value, 3), ", R² =", round(r_squared, 3))) +
        theme(
            axis.title.x = element_text(size = 20), 
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.position = "bottom",
            panel.background = element_blank(), 
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", size = 0.8)
        )
    } else {
                ggplot() +
        geom_point(data = x_data, aes(x = effect_size_d, y = y_data$effect_size_d, color = region), size = 3) +
        geom_smooth(aes(x = x_data$effect_size_d, y = y_data$effect_size_d), method = "lm", se = TRUE, color = "black", size = 1.5) +
        scale_color_manual(values = c(
            "Frontal" = "red2",
            "Temporal" = "blue1",
            "Parietal" = "green4",
            "Occipital" = "purple2",
            "Cingulate" = "orange1",
            "Insula" = "brown"
        )) +
        labs(x = "Effect size: CT", y = paste("Effect size: ", y_label, sep = "")) +
        ggtitle(paste("p < 0.0005", ", r =", round(r_value, 3), ", R² =", round(r_squared, 3))) +
        theme(
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.position = "bottom",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", size = 0.8)
        )
    }
}

plot_thickness_vs_GMWM_fd <- plot_scatter(thickness_data, GMWM_fd_data, "FD") 
ggsave(filename = "Correlation_FD_CT_GMWM.png", plot = plot_thickness_vs_GMWM_fd, dpi = 300, width = 7, height = 8)
plot_thickness_vs_L1_fd <- plot_scatter(thickness_data, L1_fd_data, "FD")
ggsave(filename = "Correlation_FD_CT_L1.png", plot = plot_thickness_vs_L1_fd, dpi = 300, width = 7, height = 8)
plot_thickness_vs_GMWM_FA <- plot_scatter(thickness_data, GMWM_FA_data, "FA") 
ggsave(filename = "Correlation_FA_CT_GMWM.png", plot = plot_thickness_vs_GMWM_FA, dpi = 300, width = 7, height = 8)
plot_thickness_vs_L1_FA <- plot_scatter(thickness_data, L1_FA_data, "FA")
ggsave(filename = "Correlation_FA_CT_L1.png", plot = plot_thickness_vs_L1_FA, dpi = 300, width = 7, height = 8)

# Correlation between surface area and FA/FD
area_data <- data %>% filter(grepl("_area$", column_name))

comparison_datasets <- list(
    "GMWM_fd_sr" = GMWM_fd_data,
    "L1_fd_sr" = L1_fd_data,
    "GMWM_FA_sr" = GMWM_FA_data,
    "L1_FA_sr" = L1_FA_data
)

unique_regions <- unique(area_data$region)
results <- list()

for (dataset_name in names(comparison_datasets)) {
    subset2_data <- comparison_datasets[[dataset_name]]
    for (region in unique_regions) {
        subset1 <- area_data %>% filter(region == !!region)
        subset2 <- subset2_data %>% filter(region == !!region)
        if (nrow(subset1) > 2 & nrow(subset2) > 2) {
            result <- calculate_correlation(subset1, subset2, region)
            results[[paste0(region, "_", dataset_name)]] <- result
            cat("Region:", region, " | Dataset:", dataset_name, "\n")
            cat("r =", result$r_value, ", p =", result$p_value, "\n\n")
        } else {
            results[[paste0(region, "_", dataset_name)]] <- list(r_value = NA, p_value = NA)
            cat("Region:", region, " | Dataset:", dataset_name, "\n")
            cat("Not enough data to calculate correlation.\n\n")
        }
    }
}  

plot_scatter <- function(x_data, y_data, y_label) {
    correlation <- cor.test(x_data$effect_size_d, y_data$effect_size_d, method = "pearson")
    p_value <- correlation$p.value
    r_value <- correlation$estimate
    r_squared <- r_value^2
    if (p_value > 0.00049){
        ggplot() +
            geom_point(data = x_data, aes(x = effect_size_d, y = y_data$effect_size_d, color = region), size = 3) +
            geom_smooth(aes(x = x_data$effect_size_d, y = y_data$effect_size_d), method = "lm", se = TRUE, color = "black", size = 1.5) +
            scale_color_manual(values = c(
                "Frontal" = "red2",
                "Temporal" = "blue1",
                "Parietal" = "green4",
                "Occipital" = "purple2",
                "Cingulate" = "orange1",
                "Insula" = "brown"
            )) +
            labs(x = "Effect size: SA", y = paste("Effect size: ", y_label, sep = "")) +
            ggtitle(paste("p =", round(p_value, 3), ", r =", round(r_value, 3), ", R² =", round(r_squared, 3))) +
            theme(
                axis.title.x = element_text(size = 20), 
                axis.title.y = element_text(size = 20), 
                axis.text.x = element_text(size = 15), 
                axis.text.y = element_text(size = 15), 
                plot.title = element_text(size = 20, hjust = 0.5),
                legend.position = "bottom",
                panel.background = element_blank(), 
                panel.grid = element_blank(),
                axis.line = element_line(color = "black", size = 0.8) 
            )
    } else {
        ggplot() +
            geom_point(data = x_data, aes(x = effect_size_d, y = y_data$effect_size_d, color = region), size = 3) +
            geom_smooth(aes(x = x_data$effect_size_d, y = y_data$effect_size_d), method = "lm", se = TRUE, color = "black", size = 1.5) +
            scale_color_manual(values = c(
                "Frontal" = "red2",
                "Temporal" = "blue1",
                "Parietal" = "green4",
                "Occipital" = "purple2",
                "Cingulate" = "orange1",
                "Insula" = "brown"
            )) +
            labs(x = "Effect size: SA", y = paste("Effect size: ", y_label, sep = "")) +
            ggtitle(paste("p < 0.0005", ", r =", round(r_value, 3), ", R² =", round(r_squared, 3))) +
            theme(
                axis.title.x = element_text(size = 20), 
                axis.title.y = element_text(size = 20), 
                axis.text.x = element_text(size = 15), 
                axis.text.y = element_text(size = 15),
                plot.title = element_text(size = 20, hjust = 0.5),
                legend.position = "bottom",
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_line(color = "black", size = 0.8)
            )
    }
}

plot_area_vs_GMWM_fd <- plot_scatter(area_data, GMWM_fd_data, "FD") 
ggsave(filename = "Correlation_FD_SA_GMWM.png", plot = plot_area_vs_GMWM_fd, dpi = 300, width = 7, height = 8)
plot_area_vs_L1_fd <- plot_scatter(area_data, L1_fd_data, "FD")
ggsave(filename = "Correlation_FD_SA_L1.png", plot = plot_area_vs_L1_fd, dpi = 300, width = 7, height = 8)
plot_area_vs_GMWM_FA <- plot_scatter(area_data, GMWM_FA_data, "FA") 
ggsave(filename = "Correlation_FA_SA_GMWM.png", plot = plot_area_vs_GMWM_FA, dpi = 300, width = 7, height = 8)
plot_area_vs_L1_FA <- plot_scatter(area_data, L1_FA_data, "FA")
ggsave(filename = "Correlation_FA_SA_L1.png", plot = plot_area_vs_L1_FA, dpi = 300, width = 7, height = 8)


# Polynominal regression on thickness and FA/FD
library(readxl)
library(dplyr)
library(tidyr)
library(viridis)
library(tibble)
library(patchwork)
library(tidyverse)
library(ggsegDesterieux)
library(ggseg)
library(ggrepel)
library(ggplot2)

data <- read_excel("SWM_results_with_demographics.xlsx")
data$interview_age <- data$interview_age / 12
pattern_list <- c("_area", "_thickness$", "_GMWM_fd$", "_L1_fd$", "_GMWM_FA$", "_L1_FA$")
results_WM <- data.frame()
results_GM <- data.frame()
data_start= # column number where data you are going to check group difference starts
data_end= # column number where data you are going to check group difference ends

for (i in data_start:data_end) {
    target_col_name <- names(data)[i]
    group1_data <- data[data$dx_cat == 0, target_col_name, drop = TRUE]
    group2_data <- data[data$dx_cat == 1, target_col_name, drop = TRUE]
    n1 <- length(na.omit(group1_data))
    n2 <- length(na.omit(group2_data))
    sd1 <- sd(group1_data, na.rm = TRUE)
    sd2 <- sd(group2_data, na.rm = TRUE)
    if (any(sapply(pattern_list, function(p) grepl(p, target_col_name)))) {
        for (A in 16:30) {        
            if (grepl("_fd", target_col_name) || grepl("_FA", target_col_name)) {
                polyformula <- as.formula(paste("`", target_col_name, "` ~ dx_cat + I(interview_age -", A, ") + dx_cat:I(interview_age -", A, ") + sex + Abs_mov + Rel_mov", sep = ""))
                polymodel <- lm(polyformula, data = data)
                summary_model <- summary(polymodel)

                beta_values_WM <- summary_model$coefficients[, "Estimate"]
                se_values_WM <- summary_model$coefficients[, "Std. Error"]
                p_values_WM <- summary_model$coefficients[, "Pr(>|t|)"]
                sigma_values_WM <- summary_model$sigma
                t_values_WM <- beta_values_WM / se_values_WM
                d_values_WM <- beta_values_WM / sigma_values_WM
                
                results_WM <- rbind(results_WM, c(target_col_name, A, beta_values_WM, se_values_WM, p_values_WM, t_values_WM, d_values_WM))
            } else {
                polyformula <- as.formula(paste("`", target_col_name, "` ~ dx_cat + I(interview_age -", A, ") + dx_cat:I(interview_age -", A, ") + sex", sep = ""))
                polymodel <- lm(polyformula, data = data)
                summary_model <- summary(polymodel)
                
                beta_values_GM <- summary_model$coefficients[, "Estimate"]
                se_values_GM <- summary_model$coefficients[, "Std. Error"]
                p_values_GM <- summary_model$coefficients[, "Pr(>|t|)"]
                sigma_values_GM <- summary_model$sigma
                t_values_GM <- beta_values_GM / se_values_GM
                d_values_GM <- beta_values_GM / sigma_values_GM
                
                results_GM <- rbind(results_GM, c(target_col_name, A, beta_values_GM, se_values_GM, p_values_GM, t_values_GM, d_values_GM))
            }
        }
    }
}

colnames(results_WM) <- c("Variable", "Age", names(beta_values_WM), paste0("SE_", names(beta_values_WM)), paste0("p_", names(beta_values_WM)), paste0("t_", names(beta_values_WM)), paste0("d_", names(beta_values_WM)))
colnames(results_GM) <- c("Variable", "Age", names(beta_values_GM), paste0("SE_", names(beta_values_GM)), paste0("p_", names(beta_values_GM)), paste0("t_", names(beta_values_GM)), paste0("d_", names(beta_values_GM)))

write.csv(results_WM, "polynomial_result_fdFA.csv", row.names = FALSE)
write.csv(results_GM, "polynomial_result_CTSA.csv", row.names = FALSE)

FD_filtered <- subset(results_WM, grepl("GMWM_fd$", Variable))
d_dx_cat_FD <- FD_filtered[, c("Variable", "Age", "dx_cat", "p_dx_cat", "SE_dx_cat", "d_dx_cat")]
FA_filtered <- subset(results_WM, grepl("GMWM_FA$", Variable))
d_dx_cat_FA <- FA_filtered[, c("Variable", "Age", "dx_cat", "p_dx_cat", "SE_dx_cat", "d_dx_cat")]
CT_filtered <- subset(results_GM, grepl("thickness$", Variable))
d_dx_cat_CT <- CT_filtered[, c("Variable", "Age", "dx_cat", "p_dx_cat", "SE_dx_cat", "d_dx_cat")]
SA_filtered <- subset(results_GM, grepl("thickness$", Variable))
d_dx_cat_SA <- SA_filtered[, c("Variable", "Age", "dx_cat", "p_dx_cat", "SE_dx_cat", "d_dx_cat")]

d_dx_cat_FD$d_dx_cat <- as.numeric(as.character(d_dx_cat_FD$d_dx_cat))
d_dx_cat_FA$d_dx_cat <- as.numeric(as.character(d_dx_cat_FA$d_dx_cat))
d_dx_cat_CT$d_dx_cat <- as.numeric(as.character(d_dx_cat_CT$d_dx_cat))
d_dx_cat_SA$d_dx_cat <- as.numeric(as.character(d_dx_cat_SA$d_dx_cat))

# Visualization
data_list <- list(
    D_fd = d_dx_cat_FD,
    D_FA = d_dx_cat_FA, 
    D_thickness = d_dx_cat_CT,
    D_area = d_dx_cat_SA
)

create_someData <- function(data, suffix) {
    tibble(
        label = gsub(suffix, "", data$column_name),
        d = c(data$d_dx_cat),
        p_value = c(data$adjusted_p_dx_cat)
    )
}

for (data_name in names(data_list)) {
    data_subset <- data_list[[data_name]]
    suffix <- ifelse(data_name == "D_fd", "_GMWM_fd",
                     ifelse(data_name == "D_FA", "_GMWM_FA",
                            ifelse(data_name == "D_thickness", "_thickness", "_area")))
    all_plots <- list()
    plots_forfig <- list()
    for (age in unique(data_subset$Age)) {
        age_data <- subset(data_subset, Age == age)
        age_data$column_name <- age_data$Variable
        data_original <- age_data %>% filter(grepl(suffix, column_name))
        data_extracted <- data_original %>%
            mutate(across(starts_with("p_"), ~ p.adjust(., method = "BH"), .names = "adjusted_{.col}"))
        #data_extracted$d_dx_cat[data_extracted$adjusted_p_dx_cat > 0.05] <- NA
        someData <- create_someData(data_extracted, suffix)

        plot <- ggplot(someData) +
          geom_brain(
            atlas = desterieux,
            position = position_brain(hemi ~ side),
            aes(fill = d),
          ) +
          scale_fill_gradient2(
            low = "red4", mid = "white", high = "green4",
            midpoint = 0, limits = c(-1, 1), na.value = "grey80",
            oob = scales::squish
          ) +
          theme_void() +
          labs(title = paste("Effect size mapping for", data_name, "Age", age)) +
          theme(plot.title = element_text(size = 10))

        if (age == 16 | age == 20 | age== 25 | age ==30 ){
            plots_forfig[[paste0("plot_", age)]] <- plot
        }
        all_plots[[paste0("plot_", age)]] <- plot
    }
    combined_plot <- wrap_plots(all_plots, ncol = length(unique(data_subset$Age)))
    combined_plot_forfig <- wrap_plots(plots_forfig, ncol = length(unique(data_subset$Age)))

    ggsave(filename = paste0("combined_brain_map_", data_name, "_GMWM_cohend.png"), plot = combined_plot, width = 49, height = 5, dpi = 300)
    ggsave(filename = paste0("combined_brain_map_", data_name, "_GMWM_cohend_forfig.png"), plot = combined_plot_forfig, width = 49, height = 5, dpi = 300)
}

# Correlation between CT and FA/FD by age based on polynomial regression model
scatter_plot_correlation <- function(data_WM, data_GM, age, WM_name, GM_name) {
    merged_data <- merge(data_GM, data_WM, by = "label", suffixes = c(paste0("_", GM_name), paste0("_", WM_name))) 
    merged_data[[paste0("d_dx_cat_", GM_name)]] <- as.numeric(merged_data[[paste0("d_dx_cat_", GM_name)]])
    merged_data[[paste0("d_dx_cat_", WM_name)]] <- as.numeric(merged_data[[paste0("d_dx_cat_", WM_name)]])
    
    cor_test <- cor.test(merged_data[[paste0("d_dx_cat_", GM_name)]], merged_data[[paste0("d_dx_cat_", WM_name)]])
    cor_value <- round(cor_test$estimate, 2)
    p_value <- round(cor_test$p.value, 3)

    ggplot(merged_data, aes(x = get(paste0("d_dx_cat_", GM_name)), y = get(paste0("d_dx_cat_", WM_name)))) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", color = "red", size = 1.5) +
        annotate("text", x = Inf, y = Inf, label = paste("r =", cor_value, "\np =", p_value), 
                 hjust = 1.1, vjust = 2, size = 5, color = "black") +
        labs(x = paste("Effect size: ", GM_name, "(Age", age, ")"),
             y = paste("Effect size: ", WM_name, "(Age", age, ")"),
             title = paste("Correlation between", GM_name, "and", WM_name, "at Age", age)) +
        theme_minimal()
}

scatter_plots_by_group_CT <- list(
    D_fd = list(),
    D_FA = list()
)
scatter_plots_by_group_SA <- list(
    D_fd = list(),
    D_FA = list()
)
correlation_results_CT <- data.frame(Age = numeric(), other_data_name = character(), corr_value = numeric(), p_value = numeric())
correlation_results_SA <- data.frame(Age = numeric(), other_data_name = character(), corr_value = numeric(), p_value = numeric())


for (age in unique(d_dx_cat_CT$Age)) {
    for (WM_data_name in c("FD", "FA")) {
        data_thickness <- subset(d_dx_cat_CT, Age == age)
        data_thickness$label <- gsub("_thickness", "", data_thickness$Variable)
        data_area <- subset(d_dx_cat_SA, Age == age)
        data_area$label <- gsub("_area", "", data_area$Variable)
        data_WM <- get(paste0("d_dx_cat_", WM_data_name))
        data_WM <- subset(data_WM, Age == age)
        data_WM$label <- gsub(ifelse(WM_data_name == "FD", "_GMWM_fd", "_GMWM_FA"), "", data_WM$Variable)
        
        scatter_plot <- scatter_plot_correlation(data_WM, data_thickness, age, WM_data_name, "thickness")
        if (!is.null(scatter_plot)) {
            scatter_plots_by_group_CT[[WM_data_name]][[paste0(WM_data_name, "_scatter_", age)]] <- scatter_plot
        }
        scatter_plot <- scatter_plot_correlation(data_WM, data_area, age, WM_data_name, "area")
        if (!is.null(scatter_plot)) {
            scatter_plots_by_group_SA[[WM_data_name]][[paste0(WM_data_name, "_scatter_", age)]] <- scatter_plot
        }

        merged_data <- merge(data_thickness, data_WM, by = "label", suffixes = c("_thickness", paste0("_", WM_data_name)))
        # Ensure numeric values for correlation
        merged_data$d_dx_cat_thickness <- as.numeric(merged_data$d_dx_cat_thickness)
        merged_data[[paste0("d_dx_cat_", WM_data_name)]] <- as.numeric(merged_data[[paste0("d_dx_cat_", WM_data_name)]])

        if (nrow(merged_data) > 1) { # Ensure enough data for correlation
            cor_test <- cor.test(merged_data$d_dx_cat_thickness, merged_data[[paste0("d_dx_cat_", WM_data_name)]])

            # Store the results
            correlation_results_CT <- rbind(correlation_results_CT, data.frame(
                Age = age,
                other_data_name = WM_data_name,
                corr_value = cor_test$estimate,
                p_value = cor_test$p.value
            ))
        }

        merged_data <- merge(data_area, data_WM, by = "label", suffixes = c("_area", paste0("_", WM_data_name)))
        # Ensure numeric values for correlation
        merged_data$d_dx_cat_area <- as.numeric(merged_data$d_dx_cat_area)
        merged_data[[paste0("d_dx_cat_", WM_data_name)]] <- as.numeric(merged_data[[paste0("d_dx_cat_", WM_data_name)]])

        if (nrow(merged_data) > 1) { # Ensure enough data for correlation
            cor_test <- cor.test(merged_data$d_dx_cat_area, merged_data[[paste0("d_dx_cat_", WM_data_name)]])

            # Store the results
            correlation_results_SA <- rbind(correlation_results_SA, data.frame(
                Age = age,
                other_data_name = WM_data_name,
                corr_value = cor_test$estimate,
                p_value = cor_test$p.value
            ))
        }
    }
}

for (other_data_name in names(scatter_plots_by_group_CT)) {
    all_scatter_plots <- scatter_plots_by_group_CT[[other_data_name]]
    if (length(all_scatter_plots) > 0) {
        combined_scatter_plot <- wrap_plots(all_scatter_plots, ncol = length(unique(d_dx_cat_CT$Age)))
        ggsave(filename = paste0("combined_scatter_plots_CT_", other_data_name, "_GMWM.png"), plot = combined_scatter_plot, width = 49, height = 24, dpi = 300)
    }
}

for (other_data_name in names(scatter_plots_by_group_SA)) {
    all_scatter_plots <- scatter_plots_by_group_SA[[other_data_name]]
    if (length(all_scatter_plots) > 0) {
        combined_scatter_plot <- wrap_plots(all_scatter_plots, ncol = length(unique(d_dx_cat_SA$Age)))
        ggsave(filename = paste0("combined_scatter_plots_SA_", other_data_name, "_GMWM.png"), plot = combined_scatter_plot, width = 49, height = 24, dpi = 300)
    }
}

for (current_other_data_name in unique(correlation_results_CT$other_data_name)) {
    data_for_plot <- subset(correlation_results_CT, other_data_name == current_other_data_name)
    p <- ggplot() +
        geom_point(
            data = data_for_plot,
            aes(
                x = as.numeric(Age),
                y = corr_value,
                fill = p_value <= 0.05,
                alpha = ifelse(p_value <= 0.05, 1, 1),
                size = ifelse(p_value <= 0.05, 5, 5)
            ),
            shape = 21,
            color = "black",
            stroke = 0.7
        ) +
        scale_fill_manual(
            values = c(
                "FALSE" = "white",
                "TRUE" = "red"
            ),
            name = "Significance"
        ) +
        labs(
            x = "Age",
            y = "Correlation",
            title = paste("Correlation vs Age for", current_other_data_name),
            alpha = "Transparency",
            size = "Significance"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_alpha_identity() +
        scale_size_identity() +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            axis.ticks.length = unit(0.25, "lines"),
            axis.line = element_line(color = "black"),
            legend.position = "right"
        )
    ggsave(
        filename = paste0("correlation_vs_age_CT_", current_other_data_name, "_wholebrain_circle_alpha.png"),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
    )
}

for (current_other_data_name in unique(correlation_results_SA$other_data_name)) {
    data_for_plot <- subset(correlation_results_SA, other_data_name == current_other_data_name)
    p <- ggplot() +
        geom_point(
            data = data_for_plot,
            aes(
                x = as.numeric(Age),
                y = corr_value,
                fill = p_value <= 0.05,
                alpha = ifelse(p_value <= 0.05, 1, 1),
                size = ifelse(p_value <= 0.05, 5, 5)
            ),
            shape = 21,
            color = "black",
            stroke = 0.7
        ) +
        scale_fill_manual(
            values = c(
                "FALSE" = "white",
                "TRUE" = "red"
            ),
            name = "Significance"
        ) +
        labs(
            x = "Age",
            y = "Correlation",
            title = paste("Correlation vs Age for", current_other_data_name),
            alpha = "Transparency",
            size = "Significance"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        scale_alpha_identity() +
        scale_size_identity() +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            axis.ticks.length = unit(0.25, "lines"),
            axis.line = element_line(color = "black"),
            legend.position = "right"
        )
    ggsave(
        filename = paste0("correlation_vs_age_SA_", current_other_data_name, "_wholebrain_circle_alpha.png"),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
    )
}


# Association with cognition / mediation analysis
library(tidyverse)
library(broom)
library(readxl)
library(lavaan)
library(mediation)
library(dplyr)

data <- read_excel("SWM_results_with_demographics.xlsx")
data <- data %>% filter(dx_cat %in% c(0,1))
data <- data %>% mutate(
    dx_cat = factor(dx_cat, levels = c(0,1)), 
    sex = factor(sex)
)
predictor <- "Cortex_GMWM_fd"
outcomes <- c("nih_patterncomp_raw", "tbx_ls", "nih_tlbx_theta_psm", "nih_tlbx_theta_orrt")
data <- data %>% mutate(across(c(predictor, outcomes), as.numeric))

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

regression_results <- data.frame()

for (outcome in outcomes) {
  for (group in unique(data$dx_cat)) {
    df_group <- data %>% filter(dx_cat == group) %>%
      mutate(across(all_of(c(outcome, predictor)), standardize))
    model <- lm(as.formula(paste(outcome, "~", predictor, "+ interview_age + sex")), data = df_group)
    summary_model <- summary(model)
    result <- data.frame(
      outcome = outcome,
      group = group,
      predictor = predictor,
      standardized_B = coef(model)[2],
      T_value = summary_model$coefficients[2, "t value"],
      p_value = summary_model$coefficients[2, "Pr(>|t|)"],
      Adjusted_R2 = summary_model$adj.r.squared
    )
    regression_results <- bind_rows(regression_results, result)
  }
}
write_csv(regression_results, "Cortex_regression_results.csv")

# region-wise analysis
regions <- c("Frontal_GMWM_fd", "Temporal_GMWM_fd", 
             "Parietal_GMWM_fd", "Occipital_GMWM_fd", 
             "Cingulate_GMWM_fd", "Insula_GMWM_fd")

all_regression_results <- data.frame()

for (region in regions) {
    for (outcome in outcomes) {
        for (group in unique(data$dx_cat)) {
            df_group <- data %>% filter(dx_cat == group) %>%
                mutate(across(all_of(c(outcome, region)), standardize))
            model <- lm(as.formula(paste(outcome, "~", region, "+ interview_age + sex")), data = df_group)
            summary_model <- summary(model)
            result <- data.frame(
                region = region,
                outcome = outcome,
                group = group,
                predictor = region,
                standardized_B = coef(model)[2],
                T_value = summary_model$coefficients[2, "t value"],
                p_value = summary_model$coefficients[2, "Pr(>|t|)"],
                Adjusted_R2 = summary_model$adj.r.squared
            )
            all_regression_results <- bind_rows(all_regression_results, result)
        }
    }
}

output_path <- "GMWM_regression_results.csv"
write_csv(all_regression_results, output_path)

GM_data <- read.csv("GM_data_with_demographics.csv")
data <- data %>%
    left_join(weighted_data[, c(1, 21:ncol(weighted_data))], by = "Subject")
data <- data %>% filter(dx_cat %in% c(0,1))
data <- data %>% mutate(
    dx_cat = factor(dx_cat, levels = c(0,1)),
    sex = factor(sex)
)
outcomes <- c("nih_patterncomp_raw", "tbx_ls")
GMvariables <- c("weighted_thickness_", "total_area_")
regions <- c("Cortex", "Frontal", "Cingulate", "Temporal", "Parietal", "Occipital", "Insula")

mediation_results <- data.frame()

for (outcome in outcomes) {
    for (GMvariable in GMvariables) {
        for (region in regions) {
            for (group in c(0, 1)) {
                data_filtered <- data %>%
                    filter(
                        !is.na(.data[[outcome]]),
                        !is.na(dx_cat),
                        !is.na(interview_age),
                        !is.na(sex),
                        dx_cat == group
                    )
                formula_str1 <- sprintf("%s_GMWM_fd ~ %s%s + interview_age + sex", region, GMvariable, region)
                my_formula1 <- as.formula(formula_str1)
                model_M <- lm(my_formula1, data = data_filtered)
                formula_str2 <- sprintf("%s ~ %s_GMWM_fd + %s%s + interview_age + sex", outcome, region, GMvariable, region)
                my_formula2 <- as.formula(formula_str2)
                model_Y <- lm(my_formula2, data = data_filtered)
                mediator <- sprintf("%s_GMWM_fd", region)

                set.seed(1234)
                med_analysis <- mediate(
                    model_M, model_Y,
                    treat = "interview_age",
                    mediator = mediator,
                    boot = TRUE, sims = 5000
                )

                model_M_summary <- summary(model_M)
                model_Y_summary <- summary(model_Y)

                a_coef <- coef(model_M)["interview_age"]
                a_p <- model_M_summary$coefficients["interview_age", "Pr(>|t|)"]
                a_se <- model_M_summary$coefficients["interview_age", "Std. Error"]
                a_ci <- confint(model_M)["interview_age", ]
                a_r2 <- model_M_summary$r.squared
                a_f <- model_M_summary$fstatistic[1]
                a_f_p <- pf(model_M_summary$fstatistic[1], model_M_summary$fstatistic[2], model_M_summary$fstatistic[3], lower.tail = FALSE)
                b_coef <- coef(model_Y)[mediator]
                b_p <- model_Y_summary$coefficients[mediator, "Pr(>|t|)"]
                b_se <- model_Y_summary$coefficients[mediator, "Std. Error"]
                b_ci <- confint(model_Y)[mediator, ]
                c_prime_coef <- coef(model_Y)["interview_age"]
                c_prime_p <- model_Y_summary$coefficients["interview_age", "Pr(>|t|)"]
                c_prime_se <- model_Y_summary$coefficients["interview_age", "Std. Error"]
                c_prime_ci <- confint(model_Y)["interview_age", ]
                c_r2 <- model_Y_summary$r.squared
                c_f <- model_Y_summary$fstatistic[1]
                c_f_p <- pf(model_Y_summary$fstatistic[1], model_Y_summary$fstatistic[2], model_Y_summary$fstatistic[3], lower.tail = FALSE)

                result <- data.frame(
                    group = group,
                    outcome = outcome,
                    mediator = mediator,
                    GMvariable = GMvariable,
                    a_coef = a_coef, a_p = round(a_p, 10), a_se = round(a_se, 10),
                    a_CI = paste0("[", round(a_ci[1], 10), ", ", round(a_ci[2], 10), "]"),
                    b_coef = b_coef, b_p = round(b_p, 10), b_se = round(b_se, 10),
                    b_CI = paste0("[", round(b_ci[1], 10), ", ", round(b_ci[2], 10), "]"),
                    c_prime_coef = c_prime_coef, c_prime_p = round(c_prime_p, 10), c_prime_se = round(c_prime_se, 10),
                    c_prime_CI = paste0("[", round(c_prime_ci[1], 10), ", ", round(c_prime_ci[2], 10), "]"),
                    a_r2 = round(a_r2, 10), a_f = round(a_f, 10), a_f_p = round(a_f_p, 10),
                    c_r2 = round(c_r2, 10), c_f = round(c_f, 10), c_f_p = round(c_f_p, 10),
                    ACME = med_analysis$d0,
                    ACME_se = round(sd(med_analysis$d0.sims), 10),
                    ACME_CIi = round(med_analysis$d0.ci[1], 10),
                    ACME_CIu = round(med_analysis$d0.ci[2], 10),
                    ADE = med_analysis$z0,
                    ADE_se = round(sd(med_analysis$z0.sims), 10),
                    ADE_CIi = round(med_analysis$z0.ci[1], 10),
                    ADE_CIu = round(med_analysis$z0.ci[2], 10),
                    Total_Effect = med_analysis$tau.coef,
                    Total_se = round(sd(med_analysis$tau.sims), 10),
                    Total_CIi = round(med_analysis$tau.ci[1], 10),
                    Total_CIu = round(med_analysis$tau.ci[2], 10)
                )
                mediation_results <- bind_rows(mediation_results, result)
            }
        }
    }
}

create_empty_rows <- function(n, cols) {
  empty_rows <- as.data.frame(matrix(NA, nrow = n, ncol = length(cols)))
  colnames(empty_rows) <- cols
  return(empty_rows)
}

mediation_results <- do.call(rbind, 
                                             lapply(1:nrow(mediation_results), function(i) {
                                               rbind(mediation_results[i, ], 
                                                     create_empty_rows(2, colnames(mediation_results)))
                                             }))

write.csv(mediation_results, "mediation_analysis_results.csv", row.names = FALSE)
