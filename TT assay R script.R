#By default:
# -Files with names containing "823" will be considered WT
# -Files with names containing "824" will be considered KO
# -If a file name contains "root" the value in "organ" column will be set to Root, otherwise to Shoot
# -Files with names containing "mock" will be considered 0h time point
# -Time point values are taken from the file names as a substring found between space and the letter "h"

#Clear the work space
rm(list=ls())

# Install and load the packages

install_and_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, dependencies = TRUE)
    library(package_name, character.only = TRUE)
  }
}

# List of required packages
required_packages <- c("dplyr", "ggplot2", "stringr", "readr", "data.table", "plyr", "purrr", "emmeans", "multcomp", "multcompView", "palmerpenguins")

# Install and load all required packages
lapply(required_packages, install_and_load)


# Ask user to select the experiment directory
Experiment_dir <- choose.dir(default = "", caption = "Select folder")

# Create a time-tamped subfolder
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(Experiment_dir, paste0("R_analysis_", timestamp))
dir.create(output_dir)

# Set the experiment folder as the working directory
setwd(Experiment_dir)

# List all .csv files in the experiment directory, excluding subfolders with results
csv_file_list <- list.files(path = Experiment_dir,
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = TRUE) %>%
  
# Filter out files that are in directories containing "R_analysis_" in the name
  .[!grepl("R_analysis_", dirname(.))]

# Merge the .csv files and save the merged results
merged_results <- csv_file_list %>%
  lapply(read_csv) %>%
  bind_rows

# Replace spaces in the column names with underscores, add columns for the time point, organ, and genotype
colnames(merged_results) <- gsub(" ","_", colnames(merged_results))
merged_results <- merged_results %>%
  mutate(Genotype = if_else(str_detect(File_name, '823'), 'WT', 'KO'),
         Organ = if_else(str_detect(File_name, 'root'), 'Root', 'Shoot'),
         Time_point_h = str_extract(File_name, "\\d+(?=h)"),
         Treatment = case_when(
           str_detect(File_name, 'AZD') ~ 'AZD',
           str_detect(File_name, '-N') ~ 'no_N',
           str_detect(File_name, '-C') ~ 'no_C',
           str_detect(File_name, 'mock') ~ 'mock',
           TRUE ~ NA_character_))
merged_results <- merged_results %>%
  mutate(Time_point_h = if_else(str_detect(File_name, 'mock'), '0', Time_point_h))

# Calculate the mean value for KO for each time point, organ, and treatment
KO_mean <- merged_results %>%
  filter(Genotype == "KO") 

# Calculate mean for each treatment/time point
KO_mean <- aggregate(Fluorescence_ratio_of_RFP_to_GFP ~ Organ + Time_point_h + Treatment,
                     data = KO_mean,
                     FUN = mean, na.rm = TRUE)
KO_mean <- KO_mean %>%
  dplyr::rename(Mean_ratio_for_KO = Fluorescence_ratio_of_RFP_to_GFP)

# Merge mean_KO data frame with the original data frame, organized based on "Organ", "Treatment" and "Time_point_h"
KO_mean_data <- merged_results %>%
  inner_join(KO_mean, by = c("Organ", "Time_point_h", "Treatment"))

# Normalize the RFP/GFP ratio to the mean value obtained for the KO in the same organ, at the same time point
normalized_data <- KO_mean_data %>%
  mutate(Normalized_ratio = Fluorescence_ratio_of_RFP_to_GFP / Mean_ratio_for_KO)

write.csv(normalized_data,file.path(output_dir, "Merged_results.csv"))

# Run Two-way-Anova and Tukey's HSD test on the normalized ratios for WT
WT_only <- normalized_data %>%
  filter(Genotype == "WT")

# Log-transform the Normalized_ratio column for Anova test
WT_only <- WT_only %>%
  mutate(Log_Normalized_ratio = log(Normalized_ratio))

# Save the modified dataframe as a CSV file
write.csv(WT_only, file.path(output_dir,"Merged_results_for_WT_only.csv"))

# Dynamically determine the levels for Time_point_h
time_points <- sort(unique(as.numeric(normalized_data$Time_point_h)))
normalized_data$Time_point_h <- factor(normalized_data$Time_point_h, levels = as.character(time_points))

# Run comparisons for each treatment (mock is included into each run)
perform_analysis <- function(data_for_analysis, trtm, output_dir) {
 
  # Run analysis of variance by Two-way Anova
  Two_way_Anova <- aov(Log_Normalized_ratio ~ Time_point_h * Organ, data = data_for_analysis)
  
  # Run means comparison by Tukey's HSD test
  Tukey <- TukeyHSD(Two_way_Anova)
  
  # Add CLD column to the summary stats
  cld <- multcompLetters4(Two_way_Anova, Tukey)  
  print(cld)
  
  # Create a table with factors and 3rd quantile
  Tukey_with_CLD <- data_for_analysis %>%
    group_by(Time_point_h, Organ) %>%
    dplyr::summarise(mean = mean(Log_Normalized_ratio), quant = quantile(Log_Normalized_ratio, probs = 0.75), .groups = 'drop') %>%
    arrange(desc(mean))
  
  # Add CLD to the Tk table
  cld_df <- as.data.frame.list(cld$`Time_point_h:Organ`)
  Tukey_with_CLD$cld <- cld_df$Letters
  write.csv(Tukey_with_CLD, file.path(output_dir,  paste0("Tukey_with_CLD_", trtm, ".csv")))
  
  # Build box plots for WT ratios
  
  # Sort organ categories
  data_for_analysis$Organ <- factor(data_for_analysis$Organ, c("Shoot", "Root"))
  
  # Ensure CLDs are sorted the same way as the organ categories
  Tukey_with_CLD$Organ <- factor(Tukey_with_CLD$Organ, levels = c("Shoot", "Root"))
  
  # Build a box plot on WT ratios
  
  #NB! CLDs were not calculated for this plot
  p1<-ggplot(data_for_analysis, aes(Time_point_h, Fluorescence_ratio_of_RFP_to_GFP, fill = Organ)) +
    geom_boxplot(width = 0.5, size = 0.15, outlier.shape = 1, outlier.size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = c("#83D40F", "#EA9A16")) +
    scale_y_continuous(limits = c(0, NA))+
    labs(title = paste0("WT, RFP/GFP ratio, ", trtm), x = "Elapsed time (h)", y = "RFP/GFP ratio")
  ggsave(file.path(output_dir,paste0("Box_plot_ratio_", trtm, ".pdf")), width = 5, height = 4)
  # Print the plot to RStudio’s Plots panel
  print(p1)
  
  # Build a box plot normalized WT ratios
  
  # Calculate the y-axis maximum for placing CLDs on the plot
  y_max <- max(data_for_analysis$Normalized_ratio, na.rm = TRUE)

  # Set a fixed y position for CLDs to be 5% above the y_max value
  Tukey_with_CLD <- Tukey_with_CLD %>%
    mutate(label_y_position = y_max*1.05)

    p2<-ggplot(data_for_analysis, aes(Time_point_h, Normalized_ratio, fill = Organ)) +
    geom_boxplot(width = 0.5, size = 0.15, outlier.shape = 1, outlier.size = 1) +
    theme_bw() + 
    geom_text(data = Tukey_with_CLD, aes(x = Time_point_h, y = label_y_position, label = cld), size = 3, vjust = 0, hjust = 0.5, position = position_dodge(width = 0.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = c("#83D40F", "#EA9A16")) +
    scale_y_continuous(limits = c(0, NA))+
    labs(title = paste0("WT, RFP/GFP ratio normalized to KO, ", trtm), x = "Elapsed time (h)", y = "Normalized RFP/GFP ratio")
  ggsave(file.path(output_dir,paste0("Box_plot_normalized_ratio_", trtm, ".pdf")), width = 5, height = 4)
  # Print the plot to RStudio’s Plots panel
  print(p2)
}

# Get the unique treatments, excluding 'mock'
unique_treatments <- unique(WT_only$Treatment[WT_only$Treatment != "mock"])

# Perform analysis for each unique treatment, including 'mock' in the analysis
for (trtm in unique_treatments) {
  if (!is.na(trtm)) {
    
    # Subset the data for the current treatment and 'mock'
    data_for_analysis <- WT_only %>%
      filter(Treatment == trtm | Treatment == "mock")
    
    # Perform analysis on the combined data (current treatment + 'mock')
    perform_analysis(data_for_analysis, trtm, output_dir)
  }
}