#By default:
# -files with names containing "823" will be considered WT
# -files with names containing "824" will be considered KO
# -If a file name contains "root" the organ will be set to Root, otherwise to Shoot
# -files with names containing "mock" will be considered 0h time point
# -time point values are taken from the file names as a value preceding letter "h"


#Install and load the packages
  #install.packages("dplyr")
  #install.packages("plyr")
  #install.packages("readr")
  #install.packages(data.table)
  #install.packages("stringr")
  #install.packages("purrr")
  #install.packages("ggplot2")
  #install.packages("palmerpenguins")
  #install.packages("emmeans")
  #install.packages("multcomp")
  #install.packages("multcompView")
  library(plyr)
  library(dplyr)
  library(readr)
  library(data.table)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(stats)
  library(palmerpenguins)
  library(emmeans)
  library(multcomp)
  library(multcompView)

# Ask user to select the experiment directory
  Experiment_dir <- choose.dir(default = "", caption = "Select folder")

# Set the experiment folder as the working directory
  setwd(Experiment_dir)

# List all .csv files in the experiment directory
  csv_file_list <- list.files(path = Experiment_dir,
                            recursive = TRUE,
                            pattern = "\\.csv$")

# merge the .csv files and save the merged results
  merged_results <- csv_file_list %>%
  lapply(read_csv) %>%
  bind_rows
  
# Replace spaces in the column names with underscores, add columns for the time point, organ and genotype
  colnames(merged_results) <- gsub(" ","_", colnames(merged_results))
  merged_results <- merged_results %>%
  mutate(Genotype = if_else(str_detect(File_name, '823'), 'WT', 'KO'))
  merged_results <- merged_results %>%
  mutate(Organ = if_else(str_detect(File_name, 'root'), 'Root', 'Shoot'))
  merged_results <- merged_results %>%
  mutate(Time_point_h = if_else(str_detect(File_name, 'mock'), '0', NA_character_))
  merged_results <- merged_results %>%
  mutate(Time_point_h = ifelse(str_detect(File_name, 'h'), sub('h.*', '', File_name), Time_point_h))


# Calculate the mean value for KO for each time point and organ
  KO_mean <- merged_results %>%
  filter(Genotype == "KO") %>%
  group_by(Organ, Time_point_h) %>%
  summarise_at(vars(Fluorescence_ratio_of_RFP_to_GFP), mean)
# Rename the column of the KO mean
  colnames(KO_mean)[colnames(KO_mean) == "Fluorescence_ratio_of_RFP_to_GFP"] <- "Mean_ratio_for_KO"
 
# Merge mean_KO data frame wit the original data frame, organized based on "Organ" and "Time_point_h"
  KO_mean_data <- merged_results %>%
  inner_join(KO_mean, by = c("Organ", "Time_point_h"))
  
# Normalize the RFP/GFP ratio to the mean value obtained for the KO in the same organ, at the same time point
  normalized_data <- KO_mean_data %>%
  mutate(Normalized_ratio = Fluorescence_ratio_of_RFP_to_GFP / Mean_ratio_for_KO)
  write.csv(normalized_data, "Merged results.csv")
  
#Run statistics
  # Run Two-way-Anova and Tukey's HSD test on the normalized ratios for WT
    WT_only <- normalized_data %>%
    filter(Genotype == "WT")
    
    # Log-transform the Normalized_ratio column for Anova test
    WT_only <- WT_only %>%
    mutate(Log_Normalized_ratio = log(Normalized_ratio))
    
    # Save the modified dataframe as a CSV file
    write.csv(WT_only, "Merged results for WT only.csv")
    
    # Run anlaysis of variance by Two-way Anova
    Two_way_Anova <- aov(Log_Normalized_ratio ~ Time_point_h * Organ, data = WT_only)
   
    # Run means comparison by Tukey's HSD test
    Tukey <- TukeyHSD(Two_way_Anova)
    
    
  # Add CLD column to the summary stats
    cld <-multcompLetters4(Two_way_Anova,Tukey)  
    print(cld)
    
  # Create a table with factors and 3rd quantile
    Tukey_with_CLD <- WT_only %>%
      group_by(Time_point_h, Organ) %>%
      summarise(mean = mean(Log_Normalized_ratio), quant = quantile(Log_Normalized_ratio, probs = 0.75)) %>%
      arrange(desc(mean))
    
    # Add CLD to the Tk table
    cld_df <- as.data.frame.list(cld$`Time_point_h:Organ`)
    Tukey_with_CLD$cld <- cld_df$Letters
    write.csv(Tukey_with_CLD, "Tukey_with_CLD.csv")
    


# Build box plots for WT ratios
  
  # Sort organ categories
    normalized_data$Organ <- factor(normalized_data$Organ,c("Shoot", "Root"))
  
  # Ensure CLDs are sorted the same way as the organ categories
    Tukey_with_CLD$Organ <- factor(Tukey_with_CLD$Organ, levels = c("Shoot", "Root"))
  
  # Ensure that time points are taken in numerical order
    normalized_data$Time_point_h <- factor(normalized_data$Time_point_h, c("0", "2", "16", "24", "36", "48"))
 
  # Ensure only WT ratios are plotted
    WT_only <- normalized_data %>%
    filter(Genotype == "WT")
    
  #build a box plot on WT ratios
      ggplot(WT_only , aes(Time_point_h, Fluorescence_ratio_of_RFP_to_GFP, fill = Organ)) +
      geom_boxplot(width = 0.5, size = 0.15,outlier.shape = 1, outlier.size = 1) +
      geom_text(data = Tukey_with_CLD, aes(x = Time_point_h, y = 3, label = cld), size = 3, vjust=-10, hjust =0.5, position = position_dodge(width = 0.5)) +
      theme_bw() +
      theme(panel.grid.major = element_line(size = 0.1),panel.grid.minor = element_line(size = 0.1)) +
      scale_fill_manual(values = c("#83D40F", "#EA9A16")) +
      labs(title = "WT, RFP/GFP ratio", x = "Elapsed time (h)", y = "RFP/GFP ratio")
      ggsave("Box_plot_ratio.pdf", width = 5, height = 4)
    
    
    
    #build a box plot normalized WT ratios
      ggplot(WT_only , aes(Time_point_h, Normalized_ratio, fill = Organ)) +
      geom_boxplot(width = 0.5, size = 0.15,outlier.shape = 1, outlier.size = 1) +
      theme_bw() + 
      geom_text(data = Tukey_with_CLD, aes(x = Time_point_h, y = 2, label = cld), size = 3, vjust=-20, hjust =0.5, position = position_dodge(width = 0.5)) +
      theme(panel.grid.major = element_line(size = 0.1),panel.grid.minor = element_line(size = 0.1))+
      scale_fill_manual(values = c("#83D40F", "#EA9A16"))+
      labs(title = "WT, RFP/GFP ratio normalized to KO", x = "Elapsed time (h)", y = "Normalized RFP/GFP ratio")
      ggsave("Box_plot_normalized_ratio.pdf", width = 5, height = 4)