# TT (Tandem Tag) Assay for Autophagic Activity

This TT assay was developed to measure autophagic activity in the roots and shoots of plants expressing GFP–RFP–ATG8 (TT–ATG8).

## TT ImageJ Macro

The ImageJ macro is designed for semi-automated analysis of confocal images in Carl Zeiss format (.czi). The macro will:
- Analyze all .czi files within a user-selected folder.
- Allow the user to select any number of ROIs (Regions of Interest) per image and adjust their size/placement.
- Measure the intensity of GFP and RFP within each ROI for each image.
- Calculate intensity/area and the RFP/GFP intensity ratio (an increased ratio indicates increased autophagic activity).
- Save measurement results in a .csv file.
- Save ROI sets for each image.
- Automatically load saved ROI sets and allow repositioning/adjustment if the macro is re-run on the same folder.

<p align="center">
  <img src="https://github.com/AlyonaMinina/TT-assay-for-autophagic-activity/blob/4da980b293e892e50ea79d038d63f10ddb5336d2/Readme%20files/TT%20example.png" height="300" title="TT example">
</p>

## TT R Script

The R script is designed to process the quantitative output from the TT macro and compare mock vs. treatment for each genotype, as well as shoots vs. roots. The R script will:

- Ask the user to choose the folder with TT macro results.
- Use information from file names to determine treatment duration, type, genotype, and imaged plant organ. (NB: Edit the script or ensure correct file naming before analysis.)
 -Obtain mean values for the ATG knockout and use them to normalize corresponding wild-type plant values.
- Run ANOVA, treating time points as factors.
- Save Tukey HSD test results and Compact Letter Displays (CLDs) in a .csv file.
- Generate box plots for RFP/GFP ratios for each treatment and save them as PDFs.
- Save results into a time-stamped subfolder, which will be excluded from analysis if the script is re-run on the same experiment folder.

<p align="center">
  <img src="https://github.com/AlyonaMinina/TT-assay-for-autophagic-activity/blob/4da980b293e892e50ea79d038d63f10ddb5336d2/Readme%20files/Box_plot_normalized_ratio_no_N.png" height="300" title="TT plot example">
</p>
