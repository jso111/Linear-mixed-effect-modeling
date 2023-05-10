# Linear-mixed-effect-modeling
Raw data and R code for performing linear mixed-effect modeling of vibratory tuning curve data

The file "linear mixed-effects example.R" can be run in RStudio. It loads in "testData.csv" and analyses it using the technique in the paper. You can use it in the same way, just format your data in a similar fashion and save as a CSV file. This is what most people would probably want to do. All the coefficients and p-values are printed out in RStudio, and then at the end, just the p-values for the Cohort comparisons for the magnitude and the phase data are given.

Column variables in the CSV file:
id - indentification code of the animal (text); 
cohort - the cohort of the animal (text); 
freq - the stimulus frequency (kHz); 
level - the stimulus level (dB SPL); 
mag - the vibratory magnitude (nm); 
phase - the vibratory phase (radians); 


If you want to recreate all figures in the paper, run the file "tuning curve statistics11.Rmd". It analyzes the data in "all_data.mat". Variables d1 and d2 are the two cohorts presented as dataset 2 in the paper. Variables d3 and d4 are the two cohorts presented as dataset 1 in the paper. 



