# Breast-Cancer-Prediction-Project

This project explores the BRCA dataset from the dslabs package, focusing on breast cancer diagnosis through biopsy samples. By analyzing the characteristics of cell nuclei, we aim to classify tumors as benign or malignant, utilizing advanced statistical techniques and machine learning methods. The dataset consists of a range of numeric features that provide insights into the shape and size of the nuclei, offering a rich ground for predictive modeling. Through this project, I will demonstrate the application of data manipulation, visualization, and model evaluation techniques in R. This work not only contributes to understanding breast cancer diagnosis but also serves as a practical example of data science in healthcare. Join me on this journey to uncover patterns and improve diagnostic accuracy!

The brca dataset from the dslabs package contains information about breast cancer diagnosis biopsy samples for tumors that were determined to be either benign (not cancer) and malignant (cancer). The brca object is a list consisting of:
  brca$y: a vector of sample classifications ("B" = benign or "M" = malignant)
  brca$x: a matrix of numeric features describing properties of the shape and size of cell nuclei extracted from biopsy microscope images
For these exercises, load the data by setting your options and loading the libraries and data as shown in the code here:
  options(digits = 3)
  library(matrixStats)
  library(tidyverse)
  library(caret)
  library(dslabs)
  data(brca)
  
Remember that you can use the R help files to learn more about functions you are less familiar with. The help files can be called using the ? symbol. For example, ?mean will load the help file for the mean() function.

IMPORTANT: Make sure your package is up to date with the command install.packages("dslabs").
