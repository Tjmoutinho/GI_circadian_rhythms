# Analyze kronos circadian rhythms 
# Written by Deborah H. Luzader (dh3bj@virginia.edu) and Thomas J. Moutinho Jr. (tjm4k@virginia.edu)

#### Packages ####
# install.packages("readxl")
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("DBI")
# install.packages("lazyeval")
# install.packages("quantmod")
# install.packages("ggplot2")
# install.packages("ggthemes")
# install.packages("xlsx")
# install.packages("rJava")

#### Libraries ####
library(rJava)
library('quantmod')
library(dplyr)
library('readxl')
library(tidyr)
# library(lazyeval)
library(ggplot2)
# library(ggthemes)
library("xlsx")

source("Functions_Circadian_rhythms.R")

#### INSTRUCTIONS ####

# Enter information as follows in the paratheses after WRITE_DATA_FILES
# File name you created for your excel data in quotations. EX: "20171003per2 metabolites.xlsx"   ***Make sure that the excel file you are trying to analyze is in the same folder as the code files ("Final_Circ_rhyth_analysis.R" and "Functions_Circadian_rhythms.R")
# Make sure that if you are taking the data from the Kronos machine, that it is detrended, noise filtered, and background is removed
# Enter the number corresponding to each sheet that you want analyzed. Make sure that the first sheet is the control in your experiment ex: c(1,2,3,4,5)
# Enter the first time point in hours that you want analyzed. EX: 13 
# Enter the last time point in hours that you want analyzed. EX: 72 
# The file name "ALL_DATA.xlsx" will output your normalized luminesence data, raw luminesence data, average period, amplitude and phase
# The file name "ALL_PERIOD_DATA.xlsx" will output all your period data, which is calculated as 1/4 of each period (so peak to zero, zero to trough, trought to zero, zero to peak, and so on)
# The file name "ALL_PHASE_DATA.xlsx" will output all phase data, which is calculated by comparing each peak, zero, and trough between the control (which is sheet 1) and your sheet
# Note that the output data will be labeled as the sheet number.
# Hit the button Run in the top right corner

WRITE_DATA_FILES("Data/20180820 ent per2 asf 492 500 502 original KRONOS data.xlsx", c(1,2,3,4,5,6,7), 13, 72,
                 "Data/20180820 ent per2 asf 492 500 502 original KRONOS data_ALL_DATA.xlsx", 
                 "Data/20180820 ent per2 asf 492 500 502 original KRONOS data_ALL_PRT_DATA.xlsx.xlsx", 
                 "Data/20180820 ent per2 asf 492 500 502 original KRONOS data_ALL_PHASE_DATA.xlsx",
                 "Data/20180820 ent per2 asf 492 500 502 original KRONOS data_ALL_PRT_DATA_CORRECTED.xlsx")

