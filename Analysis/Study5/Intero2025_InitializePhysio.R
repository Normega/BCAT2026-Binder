rm(list = ls())

# -- PATH CONFIGURATION -------------------------------------------------
# Set ROOT_DIR to the repository root before running.
# Default (".") assumes the script is run from the repo root directory.
if (!exists("ROOT_DIR")) ROOT_DIR <- "."
mainPath     <- ROOT_DIR
dataPath     <- file.path(ROOT_DIR, "data")
taskPath     <- file.path(dataPath, "Behaviour")
physioPath   <- file.path(dataPath, "Physio")
analysisPath <- file.path(ROOT_DIR, "study5_processing")
resultsPath  <- file.path(ROOT_DIR, "results")
# -----------------------------------------------------------------------

qualtricsFile     <- paste0(dataPath,    'InteroceptionSummer2025_Qualtrics.xlsx')
conditionFile     <- paste0(dataPath,    'ConditionLookup.xlsx')
questionnaireFile <- paste0(resultsPath, 'questionnaireFile.csv')
taskDataFile      <- paste0(resultsPath, 'dataFile.csv')
taskTestFile      <- paste0(resultsPath, 'testFile.csv')
qcFile            <- paste0(resultsPath, 'qcFile.xlsx')
hrFile            <- paste0(resultsPath, 'hrFile.csv')
hrReportFile      <- paste0(resultsPath, 'hrReport.csv')

pipelineFile  <- file.path(analysisPath, "breath_pipeline.R")
functionsFile <- file.path(analysisPath, "Intero2025_RespirationFunctions.R")
source(pipelineFile)
source(functionsFile)

# CONSTANTS -------------------------------------------------------
STARTDUR   <- 4     # Starting duration of each breath trial (s)
LAG        <- 0.5   # Signal transduction correction (s)
NUMBREATHS <- 4     # Number of breaths per trial


# Set up ----------------------------------------------------------
packages <- c(
  "readxl", "writexl",
  "ggplot2", "ggeffects", "ggforce", "GGally", "ggpubr",
  "patchwork",
  "corrplot", "psych",
  "lme4", "lmerTest",
  "DataCombine", "sjPlot",
  "tidyverse", "tidyr", "stringr",
  "reticulate",
  "signal",
  "pracma", "tseries"
)

if (length(setdiff(packages, rownames(installed.packages()))) > 0)
  install.packages(setdiff(packages, rownames(installed.packages())))

options(readr.num_columns = 0)
for (thispack in packages)
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)

reticulate::py_install("bioread")
bioread <- import("bioread")


# Load files ------------------------------------------------------
condData <- read_excel(conditionFile) %>%
  dplyr::distinct() %>%
  dplyr::select(id, Computer, Condition) %>%
  dplyr::rename(firstCondition = Condition)
condData$id <- as.numeric(condData$id)
table(condData$firstCondition)

longData <- read.csv(taskDataFile)
longTest <- read.csv(taskTestFile)

longData <- left_join(longData, condData, by = "id")

longData$currentCondition <-
  factor(longData$firstCondition == "Visual" & longData$ses == 1,
         labels = c("Breath", "Visual"))

condLookup <- longData %>%
  group_by(id, ses) %>%
  dplyr::select(currentCondition) %>%
  distinct()

longTest <- left_join(longTest, condLookup, by = c("id", "ses"))

testingIds <- c("1234", "5678")
longData   <- longData[!longData$id %in% testingIds, ]
longTest   <- longTest[!longTest$id %in% testingIds, ]

idList   <- unique(longData$id)
fileList <- list.files(path = physioPath)
