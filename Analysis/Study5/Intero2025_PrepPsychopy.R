rm(list=ls())

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

questionnaireFile = paste0(resultsPath, "questionnaireFile.csv")
taskDataFile = paste0(resultsPath, "dataFile.csv")
taskTestFile = paste0(resultsPath, "testFile.csv")
hrReportFile = paste0(resultsPath, "hrReport.csv")



# Set Up---------
## Load libraries ---------
packages <- c("readxl",
              "ggplot2","ggeffects","ggforce","GGally","ggpubr",
              "corrplot",
              "psych",
              "lme4", "lmerTest", 
              "DataCombine",
              "sjPlot",
              "tidyverse",
              "tidyr",
              "stringr",
              "reticulate", #for using python to read in .acq files
              "signal", #for working with timeseries data like breath
              "pracma" #for calculating timeseries rates
)

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}
options(readr.num_columns = 0)
for (thispack in packages) {
  library(thispack,character.only=TRUE,quietly=TRUE,verbose=FALSE)
}



#TASK DATA --------------------

fileList = list.files(taskPath, pattern=glob2rx('*Intero2025*.csv'))
#fileList = fileList[file.info(paste0(mypath,fileList))$size > 0]
fileList = fileList[file.info(paste0(taskPath,fileList))$size > 40000]

ids <- unname(sapply(fileList, function(thisFile) {
  substr(thisFile, 1, regexpr("_", thisFile)[1] - 1)
}))

table(ids)[table(ids) > 2] #check if there are more than 2 entries
sesnum = numeric()
sesnum[1] = 1
for (thisId in 2:length(ids)){
  lastId = thisId - 1
  if ( ids[thisId] == ids[lastId] ){
    sesnum[thisId] = 2
  } else {
    sesnum[thisId] = 1
  }
}

bigData = data.frame(id = ids,
                     ses = sesnum)
t = table(bigData)
sum(t[,1]) #222
sum(t[,2]) #160

## IF WE WANT TO CUT DO IT HERE
#cutlist = rep(0,length(ids))
#for (thisId in 2:length(ids)){
#  lastId = thisId - 1
#  if ( ids[thisId] == ids[lastId] ){
#    cutlist[lastId] = 1  
#  } 
#}

#fileList = fileList[!cutlist]
#ids = ids[!cutlist]


longData = data.frame()
longTest = data.frame()
hrRating = data.frame()

#idNum = 4
for (idNum in 1:length(ids)){
  thisId = ids[idNum]
  thisFile = fileList[idNum]
  thisSes = sesnum[idNum]
  print(paste(thisId, thisSes))  
  #Load data --------------
  mydata <- read.csv(paste0(taskPath,thisFile))
  #head(mydata)
  mydata = mydata[-1,]
  #View(mydata)
  
  varList = c("trials.label",  
              "trials.thisRepN", "level", "Condition",
              "Salience", "Direction", "DirectionLabel", "Correct", "Accuracy",
              "Response",
              "ArousalRating", "confidenceSlider.response",
              "trial.started","trial.stopped")
  trialData = mydata[!is.na(mydata$trials.thisN),varList]
  trialData$id = thisId
  trialData$ses = thisSes
  
  testVarList = c("testLabel",  
                  "testTrials.thisN", "level", "Condition",
                  "Salience", "Direction", "DirectionLabel", "Correct", "Accuracy",
                  "Response",
                  "ArousalRating", "confidenceSlider.response",
                  "TestTrial.started", "TestTrial.stopped")
  testData = mydata[!is.na(mydata$testTrials.thisN),testVarList]
  testData$id = thisId
  testData$ses = thisSes
  
  hrVarList = c("thisCondition", "countDuration", "HeartbeatLoop.thisN", 
                "countScreen.started", "countScreen.stopped",
                "HRreport")
  if (thisSes == 1 & ("HRreport" %in% names(mydata))){
    hrData = mydata[!is.na(mydata$HeartbeatLoop.thisN),hrVarList]
    hrReport = data.frame(id = thisId,
                          Onset25 = hrData$countScreen.started[hrData$thisCondition=="short"],
                          Offset25 = hrData$countScreen.stopped[hrData$thisCondition=="short"],
                          Onset35 = hrData$countScreen.started[hrData$thisCondition=="medium"],
                          Offset35 = hrData$countScreen.stopped[hrData$thisCondition=="medium"],
                          Onset55 = hrData$countScreen.started[hrData$thisCondition=="long"],
                          Offset55 = hrData$countScreen.stopped[hrData$thisCondition=="long"],
                          short25Count = hrData$HRreport[hrData$thisCondition=="short"],
                          medium35Count = hrData$HRreport[hrData$thisCondition=="medium"],
                          long55Count = hrData$HRreport[hrData$thisCondition=="long"]
    )
    hrRating = rbind(hrRating, hrReport)
  } 
  
  
  longData = rbind(longData, trialData)
  longTest = rbind(longTest, testData)
  
}

#names(longData)
#names(longTest)

names(longData)[which( names(longData) %in% "trials.label")] = "taskCondition"
names(longTest)[which( names(longTest) %in% "testLabel")] = "taskCondition"
names(longData)[which( names(longData) %in% "JudgeRating")] = "Confidence"
names(longTest)[which( names(longTest) %in% "JudgeRating")] = "Confidence"
names(longData)[which( names(longData) %in% "ArousalRating")] = "Arousal"
names(longTest)[which( names(longTest) %in% "ArousalRating")] = "Arousal"
names(longData)[which( names(longData) %in% "condition")] = "Group"
names(longTest)[which( names(longTest) %in% "condition")] = "Group"
names(longData)[which( names(longData) %in% "trials.thisRepN")] = "Trial"
names(longTest)[which( names(longTest) %in% "trials.thisN")] = "Trial"
names(longData)[which( names(longData) %in% "confidenceSlider.response")] = "Confidence"
names(longTest)[which( names(longTest) %in% "confidenceSlider.response")] = "Confidence"
#names(longData)[which( names(longData) %in% "trials.response")] = "Response"
#names(longTest)[which( names(longTest) %in% "trials.response")] = "Response"

longData$Change = longData$level * longData$Direction
longData$Direction = factor(longData$Direction, labels = c("Faster", "NoChange", "Slower"))
longData$Salience = factor(longData$Salience, labels = c("Low", "High"))


View(longData)

write.csv(longData,file = taskDataFile, row.names = FALSE)
write.csv(longTest,file = taskTestFile, row.names = FALSE)
write.csv(hrRating,file = hrReportFile, row.names = FALSE)

