rm(list=ls())

#Study folders
# winpath = 'I:/Shared drives/Interoception 2025/Data/'
# taskPath = paste0(winpath, "Behaviour/") 
# qualtricsPath = paste0(winpath, "InteroceptionSummer2025_Qualtrics.xlsx") 
# conditionPath = paste0(winpath, "ConditionLookup.xlsx")
# resultsPath = "I:/Shared drives/Interoception 2025/Results/"

winpath      <- "i:/Shared drives/Interoception 2025/"
dataPath     <- file.path(winpath, "Data")
resultsPath  <- file.path(winpath, "Results")
taskPath      <- file.path(dataPath, "Behaviour")
qualtricsPath <- file.path(dataPath, "InteroceptionSummer2025_Qualtrics.xlsx")
conditionPath <- file.path(dataPath, "ConditionLookup.xlsx")

#specific data files 
questionnaireFile = file.path(resultsPath, "questionnaireFile.csv")
taskDataFile = file.path(resultsPath, "dataFile.csv")
taskTestFile = file.path(resultsPath, "testFile.csv")

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
              "stringr"
)

#Checks if libraries are installed and install as needed
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}
options(readr.num_columns = 0)

#Loads the libraries 
for (thispack in packages) {
  library(thispack,character.only=TRUE,quietly=TRUE,verbose=FALSE)
}


#Load in self-report data -------
qData <- read_excel(qualtricsPath, sheet = "Data")
qData <- qData[-1, ]                          # drop Qualtrics question-text row
qData <- qData[qData$Finished == "True", ]   # completed responses only


#Replace column names 
varNames = c("StartDate", "EndDate", "Status", "IPAddress", "Progress", "Duration", 
             "Finished", "RecordedDate", "ResponseID", "RecipientLastName", "RecipientFirstName", 
             "RecipientEmail", "ExternalReference", "LocationLatitude", "LocationLongitude", 
             "DistributionChannel", "UserLanguage", "sonaID_1", "sonaID_2", "BFI_Reserved", 
             "BFI_Trusting", "BFI_Lazy", "BFI_Relaxed", "BFI_Artistic", "BFI_Outgoing", 
             "BFI_FaultWithOthers", "BFI_ThoroughJob", "BFI_NervousEasily", 
             "BFI_Activeimagination", "SWLS", "WBLadder", "SPANE_Positive", "SPANE_Negative", 
             "SPANE_Good", "SPANE_Bad", "SPANE_Pleasant", "SPANE_Unpleasant", "SPANE_Happy", 
             "SPANE_Sad", "SPANE_Afraid", "SPANE_Joyful", "SPANE_Angry", "SPANE_Contented", 
             "BIPS_Things", "BIPS_Hurry", "BIPS_Pressure", 
             "Catch_1", "BIPS_Conflict", "BIPS_Hadto", "BIPS_Criticized", "BIPS_Difficulties", 
             "BIPS_OnTop", "BIPS_Worries", "BriefMAIA_Uncomfortable", "BriefMAIA_Comfortable", 
             "BriefMAIA_Changes", "BriefMAIA_Distract", "BriefMAIA_Discomfort", 
             "BriefMAIA_Unpleasant", "BriefMAIA_PhysicalPain", "BriefMAIA_Worry", 
             "BriefMAIA_DiscomfortPain", "BriefMAIA_MaintainAwareness", 
             "BriefMAIA_ReturnAwareness", "BriefMAIA_RefocusAttention", 
             "BriefMAIA_NoticeBodyDifferent", "BriefMAIA_NoticeBreathingChanges", 
             "BriefMAIA_NoticeBodyChanges", "BriefMAIA_AwarenessBody", 
             "BriefMAIA_BreathTension", "BriefMAIA_Thoughts", "BriefMAIA_Information", 
             "BriefMAIA_Upset", "BriefMAIA_Listen", "BriefMAIA_HomeBody", "BriefMAIA_BodySafe", 
             "BriefMAIA_TrustBody", "BARQR_Tense", "BARQR_Affected", "BARQR_NotAwareBreathe", 
             "BARQR_AttentionMove", "BARQR_StruggleToRelax", "BARQR_BodyTense", 
             "BARQR_Feeling", "BARQR_Digestion", "BARQR_Comfortable", "BARQR_Unpredictable", 
             "BARQR_AttentionBody", "BARQR_Touched", "Catch_2", "PHQ4_Nervous", "PHQ4_Worry", 
             "PHQ4_Depressed", "PHQ4_LittleInterest", "MSES12_Doubt", "MSES12_Unhappy", 
             "MSES12_SelfConscious", "MSES12_Uncomfortable", "MSES12_Loser", "MSES12_Criticize", 
             "MSES12_GoodJob", "MSES12_Convinced", "MSES12_Ashamed", "MSES12_Attractive", 
             "MSES12_SportsPerformance", "MSES12_SportsNervous", "PAQS_BadRightWords", 
             "PAQS_BadSadAngryScared", "PAQS_Ignore", "PAQS_GoodRightWords", 
             "PAQS_GoodHappyExcitedAmused", "PAQS_AttentionToEmotions", "Age", 
             "YearofStudy", "Gender_Male", "Gender_Female", "Gender_NB", "Gender_Transgender", 
             "Gender_Other", "Gender_PreferNotToAnswer", "Gender_Other_Text", 
             "Identity_Indigenous", "Identity_POC", "Cultural_Identity_Arab", 
             "Cultural_Identity_Black", "Cultural_Identity_Chinese", "Cultural_Identity_Filipino", 
             "Cultural_Identity_Indigenous", "Cultural_Identity_Japanese", 
             "Cultural_Identity_Korean", "Cultural_Identity_LatinAmerican", 
             "Cultural_Identity_Mixed", "Cultural_Identity_SouthAsian", 
             "Cultural_Identity_SouthEastAsian", "Cultural_Identity_WestAsian", 
             "Cultural_Identity_White", "Cultural_Identity_ThirdCulturalIdentity", 
             "Cultural_Identity_SelfIdentify", "Cultural_Identity_PreferNotToAnswer", 
             "Cultural_Identity_SelfIdentifySpecify", "CountryofBirth", "CommunityLadder", 
             "English_Primary", "Languages_Other", "ParentalEducation_Select", 
             "ParentalEducation_Other", "FirstGen_Select", "FirstGen_Other", "Program_Select", 
             "Program_Other", "MaritalStatus", "EmploymentStatus", "JobStatus", 
             "HouseholdIncome", "Residence", "Off_campusHousing", "LivingwithFamily", 
             "CatchTotal")


names(qData) = varNames
#View(qData)

#List which columns go with which questionnaire
MAIAcols = c("BriefMAIA_Uncomfortable", "BriefMAIA_Comfortable", 
             "BriefMAIA_Changes", "BriefMAIA_Distract", "BriefMAIA_Discomfort", 
             "BriefMAIA_Unpleasant", "BriefMAIA_PhysicalPain", "BriefMAIA_Worry", 
             "BriefMAIA_DiscomfortPain", "BriefMAIA_MaintainAwareness", 
             "BriefMAIA_ReturnAwareness", "BriefMAIA_RefocusAttention", 
             "BriefMAIA_NoticeBodyDifferent", "BriefMAIA_NoticeBreathingChanges", 
             "BriefMAIA_NoticeBodyChanges", "BriefMAIA_AwarenessBody", 
             "BriefMAIA_BreathTension", "BriefMAIA_Thoughts", "BriefMAIA_Information", 
             "BriefMAIA_Upset", "BriefMAIA_Listen", "BriefMAIA_HomeBody", "BriefMAIA_BodySafe", 
             "BriefMAIA_TrustBody"
)

BIPScols = c("BIPS_Things", "BIPS_Hurry", "BIPS_Pressure", 
             "BIPS_Conflict", "BIPS_Hadto", "BIPS_Criticized", 
             "BIPS_Difficulties", "BIPS_OnTop", "BIPS_Worries"
)

MSEScols = c("MSES12_Doubt", "MSES12_Unhappy", 
             "MSES12_SelfConscious", "MSES12_Uncomfortable", "MSES12_Loser", "MSES12_Criticize", 
             "MSES12_GoodJob", "MSES12_Convinced", "MSES12_Ashamed", "MSES12_Attractive", 
             "MSES12_SportsPerformance", "MSES12_SportsNervous")

AlexithymiaCols = c( "PAQS_BadRightWords", 
                 "PAQS_BadSadAngryScared", "PAQS_Ignore", "PAQS_GoodRightWords", 
                 "PAQS_GoodHappyExcitedAmused", "PAQS_AttentionToEmotions")

PHQcols = c("PHQ4_Nervous", "PHQ4_Worry", 
            "PHQ4_Depressed", "PHQ4_LittleInterest")
BARQcols = c("BARQR_Tense", "BARQR_Affected", "BARQR_NotAwareBreathe", 
             "BARQR_AttentionMove", "BARQR_StruggleToRelax", "BARQR_BodyTense", 
             "BARQR_Feeling", "BARQR_Digestion", "BARQR_Comfortable", "BARQR_Unpredictable", 
             "BARQR_AttentionBody", "BARQR_Touched")



SPANE_PAcols = c("SPANE_Positive", "SPANE_Good","SPANE_Pleasant",  "SPANE_Happy", 
                 "SPANE_Joyful",  "SPANE_Contented")

SPANE_NAcols = c("SPANE_Negative", 
                 "SPANE_Bad",  "SPANE_Unpleasant", "SPANE_Sad", "SPANE_Afraid", 
                 "SPANE_Angry")  

SPANEcols = c(SPANE_PAcols, SPANE_NAcols)

library(readxl)
library(dplyr)

#Converst label answers to numbers
qData <- qData %>%
  mutate(across(all_of(PHQcols), 
                ~ case_when(
                  . == "Not at all" ~ 0,
                  . == "Several days" ~ 1,
                  . == "More than half the days" ~ 2,
                  . == "Nearly every day" ~ 3,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>% 
  mutate(across(all_of(MAIAcols), 
                ~ case_when(
                  . == "Never" ~ 0,
                  . == "Sometimes" ~ 1,
                  . == "About half the time" ~ 2,
                  . == "Most of the time" ~ 3,
                  . == "Always" ~ 4,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>% 
  mutate(across(all_of(BARQcols), 
                ~ case_when(
                  . == "Completely disagree" ~ 0,
                  . == "Somewhat disagree" ~ 1,
                  . == "Somewhat agree" ~ 2,
                  . == "Completely agree" ~ 3,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>%
  mutate(across(all_of(SPANEcols), 
                ~ case_when(
                  . == "Very rarely or never" ~ 1,
                  . == "Rarely" ~ 2,
                  . == "Sometimes" ~ 3,
                  . == "Often" ~ 4,
                  . == "Very often or always" ~ 5,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>%
  mutate(across(all_of(BIPScols), 
                ~ case_when(
                  . == "Very often" ~ 5,
                  . == "Fairly often" ~ 4,
                  . == "Sometimes" ~ 3,
                  . == "Almost never" ~ 2,
                  . == "Never" ~ 1,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>%
  mutate(across(all_of(AlexithymiaCols), 
                ~ case_when(
                  . == "Strongly agree" ~ 7,
                  . == "Agree" ~ 6,
                  . == "Slightly agree" ~ 5,
                  . == "Neither agree nor disagree" ~ 4,
                  . == "Slightly disagree" ~ 3,
                  . == "Disagree" ~ 2,
                  . == "Strongly disagree" ~ 1,
                  TRUE ~ NA_real_  # handles any unexpected values
                )
  )) %>%
  mutate(across(all_of(MSEScols), as.numeric))


#View(qData)


# Each scale is scored separately with explicit min/max so that reversed items
# are reflected against the correct bounds, not the global data range.
# Without min/max, psych::scoreItems() infers bounds from the entire qData
# frame — here the MSES items (1–7) would inflate the max used to reverse
# MAIA items (0–4), producing impossible values for NotDistract & NotWorry.

## MAIA-2  (0–4 per item; items 4–9 reverse-scored via '-' prefix)
maia_keys <- list(
  MAIA           = c("BriefMAIA_Uncomfortable", "BriefMAIA_Comfortable",
                     "BriefMAIA_Changes", "-BriefMAIA_Distract", "-BriefMAIA_Discomfort",
                     "-BriefMAIA_Unpleasant", "-BriefMAIA_PhysicalPain", "-BriefMAIA_Worry",
                     "-BriefMAIA_DiscomfortPain", "BriefMAIA_MaintainAwareness",
                     "BriefMAIA_ReturnAwareness", "BriefMAIA_RefocusAttention",
                     "BriefMAIA_NoticeBodyDifferent", "BriefMAIA_NoticeBreathingChanges",
                     "BriefMAIA_NoticeBodyChanges", "BriefMAIA_AwarenessBody",
                     "BriefMAIA_BreathTension", "BriefMAIA_Thoughts", "BriefMAIA_Information",
                     "BriefMAIA_Upset", "BriefMAIA_Listen", "BriefMAIA_HomeBody",
                     "BriefMAIA_BodySafe", "BriefMAIA_TrustBody"),
  MAIA_Notice      = c("BriefMAIA_Uncomfortable", "BriefMAIA_Comfortable", "BriefMAIA_Changes"),
  MAIA_NotDistract = c("-BriefMAIA_Distract", "-BriefMAIA_Discomfort", "-BriefMAIA_Unpleasant"),
  MAIA_NotWorry    = c("-BriefMAIA_PhysicalPain", "-BriefMAIA_Worry", "-BriefMAIA_DiscomfortPain"),
  MAIA_AttentionReg = c("BriefMAIA_MaintainAwareness", "BriefMAIA_ReturnAwareness",
                        "BriefMAIA_RefocusAttention"),
  MAIA_EmoAware    = c("BriefMAIA_NoticeBodyDifferent", "BriefMAIA_NoticeBreathingChanges",
                       "BriefMAIA_NoticeBodyChanges"),
  MAIA_SelfReg     = c("BriefMAIA_AwarenessBody", "BriefMAIA_BreathTension", "BriefMAIA_Thoughts"),
  MAIA_BodyListen  = c("BriefMAIA_Information", "BriefMAIA_Upset", "BriefMAIA_Listen"),
  MAIA_Trust       = c("BriefMAIA_HomeBody", "BriefMAIA_BodySafe", "BriefMAIA_TrustBody")
)
# min=0, max=4: reversal = 4+0-x = 4-x; subscale range 0–12; total range 0–96
maia_scores <- scoreItems(maia_keys, qData[, MAIAcols], totals = TRUE, min = 0, max = 4)

## PHQ-4  (0–3; no reversal needed — just sum pairs directly)
# Anxiety = items 1+2;  Depression = items 3+4
phq_scores <- data.frame(
  Anxiety    = as.numeric(qData$PHQ4_Nervous) + as.numeric(qData$PHQ4_Worry),
  Depression = as.numeric(qData$PHQ4_Depressed) + as.numeric(qData$PHQ4_LittleInterest)
)

## BIPS  (1–5; item 8 "-BIPS_OnTop" reverse-scored)
# Reversal = 5+1-x = 6-x;  sum range: 9 items × 1–5 = 9–45 (item 8 reversed: 1–5 → 1–5)
bips_keys <- list(
  Stress = c("BIPS_Things", "BIPS_Hurry", "BIPS_Pressure",
             "BIPS_Conflict", "BIPS_Hadto", "BIPS_Criticized",
             "BIPS_Difficulties", "-BIPS_OnTop", "BIPS_Worries")
)
bips_scores <- scoreItems(bips_keys, qData[, BIPScols], totals = TRUE, min = 1, max = 5)

## BARQ  (0–3; no reversal)
barq_keys <- list(BARQ = BARQcols)
barq_scores <- scoreItems(barq_keys, qData[, BARQcols], totals = TRUE, min = 0, max = 3)

## SPANE  (1–5; no reversal — positive and negative items scored separately)
spane_keys <- list(Pos = SPANE_PAcols, Neg = SPANE_NAcols)
spane_scores <- scoreItems(spane_keys, qData[, SPANEcols], totals = TRUE, min = 1, max = 5)

## MSES-12  (1–7; 10 items reverse-scored via '-' prefix)
# Reversal = 7+1-x = 8-x; composite and self-doubt subscale
mses_keys <- list(
  MSES = c("MSES12_GoodJob", "MSES12_Convinced",
            "-MSES12_Doubt", "-MSES12_Unhappy",
            "-MSES12_SelfConscious", "-MSES12_Uncomfortable",
            "-MSES12_Loser", "-MSES12_Criticize",
            "-MSES12_Ashamed", "-MSES12_Attractive",
            "-MSES12_SportsPerformance", "-MSES12_SportsNervous"),
  MSES_selfdoubt = c("MSES12_Doubt", "MSES12_Unhappy")  # raw, not reversed
)
mses_scores <- scoreItems(mses_keys, qData[, MSEScols], totals = TRUE, min = 1, max = 7)

## Alexithymia / PAQS-S  (1–7; no reversal)
alex_keys <- list(Alexithymia = AlexithymiaCols)
alex_scores <- scoreItems(alex_keys, qData[, AlexithymiaCols], totals = TRUE, min = 1, max = 7)

## Combine all scored scales into a single data frame
pData <- data.frame(
  maia_scores$scores,
  phq_scores,
  bips_scores$scores,
  barq_scores$scores,
  spane_scores$scores,
  mses_scores$scores,
  alex_scores$scores
)

# Sanity check: print observed ranges for all scales
message("=== Scored scale ranges (should all be within theoretical bounds) ===")
for (col in names(pData)) {
  vals <- pData[[col]][!is.na(pData[[col]])]
  message(sprintf("  %-20s min=%5.1f  max=%5.1f  n=%d", col, min(vals), max(vals), length(vals)))
}
#pData$SWLS = as.numeric(qData$SWLS) 
pData$id = factor(qData$sonaID_1)
names(qData)

##Demographics ------------------
###Gender ------------------------
sum(!is.na(qData$Gender_Male)) #42
sum(!is.na(qData$Gender_Female)) #179
sum(!is.na(qData$Gender_NB)) #4
sum(!is.na(qData$Gender_Transgender)) #0
sum(!is.na(qData$Gender_Other)) #2
sum(!is.na(qData$Gender_PreferNotToAnswer)) #0

N = (42+179+4+2)
42/N *100
179/N *100
4/N *100
2/N *100

qData$Gender <- ifelse(!is.na(qData$Gender_Female), "Female",
                       ifelse(!is.na(qData$Gender_Male), "Male",
                              ifelse(!is.na(qData$Gender_NB) | !is.na(qData$Gender_Other), "Other", NA)))
pData$Gender = qData$Gender
#Age ------------------------------------------
pData$Age = as.numeric(qData$Age)
pData$Age


write.csv(pData,file = questionnaireFile, row.names = FALSE)

