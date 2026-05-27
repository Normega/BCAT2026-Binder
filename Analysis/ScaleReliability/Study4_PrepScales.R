winpath = file.path("I:/Shared drives/Interoception 2025/Paper/SourceScales")
s4scales = file.path(winpath, 'Study4SourceScales.csv')

# Set Up---------
## Load packages ---------
packages <- c("readxl", #read from excel file
              "ggplot2","ggeffects","ggforce","GGally","ggpubr", #graphing stuff
              "corrplot", #graphing correlations
              "psych", #factor analysis / data reduction / questionnaire scoring
              "lme4", "lmerTest", #testing multilevel models
              "DataCombine", 
              "sjPlot", #Making awesome results tables 
              "tidyverse", "tidyr", #reorganizing data
              "stringr", #manipulating strings / letter arrays
              "BayesFactor"
)

#Check if packages above are installed, if not, install them
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

#Load all the packages into memory
options(readr.num_columns = 0)
for (thispack in packages) {
  library(thispack,character.only=TRUE,quietly=TRUE,verbose=FALSE)
}


#Questionnaire Stuff --------------------

# Load individual difference data ------------------
qData = read.csv(s4scales)
qData = qData[-1,]
qData = qData[qData$Finished == "True",]

#Get rid of test id
qData = qData[!qData$SONAid %in% 3055,]
qData = qData[!qData$SONAid %in% 987654321,]
qData = qData[!is.na(qData$Q7),]
qData = qData[order(qData$SONAid, qData$StartDate),]


cutlist = rep(0,nrow(qData))
for (thisRow in 2:nrow(qData)){
  thisId = qData$SONAid[thisRow ]
  lastId = qData$SONAid[thisRow -1]
  if ( thisId  == lastId ){
    cutlist[thisRow] = 1  
  } 
}
qData = qData[!cutlist,]



#recode string answers into numbers
# Assuming df is your dataframe
MAIAcols = c("MAIA24_1", "MAIA24_2", "MAIA24_3", "MAIA24_4", 
             "MAIA24_5", "MAIA24_6", "MAIA24_7", "MAIA24_8", 
             "MAIA24_9", "MAIA24_10", "MAIA24_11", "MAIA24_12", 
             "MAIA24_13", "MAIA24_14", "MAIA24_15", "MAIA24_16", 
             "MAIA24_17", "MAIA24_18", "MAIA24_19", "MAIA24_20", 
             "MAIA24_21", "MAIA24_22", "MAIA24_23", "MAIA24_24"
)

PHQcols = c("PHQ4_1", "PHQ4_2", "PHQ4_3", "PHQ4_4")
BARQcols = c("BARQ_1", "BARQ_2", "BARQ_3", "BARQ_4", 
             "BARQ_5", "BARQ_6", "BARQ_7", "BARQ_8", 
             "BARQ_9", "BARQ_10", "BARQ_11", "BARQ_12")


SPANE_PAcols = c("SPANE_1",
                 "SPANE_3",
                 "SPANE_5",
                 "SPANE_7",
                 "SPANE_10",
                 "SPANE_12")
SPANE_NAcols = c("SPANE_2",
                 "SPANE_4",
                 "SPANE_6",
                 "SPANE_8",
                 "SPANE_9",           
                 "SPANE_11")

PWBcols = c("PWB_1",
            "PWB_2",
            "PWB_3",
            "PWB_4",
            "PWB_5",
            "PWB_6",
            "PWB_7",
            "PWB_8")

SPANEcols = c(SPANE_PAcols, SPANE_NAcols)


BIPScols = c("BIPS_1",
             "BIPS_2",
             "BIPS_3",
             "BIPS_4",
             "BIPS_5",
             "BIPS_6",
             "BIPS_7",
             "BIPS_8",
             "BIPS_9"
)

bqData = qData

######################################################
qData = bqData

qData <- qData |>
  dplyr::mutate(dplyr::across(dplyr::all_of(PHQcols),
    ~ dplyr::case_when(
      . == "Not at all"              ~ 0,
      . == "Several days"            ~ 1,
      . == "More than half the days" ~ 2,
      . == "Nearly everyday"         ~ 3,
      TRUE                           ~ NA_real_
    ))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(MAIAcols),
    ~ dplyr::case_when(
      . == "Never"               ~ 0,
      . == "Sometimes"           ~ 1,
      . == "About half the time" ~ 2,
      . == "Most of the time"    ~ 3,
      . == "Always"              ~ 4,
      TRUE                       ~ NA_real_
    ))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(BARQcols),
    ~ dplyr::case_when(
      . == "Completely disagree" ~ 0,
      . == "Somewhat disagree"   ~ 1,
      . == "Somewhat agree"      ~ 2,
      . == "Completely agree"    ~ 3,
      TRUE                       ~ NA_real_
    ))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(SPANEcols),
    ~ dplyr::case_when(
      . == "Very rarely or never" ~ 1,
      . == "Rarely"               ~ 2,
      . == "Sometimes"            ~ 3,
      . == "Often"                ~ 4,
      . == "Very often or always" ~ 5,
      TRUE                        ~ NA_real_
    ))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(PWBcols),
    ~ dplyr::case_when(
      . == "Agree"                      ~ 7,
      . == "Somewhat agree"             ~ 6,
      . == "Somewhat disagree"          ~ 5,
      . == "Neither agree nor disagree" ~ 4,
      . == "Strongly agree"             ~ 3,
      . == "Disagree"                   ~ 2,
      . == "Strongly disagree"          ~ 1,
      TRUE                              ~ NA_real_
    ))) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(BIPScols),
    ~ dplyr::case_when(
      . == "Very Often"   ~ 5,
      . == "Fairly Often" ~ 4,
      . == "Sometimes"    ~ 3,
      . == "Almost never" ~ 2,
      . == "Never"        ~ 1,
      TRUE                ~ NA_real_
    )))


#View(qData)

# ── Scoring keys ──────────────────────────────────────────────────────────
# Scales with reversed items are scored in separate scoreItems calls with
# explicit min/max to prevent global auto-detection corrupting reversals.
# Passing all columns to one call would use the global min/max (0-7 from
# the PWB scale), reversing MAIA items as 7-x instead of 4-x and BIPS
# item 8 as 7-x instead of 6-x.
#
#   maia_scales  : MAIA total + subscales (0-4; items 4-9 reversed as 4-x)
#   bips_scales  : BIPS (1-5; item 8 reversed as 6-x)
#   other_scales : PHQ-4, BARQ-R, SPANE, Wellbeing (no reversed items)

maia_keys <- list(
  MAIA = c(
    "MAIA24_1",  "MAIA24_2",  "MAIA24_3",
    "-MAIA24_4", "-MAIA24_5", "-MAIA24_6",
    "-MAIA24_7", "-MAIA24_8", "-MAIA24_9",
    "MAIA24_10", "MAIA24_11", "MAIA24_12",
    "MAIA24_13", "MAIA24_14", "MAIA24_15",
    "MAIA24_16", "MAIA24_17", "MAIA24_18",
    "MAIA24_19", "MAIA24_20", "MAIA24_21",
    "MAIA24_22", "MAIA24_23", "MAIA24_24"
  ),
  MAIA_Noticing     = c("MAIA24_1",  "MAIA24_2",  "MAIA24_3"),
  MAIA_NotDistract  = c("-MAIA24_4", "-MAIA24_5",  "-MAIA24_6"),
  MAIA_NotWorry     = c("-MAIA24_7", "-MAIA24_8",  "-MAIA24_9"),
  MAIA_AttentionReg = c("MAIA24_10", "MAIA24_11", "MAIA24_12"),
  MAIA_EmoAware     = c("MAIA24_13", "MAIA24_14", "MAIA24_15"),
  MAIA_SelfReg      = c("MAIA24_16", "MAIA24_17", "MAIA24_18"),
  MAIA_BodyListen   = c("MAIA24_19", "MAIA24_20", "MAIA24_21"),
  MAIA_Trust        = c("MAIA24_22", "MAIA24_23", "MAIA24_24")
)

bips_keys <- list(
  Stress = c("BIPS_1", "BIPS_2", "BIPS_3", "BIPS_4",
             "BIPS_5", "BIPS_6", "BIPS_7", "-BIPS_8", "BIPS_9")
)

other_keys <- list(
  Anxiety    = c("PHQ4_1", "PHQ4_2"),
  Depression = c("PHQ4_3", "PHQ4_4"),
  BARQ       = c("BARQ_1",  "BARQ_2",  "BARQ_3",  "BARQ_4",
                 "BARQ_5",  "BARQ_6",  "BARQ_7",  "BARQ_8",
                 "BARQ_9",  "BARQ_10", "BARQ_11", "BARQ_12"),
  Pos        = SPANE_PAcols,
  Neg        = SPANE_NAcols,
  Wellbeing  = PWBcols
)

maia_scales  <- psych::scoreItems(maia_keys,  qData[, MAIAcols], min = 0, max = 4)
bips_scales  <- psych::scoreItems(bips_keys,  qData[, BIPScols], min = 1, max = 5)
other_scales <- psych::scoreItems(other_keys, qData[, c(PHQcols, BARQcols, SPANEcols, PWBcols)])

pData        <- data.frame(maia_scales$scores,
                           bips_scales$scores,
                           other_scales$scores)
pData$SES    <- as.numeric(qData$Q12_1)
pData$id     <- factor(qData$SONAid)
names(pData)

# ── Scale Reliability and Descriptive Statistics ──────────────────────────
message("\n=== Study 4: Scale Reliability and Descriptive Statistics ===")

resultsPath_s4 <- winpath

# Helper: extract alpha + descriptives from scoreItems alpha matrix
.scale_rel_s4 <- function(score_obj, pdata_df) {
  alphas <- score_obj$alpha   # matrix: row = "alpha", cols = scale names
  dplyr::bind_rows(lapply(colnames(alphas), function(sub) {
    col  <- pdata_df[[sub]]
    vals <- col[!is.na(col)]
    tibble::tibble(
      subscale = sub,
      alpha    = round(as.numeric(alphas["alpha", sub]), 3),
      n        = length(vals),
      M        = round(mean(vals), 2),
      SD       = round(sd(vals),   2),
      min_obs  = round(min(vals),  2),
      max_obs  = round(max(vals),  2)
    )
  }))
}

.interitem_r <- function(x, y) {
  round(cor(as.numeric(x), as.numeric(y), use = "complete.obs"), 3)
}

table_reliability_s4 <- dplyr::bind_rows(
  .scale_rel_s4(maia_scales,  pData),
  .scale_rel_s4(bips_scales,  pData),
  .scale_rel_s4(other_scales, pData)
) |>
  dplyr::mutate(
    scale = dplyr::case_when(
      subscale == "MAIA" | grepl("^MAIA_", subscale) ~ "MAIA-2",
      subscale %in% c("Anxiety", "Depression")        ~ "PHQ-4",
      subscale == "BARQ"                              ~ "BARQ-R",
      subscale %in% c("Pos", "Neg")                  ~ "SPANE",
      subscale == "Wellbeing"                         ~ "Wellbeing",
      subscale == "Stress"                            ~ "BIPS",
      TRUE                                            ~ subscale
    ),
    n_items = dplyr::case_when(
      subscale == "MAIA"                              ~ 24L,
      grepl("^MAIA_", subscale)                       ~ 3L,
      subscale == "BARQ"                              ~ 12L,
      subscale %in% c("Pos", "Neg")                  ~ 6L,
      subscale == "Wellbeing"                         ~ 8L,
      subscale == "Stress"                            ~ 9L,
      subscale %in% c("Anxiety", "Depression")        ~ 2L,
      TRUE                                            ~ NA_integer_
    ),
    interitem_r = dplyr::case_when(
      subscale == "Anxiety"    ~ .interitem_r(qData$PHQ4_1, qData$PHQ4_2),
      subscale == "Depression" ~ .interitem_r(qData$PHQ4_3, qData$PHQ4_4),
      TRUE                     ~ NA_real_
    ),
    alpha = dplyr::if_else(n_items == 2L, NA_real_, alpha),
    note = dplyr::case_when(
      n_items == 2L ~ "2-item subscale; alpha not reported, see interitem_r",
      n_items == 3L ~ "3-item subscale; alpha interpreted with caution",
      TRUE          ~ ""
    )
  ) |>
  dplyr::select(scale, subscale, n_items, n, M, SD,
                min_obs, max_obs, alpha, interitem_r, note)

cat("\n--- Study 4 reliability summary ---\n")
print(
  table_reliability_s4 |>
    dplyr::select(scale, subscale, n_items, alpha, interitem_r, n) |>
    as.data.frame(),
  row.names = FALSE
)

readr::write_csv(
  table_reliability_s4,
  file.path(resultsPath_s4, "table_scale_reliability_study4.csv")
)
message(sprintf(
  "Saved: table_scale_reliability_study4.csv (%d rows)",
  nrow(table_reliability_s4)
))
