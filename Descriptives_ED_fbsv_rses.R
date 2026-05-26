# Load library
library(dplyr)
library(Hmisc)
library(moments)
library(EnvStats)
library(car)

# Read in data
EDData <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses.rds")

# Calculating descriptives on continuous data
mean(EDData$age, na.rm=TRUE)
sd(EDData$age, na.rm=TRUE)

mean(EDData$edeqs_total)
sd(EDData$edeqs_total)
range(EDData$edeqs_total)

mean(EDData$fbsv_total)
sd(EDData$fbsv_total)

mean(EDData$ptq_total)
sd(EDData$ptq_total)

mean(EDData$rses_sc_score)
sd(EDData$rses_sc_score)

# Percentage of the sample that meets the clinical cutoff on the edeqs
clinical_edeqs <- sum(EDData$edeqs_total > 15, na.rm = TRUE)
clinical_perc <- (clinical_edeqs / nrow(EDData)) * 100

# Calculating descriptives on categorical variables
count(EDData,sex) 
count(EDData, race_ethnicity)
count(EDData, household_income)
count(EDData, education)

# Percentages
# Male
(617/nrow(EDData)) * 100
# Female
(1514/nrow(EDData)) * 100

# Race_ethnicity: American Indian or Alaskan Native
(9/nrow(EDData))*100
# Race_ethnicity: Asian
(580/nrow(EDData)) * 100 
# Race_ethnicity: Black_AfAm
(240/nrow(EDData)) * 100
# Race_ethnicity: Hispanic
(228/nrow(EDData)) * 100
# Race_ethnicity: Native Hawaiian or Other Pacific Islander
(8/nrow(EDData)) * 100
# Race_ethnicity: White
(946/nrow(EDData)) * 100
# Race_ethnicity: Multiracial
(109/nrow(EDData)) * 100
# Race_ethnicity: Other
(11/nrow(EDData)) * 100

# Income: $31000 or less
(144/nrow(EDData)) * 100
# Income: $31,001-$42,000
(177/nrow(EDData)) * 100
# Income: $42,001-$126,000
(701/nrow(EDData)) * 100
# Income: $126,001-$188,000
(405/nrow(EDData)) * 100
# Income: $188,001 or more
(704/nrow(EDData)) * 100

# Education: 0=Currently enrolled in elementary school
(0/nrow(EDData)) * 100
# Education: 1=Currently enrolled in middle/junior high school
(0/nrow(EDData)) * 100
# Education 2=Currently enrolled in high school
(12/nrow(EDData)) * 100
# Education 3=Didn't finish high school
(0/nrow(EDData)) * 100
# Education 4=Didn't finish high school, but completed a technical/ vocational program
(1/nrow(EDData)) * 100
# Education 5=High school graduate or GED (General Education Diploma)
(318/nrow(EDData)) * 100
# Education 6=Completed high school and a technical/vocational program
(18/nrow(EDData)) * 100
# Education 7=Less than 2 years of college
(1033/nrow(EDData)) * 100
# Education 8=2 years of college or more/ including associate degree or equivalent
(376/nrow(EDData)) * 100
# Education 9=College graduate (4 or 5 year program)
(329/nrow(EDData)) * 100
# Education 10=Master's degree (or other post-graduate training)
(31/nrow(EDData)) * 100
# Education 11=Doctoral degree (PhD, MD, EdD, DVM, DDS, JD, etc)
(13/nrow(EDData)) * 100

# Correlation matrix
# Trim to just the variables that should be in the correlation matrix
EDData_subset <- EDData[c(grep("edeqs_total|fbsv_total|ptq_total|rses_sc_score",names(EDData)))]
names(EDData_subset)

# Run correlation matrix
EDData_corr <- rcorr(as.matrix(EDData_subset))
EDData_corr

# Check individual correlation p-values
cor1 <- cor.test(EDData$edeqs_total, EDData$fbsv_total, method = "pearson")
cor1

cor2 <- cor.test(EDData$edeqs_total, EDData$ptq_total, method = "pearson")
cor2

cor3 <- cor.test(EDData$edeqs_total, EDData$rses_sc_score, method = "pearson")
cor3

cor4 <- cor.test(EDData$fbsv_total, EDData$ptq_total, method = "pearson")
cor4

cor5 <- cor.test(EDData$fbsv_total, EDData$rses_sc_score, method = "pearson")
cor5

cor6 <- cor.test(EDData$ptq_total, EDData$rses_sc_score, method = "pearson")
cor6

# Check histograms, boxplots, density plots
# Histograms
hist(EDData$edeqs_total,
     xlab = "Eating Disorder Symptoms",
     breaks = sqrt(length(EDData$edeqs_total)) # set number of bins
)

hist(EDData$fbsv_total,
     xlab = "Bullying",
     breaks = sqrt(length(EDData$fbsv_total)) # set number of bins
)

hist(EDData$ptq_total,
     xlab = "Rumination",
     breaks = sqrt(length(EDData$ptq_total)) # set number of bins
)

hist(EDData$rses_sc_score,
     xlab = "Self-esteem",
     breaks = sqrt(length(EDData$rses_sc_score)) # set number of bins
)

hist(EDData$age,
     xlab = "Age",
     breaks = sqrt(length(EDData$age)) # set number of bins
)

# Boxplots
boxplot(EDData$edeqs_total,
        ylab = "Eating Disorder Symptoms"
)

boxplot(EDData$fbsv_total,
        ylab = "Bullying"
)

boxplot(EDData$ptq_total,
        ylab = "Rumination"
)

boxplot(EDData$rses_sc_score,
        ylab = "Self Esteem"
)

# Density plots
plot(density(EDData$edeqs_total), main="Density Plot of Eating Disorder Symptoms", xlab="Eating Disorder Symptoms")

plot(density(EDData$fbsv_total), main="Density Plot of Bullying Experiences", xlab="Bullying Experiences")

plot(density(EDData$ptq_total), main="Density Plot of Rumination", xlab="Rumination")

plot(density(EDData$rses_sc_score), main="Density Plot of Self-esteem", xlab="Self-esteem")

# Calculate Skewness
skewness(EDData$edeqs_total)

skewness(EDData$fbsv_total)

skewness(EDData$ptq_total)

skewness(EDData$rses_sc_score)

# Calculate Kurtosis
kurtosis(EDData$edeqs_total)

kurtosis(EDData$fbsv_total)

kurtosis(EDData$ptq_total)

kurtosis(EDData$rses_sc_score)

# Rosner's tests for outliers
rosner_result <- rosnerTest(EDData$edeqs_total, k = 5)
print(rosner_result)

rosner_result <- rosnerTest(EDData$fbsv_total, k = 5)
print(rosner_result)

rosner_result <- rosnerTest(EDData$ptq_total, k = 5)
print(rosner_result)

rosner_result <- rosnerTest(EDData$rses_sc_score, k = 5)
print(rosner_result)

# Calculate the variance inflation factor (VIF)
model <- lm(edeqs_total ~ fbsv_total + ptq_total + rses_sc_score, data = EDData)
vif(model)

# Cronbach's Alpha
# Function to extract Cronbach's alpha
get_alpha <- function(data_items) {
  result <- alpha(data_items, check.keys = FALSE, use = "complete.obs")
  return(result$total$raw_alpha)
}

# EDEQ
edeq_items <- EDData[c("edeqs_1", "edeqs_2", "edeqs_3", "edeqs_4",
                       "edeqs_5", "edeqs_6", "edeqs_7", "edeqs_8",
                       "edeqs_9", "edeqs_10", "edeqs_11", "edeqs_12")]
edeq_alpha <- get_alpha(edeq_items)

# FBSV
fbsv_items <- EDData[c("fbsv_1", "fbsv_2", "fbsv_3", "fbsv_4",
                       "fbsv_5", "fbsv_6", "fbsv_7", "fbsv_8",
                       "fbsv_9", "fbsv_10")]
fbsv_alpha <- get_alpha(fbsv_items)

# RSES
rses_items <- EDData[c("rses_sc_1", "rses_sc_2", "rses_sc_3", "rses_sc_4",
                       "rses_sc_5", "rses_sc_6", "rses_sc_7", "rses_sc_8",
                       "rses_sc_9", "rses_sc_10")]
rses_alpha <- get_alpha(rses_items)

# PTQ
ptq_items <- EDData[c("ptq_1", "ptq_2", "ptq_3", "ptq_4", "ptq_5",
                      "ptq_6", "ptq_7", "ptq_8", "ptq_9", "ptq_10",
                      "ptq_11", "ptq_12", "ptq_13", "ptq_14", "ptq_15")]
ptq_alpha <- get_alpha(ptq_items)

# Print results
cat("Cronbach's Alpha Values:\n")
cat("EDEQ:", edeq_alpha, "\n")
cat("FBSV:", fbsv_alpha, "\n")
cat("RSES:", rses_alpha, "\n")
cat("PTQ:", ptq_alpha, "\n")

# Bullying Skew Diagnostics
summary(EDData$fbsv_total)
mean(EDData$fbsv_total == 0, na.rm = TRUE)
mean(EDData$fbsv_total == min(EDData$fbsv_total, na.rm = TRUE), na.rm = TRUE)
hist(EDData$fbsv_total)


