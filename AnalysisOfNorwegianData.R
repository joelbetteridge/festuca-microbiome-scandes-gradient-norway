# The effect of environmental conditions on fungal/bacterial biodiversity ------


# Experimental overview --------------------------------------------------------
#
# This study investigated the effects of different temperature and precipitation
# averages on the biodiversity of bacteria and fungi that live alongside Festuca
# rubra.
#
# Data was collected in a factorial design of 12 sampling sites across the
# scandes mountains. Festuca rubra grass roots were collected and
# metagenomically sequenced to see associated fungi and bacteria.


# Description of data ----------------------------------------------------------
#
# All data is stored in the data-raw folder.
#
# The numerical data are named Tidy_Norway_ASVS... These are .txt files where
# columns represent different environmental samples and rows represent unique 
# ASVs. Values are the counts of each ASV in each sample.
#
# The character data are named Tidy_Norway_taxonomy... These are .txt files
# where the columns show the phylogentic classification and the rows contain
# each unique ASV
#
# The location and environmental data is stored in Tidy_sampleData...txt file.


# Analysis overview ------------------------------------------------------------
#
# This script imports, explores and analyses the relationship and significance
# between bacteral/fungal species and precipitation/temperature. PCA was used to
# observe patterns followed by shannon index calculatons and a two-way ANOVA.
# Figures are created for these results, some of which are in the report.



# Packages required ------------------------------------------------------------
# ------------------------------------------------------------------------------

# for import, manipulation and plotting 
library(tidyverse)

# for shannon index
library(vegan)

# for PCA visualisation
library(factoextra)

# Data import ------------------------------------------------------------------
# ------------------------------------------------------------------------------

about <- read_table(
  "data-raw/Tidy_sampleData_Norway_root_endophytes_York.txt",
                    col_names = TRUE)

from16S <- read_table(
  "data-raw/Tidy_Norway_ASVs_16S_root_endophytes_200_York.txt",
                   col_names = TRUE)

fromITS <- read_table(
  "data-raw/Tidy_Norway_ASVs_ITS_root_endophytes_200_York.txt",
                      col_names = TRUE)

taxonomy16S <- read_table(
  "data-raw/Tidy_Norway_taxonomy_16S_York.txt",
                          col_names = TRUE)

taxonomyITS <- read_table(
  "data-raw/Tidy_Norway_taxonomy_ITS_York.txt",
                          col_names = TRUE)

# Data exploration -------------------------------------------------------------
# ------------------------------------------------------------------------------

# Plot the predicted summer temperatures & precipitation values for all sites
ggplot(data = about, aes(x = temperature, y = precipitation)) +
  geom_point() +
  xlab("Expected Temperature (°C)") +
  ylab("Expected Precipitation (mm)")

# Plot the actual average annual temperature and precipitation for all sites
ggplot(data = about, aes(x = T.grid, y = P.grid)) +
  geom_point() +
  xlab("Average Summer Temperature (°C)") +
  ylab("Average Summer Precipitation (mm)")

# Interpretation: Each temperature is paired with each precipitation level
# across sites. This confirms the factorial design is balanced across treatment,
# making it suitable for ANOVA.

# Group ASVs by taxonomy -------------------------------------------------------
# ------------------------------------------------------------------------------

vec <- as.vector(from16S$ASV) # extract ASVs from 'from16S' into a vector
taxonomy16S <- taxonomy16S %>%
  filter(ASV %in% vec) # extract ASVs from 'taxonomy16S' which are in 'from16S'

# Create a function for grouping
groupBy <- function(countTable, taxonomicTable, taxonomicGrouping){
  
  # find the unique taxonomic groups in the taxonomic table
  groups<- unique(taxonomicTable[[taxonomicGrouping]])
  groups<- groups[-which(is.na(groups))]
  # loop through each group:
  out<- sapply(groups, function(group){
    
    # find the rows in countTable that correspond to Asvs that are part of group
    ids<- taxonomicTable$ASV[which(taxonomicTable[[taxonomicGrouping]]==group)]
    
    # extract the relevant rows from countTable and sum of all the those ASVs
    idNums=which(countTable[[1]] %in% ids)
    colSums(countTable[idNums,-1], na.rm=T)
    
  })
  
  # remove columns that are all empty
  out=out[,which(colSums(out, na.rm=T)>0)]
  
  out
}

# Apply the function to group the ASV data by phylum
phylum16S <- groupBy(from16S[, -dim(from16S)[2]], taxonomy16S, "phylum")
phylumITS <- groupBy(fromITS[, -dim(fromITS)[2]], taxonomyITS, "phylum")


# Fungal PCA -------------------------------------------------------------------
# PCA used to see whether samples from similar environments group together based
# on their fungal communities. A useful first look before doing formal
# statistics as it can show patterns that might not be obvious otherwise.

# Run PCA on fungal phylum data
pca_ITS <- prcomp(phylumITS, scale = TRUE)

# Match environmental data to samples
ids <- rownames(phylumITS)
idNums <- which(about$site %in% ids)

# Add environmental variables as a factor for colouring
# Using precipitation as discrete groups
precip_groups <- as.factor(about$precipitation[idNums])
temp_groups <- as.factor(about$temperature[idNums])

# Plot coloured by precipitation
fviz_pca_ind(pca_ITS,
             col.ind = precip_groups,
             legend.title = "Precipitation (mm)",
             label = "none",
             pointshape = 19,
             mean.point = FALSE) +
  theme_classic() +
  labs(title = "PCA of Fungal Community Composition")

# Interpretation : no clear grouping by precipitation . Fungal
# communities at phylum level are mixed across all sites, which infers the
# precipitation differences are not strong enough to separate communities
# at this broad taxonomic level.

# Plot coloured by temperature
fviz_pca_ind(pca_ITS,
             col.ind = temp_groups,
             legend.title = "Temperature (°C)",
             label = "none",
             pointshape = 19,
             mean.point = FALSE) +
  theme_classic() +
  labs(title = "PCA of Fungal Community Composition")

# Interpretation : no clear grouping by temperature . Fungal
# communities at phylum level are mixed across all sites, which infers the
# temperature differences are not strong enough to separate communities
# at this broad taxonomic level.

# Bacterial PCA ----------------------------------------------------------------
# ------------------------------------------------------------------------------

pca_16S <- prcomp(phylum16S, scale = TRUE)

ids <- rownames(phylum16S)
idNums <- which(about$site %in% ids)

precip_groups <- as.factor(about$precipitation[idNums])
temp_groups <- as.factor(about$temperature[idNums])

# Plot coloured by precipitation
fviz_pca_ind(pca_16S,
             col.ind = precip_groups,
             legend.title = "Precipitation (mm)",
             label = "none",
             pointshape = 19,
             mean.point = FALSE) +
  theme_classic() +
  labs(title = "PCA of Bacterial Community Composition")

# Interpretation : Some of the lower precipitation values (pink + green) appear
# higher on the y axis. However, Dim1 explains 21.1% and Dim2 explains 14.2% of
# the variation, meaning together they only capture about a third of the total
# variation in community composition.

# Plot coloured by temperature
fviz_pca_ind(pca_16S,
             col.ind = temp_groups,
             legend.title = "Temperature (°C)",
             label = "none",
             pointshape = 19,
             mean.point = FALSE) +
  theme_classic() +
  labs(title = "PCA of Bacterial Community Composition")

# Interpretation : Most high temperature values (blue) appear above 0,0
# suggesting temperature beinga driver of community composition. However, as 
# stated above only around a third of the variation is displayed.

pca_scores <- as.data.frame(pca_16S$x)
pca_scores$temperature <- temp_groups
pca_scores$precipitation <- precip_groups

# Putting temperature and precipitation on the same plot (colour and shape) so
# we can see both variables at once and spot any patterns where the two interact
pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, 
                       colour = temperature,
                       shape = precipitation)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c(
    "6.5"  = "#74ADD1",
    "8.5"  = "#FEE090",
    "10.5" = "#F46D43"
  )) +
  scale_shape_manual(values = c(
    "600"  = 15,
    "1200" = 17,
    "2000" = 18,
    "2700" = 19
  )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = paste0("Dim1 (", round(summary(pca_16S)$importance[2,1]*100, 1), "%)"),
    y = paste0("Dim2 (", round(summary(pca_16S)$importance[2,2]*100, 1), "%)"),
    colour = "Temperature (°C)",
    shape = "Precipitation (mm)"
  )

pca

# Interpretation : High temperature , low precipitation values appear much 
# higher on the y axis compared to low temperature , any precipitation.
# Potential relationship between temperature and precipitation for bacterial 
# communities.

# Save it to file
ggsave("./figures/pca.tif", 
       plot = pca, 
       device = "tiff",
       width = 6, 
       height = 6,
       units = "in",
       dpi = 300)


# Assess correlations between temperature/precipitation on fungal biodiversity -
# ------------------------------------------------------------------------------

# Calculate shannon index of fungal diversity
fungal_diversity <- diversity(phylumITS, index = "shannon")

# Extract sample IDs
ids <- rownames(phylumITS)

# Identify rows in the data that match the sample IDs
idNums <- which(about$site %in% ids)

# Create a dataframe combining fungal diversity with environmental variables
fungal_env <- data.frame(
  diversity = fungal_diversity,
  precip = about[idNums, "precipitation"],
  temp = about[idNums, "temperature"]
)

# Creating a grouped boxplot of fungal shannon diversity and expected
# environmental variables

fungal_diversity <- ggplot(fungal_env, aes(x = as.factor(precipitation), y = diversity,
                       fill = as.factor(temperature))) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.75), 
    alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("#74ADD1", "#FEE090", "#F46D43")) +
  theme_classic() +
  labs(x = "Expected Summer Precipitation (mm)" ,
       y = "Fungal Shannon Diversity Index") +
  guides(fill = guide_legend(
    title = "Expected Summer\nTemperature (°C)"
  ))

fungal_diversity

# Save it to file
ggsave("./figures/fungal_diversity.tif", 
       plot = fungal_diversity, 
       device = "tiff",
       width = 6, 
       height = 6,
       units = "in",
       dpi = 300)

# Two-way ANOVA ----------------------------------------------------------------
# We want to test temperature and precipitation both separately and together.
# I use as factors rather than continuous variables because the experiment used
# set discrete treatment levels.

# Compute the analysis of variance
fungal.aov <- aov(diversity ~ as.factor(precipitation) * as.factor(temperature),
                  data = fungal_env)

# Summary of the analysis
summary(fungal.aov)

#                                    Df Sum Sq Mean Sq F value  Pr(>F)   
# as.factor(precipitation)           3  0.774 0.25803   4.726 0.00431 **
#   as.factor(temperature)           2  0.157 0.07850   1.438 0.24333   
# as.factor(precip):as.factor(temp)  6  0.880 0.14664   2.686 0.01976 * 
#   Residuals                                       82  4.477 0.05459                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interpretation : Precipitation alone is significant (0.01) and 
# the relationship between precipitation and temperature is significant (0.05).

# Tukey HSD
# Acts as follow-up comparisons between precipitation groups because
# it accounts for us are making multiple comparisons at once, which would
# otherwise increase the chance of a false positive
TukeyHSD(fungal.aov, "as.factor(precipitation)")
plot(TukeyHSD(fungal.aov, "as.factor(precipitation)"))

# Interpretation :
# No significant difference between:
# 1200 and 600 (p-value = 0.084)
# 2000 and 600 (p-value = 0.992)
# 2000 and 1200 (p-value = 0.162)
# 2700 and 1200 (p-value = 0.906)
#
# Significant difference between:
# 2700 and 600 (p-value = 0.013)
# 2700 and 2000 (p-value = 0.032)


# Assumption testing for ANOVA -------------------------------------------------
# Assumption checks done here to confirm the ANOVA results above are reliable.
# Normality and equal variance are both met.

# Shapiro-Wilk normality test
shapiro.test(residuals(fungal.aov))

# Interpretation : Normality not violated (w = 0.982 , p = 0.231)

# Bartlett test for equal variance
bartlett.test(diversity ~ as.factor(precipitation), data = fungal_env)
bartlett.test(diversity ~ as.factor(temperature), data = fungal_env)

# Interpretation : Homogeneity of variance not violated (k-squared = 1.16 ,
# df = 3 , p = 0.7625)

# Residuals vs fitted
plot(fungal.aov, which = 1) # Interpretation : Acceptable

# Q-Q plot
qqnorm(residuals(fungal.aov))
qqline(residuals(fungal.aov), col="red") # Interpretation : Acceptable

# Histogram of residuals
hist(residuals(fungal.aov)) # Interpretation : Acceptable


# What specific groups of fungi are environmentally associated -----------------
# ------------------------------------------------------------------------------

# Group fungi by class
classITS <- groupBy(fromITS[, -dim(fromITS)[2]], taxonomyITS, "class")

# Extract sample IDs from class16S
ids <- rownames(classITS)

# Match samples to the respective rows 
idNums <- which(about$site %in% ids)

# Create loop
pvals = apply(classITS, 2, function(class){
  
  # create a data frame for the linear regression
  class_lr = data.frame(abundance = class, precip = about[idNums,"P.grid"],
                        temp = about[idNums, "T.grid"])
  
  
  # perform a linear regression, with precipitation and temperature as input
  # variables and class abundance as the output variable.
  out = lm(abundance ~ P.grid + T.grid, data = class_lr)
  
  # extract the p-values of these fits.
  sumres = summary(out)
  pf(sumres$fstatistic[1L],
     sumres$fstatistic[2L],
     sumres$fstatistic[3L],
     lower.tail = FALSE)
  
})

# adjust the p-values for multiple comparisons
p.adjust(pvals)

# Interpretation : No specific groups are significant , suggesting that at this 
# taxonomic level gradient doesn't favour certain groups.


# Assess correlations between temperature and precipitation on bacterial
# biodiversity -----------------------------------------------------------------


# Calculate shannon index of bacterial diversity
bacterial_diversity <- diversity(phylum16S, index = "shannon")

# Extract sample IDs
ids <- rownames(phylum16S)

# Identify rows in the data that match the sample IDs
idNums <- which(about$site %in% ids)

# Create a dataframe combining bacterial diversity with environmental variables
bacterial_env <- data.frame(
  diversity = bacterial_diversity,
  precip = about[idNums, "precipitation"],
  temp = about[idNums, "temperature"]
)


# Creating a grouped boxplot of shannon diversity and expected environmental
# variables

bacterial_diversity <- ggplot(bacterial_env, aes(x = as.factor(precipitation), y = diversity,
                          fill = as.factor(temperature))) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.75), 
    alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("#74ADD1", "#FEE090", "#F46D43")) +
  theme_classic() +
  labs(x = "Expected Summer Precipitation (mm)" ,
       y = "Bacterial Shannon Diversity Index") +
  guides(fill = guide_legend(
    title = "Expected Summer\nTemperature (°C)"
  ))

bacterial_diversity


# Save it to file
ggsave("./figures/bacterial_diversity.tif", 
       plot = bacterial_diversity, 
       device = "tiff",
       width = 6, 
       height = 6,
       units = "in",
       dpi = 300)

# Two-way ANOVA ----------------------------------------------------------------
# Same approach as the fungal ANOVA. A two-way with an interaction term to test
# whether temperature and precipitation affect bacterial diversity individually
# and in combination

# Compute the analysis of variance
bacterial.aov <- aov(diversity ~
                       as.factor(precipitation) * as.factor(temperature),
                       data = bacterial_env)

# Summary of the analysis
summary(bacterial.aov)

#                                   Df Sum Sq Mean Sq F value  Pr(>F)   
# as.factor(precipitation)           3  0.328 0.10927   2.032 0.11578   
# as.factor(temperature)             2  0.079 0.03942   0.733 0.48352   
# as.factor(precip):as.factor(temp)  6  1.202 0.20037   3.726 0.00251 **
#   Residuals                                       82  4.409 0.05377                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interpretation : Relationship between precipitation and temperature is 
# significant (0.01).


# Interaction plot - visualise diversity changes across environmental
# combinations

interaction.plot(bacterial_env$temp,
                 bacterial_env$precip,
                 bacterial_env$diversity,
                 xlab = "Temperature",
                 ylab = "Shannon Diversity",
                 trace.label = "Precipitation")

# Interpretation : Lines cross and behave differently suggesting interaction.


# Assumption testing for bacterial ANOVA ---------------------------------------
# ------------------------------------------------------------------------------

# Shapiro-Wilk normality test
shapiro.test(residuals(bacterial.aov))

# Interpretation : Normality violated (w = 0.862 , p = 7.091e-08)

# Bartlett test for equal variance
bartlett.test(diversity ~ as.factor(precipitation), data = bacterial_env)
bartlett.test(diversity ~ as.factor(temperature), data = bacterial_env)

# Interpretation : Homogeneity normality violated (k-squared = 15.416 , 
# df = 3 , p = 0.0015)

# Residuals vs fitted
plot(bacterial.aov, which = 1) # Interpretaton : Unacceptable deviation

# Q-Q plot
qqnorm(residuals(bacterial.aov))
qqline(residuals(bacterial.aov), col="red") # Interpretation : Unacceptable

# Histogram of residuals
hist(residuals(bacterial.aov)) # Interpretation : Unacceptable distribution


# Kruskall-Wallis test ---------------------------------------------------------
# This test is used instead as it is the non-parametric equivalent of a one-way
# ANOVA and does not assume normality or equal variance.

# Testing effect of precipitation on bacterial diversity
kruskal.test(diversity ~ as.factor(precipitation), data = bacterial_env)

# Testing effect of temperature on bacterial diversity
kruskal.test(diversity ~ as.factor(temperature), data = bacterial_env)

# Interpretation : neither precipitation (p = 0.246) nor temperature (p = 0.730)
# had a significant effect on bacterial diversity individually. This contrasts
# with the interaction term in the original ANOVA (p = 0.003) suggesting that
# temperature and precipitation only affect bacterial diversity when considered
# together rather than in isolation. This pattern is visible in the interaction
# plot above where the lines cross.


# What specific groups of bacteria are environmentally associated --------------
# ------------------------------------------------------------------------------
# Group bacteria by class
class16S <- groupBy(from16S[, -dim(from16S)[2]], taxonomy16S, "class")

# Extract sample IDs from class16S
ids <- rownames(class16S)

# Match samples to the respective rows 
idNums <- which(about$site %in% ids)

# Create loop
pvals = apply(class16S, 2, function(class){
  
  # create a data frame for the linear regression
  class_lr = data.frame(abundance = class,
                        precip = about[idNums,"P.grid"],
                        temp = about[idNums, "T.grid"])
  
  
  # perform a linear regression, with precipitation and temperature as input 
  # variables and class abundance as the output variable.
  out = lm(abundance ~ P.grid + T.grid, data = class_lr)
  
  # extract the p-values of these fits.
  sumres = summary(out)
  pf(sumres$fstatistic[1L],
     sumres$fstatistic[2L],
     sumres$fstatistic[3L],
     lower.tail = FALSE)
  
})

# adjust the p-values for multiple comparisons
p.adjust(pvals)

# Bacilli and Alphaproteobacteria analysed individually as they both had
# significant p-values, suggesting environmental relationships.

# Create data frame for linear regression
bacilli<- class16S[,"Bacilli"] # extract bacilli abundance
alpha<- class16S[,"Alphaproteobacteria"] # extract alphaproteobacteria abundance
df<-data.frame(bacilli=bacilli, alpha=alpha,
               precip=about[idNums,"P.grid"],
               temp=about[idNums, "T.grid"])

# Perform linear regressions for Bacilli and Alphaproteobacteria
outBacilli<-lm(bacilli ~ P.grid + T.grid, data=df)
outAlpha<-lm(alpha ~ P.grid + T.grid, data=df)

summary(outBacilli)

#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -43.051874  31.509914  -1.366  0.175
# P.grid       -0.015332   0.007134  -2.149  0.034 *
# T.grid       10.639622   3.384146   3.144  0.002 **
#
# Residual standard error: 57.96 on 91 degrees of freedom
# Multiple R-squared: 0.1322, Adjusted R-squared: 0.1131
# F-statistic: 6.931 on 2 and 91 DF, p-value: 0.001579

# Interpretation : Bacilli abundance increases with temperature (p = 0.002) and
# decreases slightly with precipitation (p = 0.034). The model is significant
# overall (p = 0.002) but only explains about 13% of the variation in abundance
# suggesting other factors are also important.

summary(outAlpha)

#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) 533.17871   89.35954   5.967 4.58e-08 ***
# P.grid       -0.06529    0.02023  -3.227  0.002 **
# T.grid      -14.88493    9.59716  -1.551  0.124
#
# Residual standard error: 164.4 on 91 degrees of freedom
# Multiple R-squared: 0.128, Adjusted R-squared: 0.1088
# F-statistic: 6.678 on 2 and 91 DF, p-value: 0.001968

# Interpretation : Alphaproteobacteria abundance decreases significantly with
# precipitation (p = 0.002) but temperature had no significant effect
# (p = 0.124). The model is significant overall (p = 0.002) but explains only 
# about 13% of the variation. This reacts oppositely to Bacilli, suggesting 
# different responses to the same conditions.

# Check whether the model predictions are reasonable for Bacilli and Alpha

predicted_bacilli <- data.frame("Predicted" = predict(outBacilli),
                                "Actual" = bacilli)
ggplot(predicted_bacilli, aes(x = Predicted, y = Actual)) +
  geom_point() +
  xlab("Predicted Bacilli Abundance") +
  ylab("Actual Bacilli Abundance") +
  geom_abline(slope = 1, intercept = 0)

# Interpretation : weak trend with lots of scatter

predicted_alpha <- data.frame("Predicted" = predict(outAlpha),
                                "Actual" = alpha)

ggplot(predicted_alpha, aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  xlab("Predicted Alphaproteobacteria Abundance") +
  ylab("Actual Alphaproteobacteria Abundance") +
  geom_abline(slope = 1, intercept = 0)

# Interpretation : better trend but still lots of scatter

              
# Order-level analysis of Alphaproteobacteria ----------------------------------
# Targeting Alphaproteibacteria as its model prediction plot is far stronger
# than Bacillis. Breaking down into order-level as it might show which
# specific orders are responsible.

alpha_only <- from16S[which(taxonomy16S[rownames(from16S),
                                        "class"] == "Alphaproteobacteria"),
                                        -dim(from16S)[2]]

order16S_alpha <- groupBy(alpha_only, taxonomy16S, "order")

# Extract orders that have a high enough abundance
threshold <- dim(order16S_alpha)[1]
order16S_alpha <- order16S_alpha[
  ,which(colSums(order16S_alpha) > threshold)
  ]

order16S_alpha

# Removing orders with very low total counts because there are not enough
# observations to fit a reliable regression . The cutoff used is at least one
# count per sample on average

# Check for correlations with temperature and precipitation

df <- data.frame(order16S_alpha,
                 precip = about[idNums, "P.grid"],
                 temp = about[idNums, "T.grid"])


outRhizobiales = lm(Rhizobiales ~ P.grid + T.grid, data = df)
outSphingomonadales = lm(Sphingomonadales ~ P.grid + T.grid, data = df)

summary(outRhizobiales)

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) 449.60673   82.22151   5.468 3.95e-07 ***
# P.grid       -0.06370    0.01862  -3.422  0.001 ***
# T.grid       -9.11768    8.83054  -1.033  0.305
#
# Residual standard error: 151.2 on 91 degrees of freedom
# Multiple R-squared: 0.1263, Adjusted R-squared: 0.1071
# F-statistic: 6.58 on 2 and 91 DF, p-value: 0.002142

# Interpretation : Rhizobiales abundance decreases significantly with
# precipitation (p = 0.001) but temperature had no significant effect
# (p = 0.305). The model is significant overall (p = 0.002) but only
# explains about 13% of the variation, consistent with the scatter seen
# in the predicted vs actual plot

summary(outSphingomonadales)

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) 55.633821  21.527350   2.584   0.011 *
# P.grid       0.001696   0.004874   0.348   0.729
# T.grid      -4.638129   2.312025  -2.006   0.048 *
#
# Residual standard error: 39.6 on 91 degrees of freedom
# Multiple R-squared: 0.04297, Adjusted R-squared: 0.02194
# F-statistic: 2.043 on 2 and 91 DF, p-value: 0.1356

# Interpretation : Sphingomonadales abundance decreases with temperature
# (p = 0.048) but precipitation had no significant effect (p = 0.729).
# However the overall model is not significant (p = 0.136) and explains
# only 4% of the variation.

# Rhizobiales responds to precipitation while Sphingomondales responds
# weakly to temperature, suggesting different environmental drivers within the
# same order. However, both models explain very little of the variation
# so it is likely other factors are influencing this.

# Check whether model predictions are reasonable for Rhizobiales and
# Sphingomonadales

dfPlotting <- data.frame("Predicted" = predict(outRhizobiales),
                         "Actual" = df[,"Rhizobiales"])

rhizobiales <- ggplot(dfPlotting, aes(x = Predicted, y = Actual)) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  xlab('Predicted Rhizobiales abundance') +
  ylab('Actual Rhizobiales abundance') +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed",
              colour = "red",
              linewidth = 0.8)

rhizobiales

# Save it to file
ggsave("./figures/rhizobiales.tif", 
       plot = rhizobiales, 
       device = "tiff",
       width = 6, 
       height = 6,
       units = "in",
       dpi = 300)

# Interpretation : the model predictions for Rhizobiales are reasonable with a 
# general 1:1 relationship, however there is still quite a lot of scatter

dfPlotting <- data.frame("Predicted" = predict(outSphingomonadales),
                         "Actual" = df[,"Sphingomonadales"])

ggplot(dfPlotting, aes(x = Predicted, y = Actual)) +
  geom_point() +
  theme_classic() +
  xlab('Predicted Sphingomonadales abundance') +
  ylab('Actual Sphingomonadales abundance') +
  geom_abline(slope=1, intercept=0)

# Interpretation : the model predictions for Sphingomonadales are poor as almost
# all values sit near zero apart from one extreme outlier which pulls the line.


# Citations --------------------------------------------------------------------

# R Core Team (2024). R: A Language and Environment for Statistical Computing.
# R Foundation for Statistical Computing, Vienna, Austria.
# https://www.R-project.org/

# Wickham H, et al. (2019). Welcome to the tidyverse. Journal of Open Source
# Software, 4(43), 1686. https://doi.org/10.21105/joss.01686

# Oksanen J, et al. (2025). vegan: Community Ecology Package. R package
# version 2.7-2. https://CRAN.R-project.org/package=vegan

# Kassambara A, Mundt F (2026). factoextra: Extract and Visualize the Results
# of Multivariate Data Analyses. R package version 2.0.0.
# https://CRAN.R-project.org/package=factoextra