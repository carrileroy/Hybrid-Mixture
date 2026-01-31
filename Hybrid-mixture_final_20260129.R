# Script for statistical analysis for the manuscript,"Do riparian plant hybrids
# mimic leaf mixtures for in-stream litter # dynamics?" by Andrews, LeRoy, and 
# Fischer. The following analysis is reported in the results section of the 
# manuscript. 

# Compiled by Walt Andrews
# Last updated 1/28/2026

# Significant p-values (p-value < 0.05) are shown in the comments.

rm(list=ls())
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("car")) install.packages("car")
if (!require("DataExplorer")) install.packages("DataExplorer")
if (!require("stats")) install.packages("stats")
if (!require("vegan")) install.packages("vegan")
if (!require("labdsv")) install.packages("labdsv")
if (!require("ecodist")) install.packages("ecodist")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("scales")) install.packages("scales")
if (!require("goeveg")) install.packages("goeveg")
if (!require("lmPerm")) install.packages("lmPerm")

# Libraries
library(tidyverse)
library(car)
library(DataExplorer)
library(stats)
library(vegan)
library(labdsv)
library(ecodist)
library(ggpubr)
library(scales)
library(goeveg)
library(lmPerm)


# Data for analysis
lnp_afdmr <- read_csv("AFDMcjl.csv") 
p_afdmr <- read_csv("percent_AFDM-wma.csv")
litter.chem <- read_csv("Litter_chem_C_N-wma.csv")
litter.CT <- read_csv("Litter_chem_CT-wma.csv")
h1 <- read_csv("Harvest_1_inverts_from_original_taxa_list-wma.csv") 
h3 <- read_csv("harvest_3_inverts_from_original_taxa_list-wma.csv") 
h1.ept <- read_csv("Harvest_1_percent_EPT-wma.csv")
h3.ept <- read_csv("Harvest_3_percent_EPT-wma.csv")
h1.func <- read_csv("Harvest_1_invert_functional_groups-wma.csv") 
h3.func <- read_csv("Harvest_3_invert_functional_groups-wma.csv") 


# Litter chemistry
# ANOVAs for litter quality: C, N, and C:N
litter.chem
Treatment <- as_factor(litter.chem$`Litter treatment`)
levels(Treatment)

lv <- ~leveneTest(.~ Treatment, center = mean)
lv <- map(litter.chem[,2:4], lv)
sw <- ~shapiro.test(resid(lm(.~ Treatment)))
sw <- map(litter.chem[,2:4], sw)
mod <- ~Anova(lm(.~ Treatment), type = "II")
mod.results <- map(litter.chem[,2:4], mod)
c(lv, sw, mod.results)

hist((litter.chem$`% C`)^2)
leveneTest((`% C`)^2 ~ Treatment, data = litter.chem)
# fails assumptions

hist((-1/(litter.chem$`% N`)^2))
leveneTest((-1/((`% N`)^2)) ~ Treatment, data = litter.chem)
# fails assumptions

colnames(litter.chem)
hist(log(litter.chem$`C:N`))
leveneTest(log(`C:N`) ~ Treatment, data = litter.chem)
# fails assumptions

# lmperm() ANOVAs for % C, % N, and C:N
set.seed(1422)
fit <- aovp(`% C` ~ Treatment, data = litter.chem, perm = "")
summary(fit)
# Degrees of freedom = 2,3
# F = 26.52
# p-value = 0.0124*

tkh <- (TukeyHSD(aov(litter.chem$`% C` ~ Treatment), data = litter.chem))
tkh

set.seed(1473)
fit <- aovp(`% N` ~ Treatment, data = litter.chem, perm = "")
summary(fit)
# Degrees of freedom = 2,3 
# F = 35.88
# p-value = 0.00804**

tkh <- (TukeyHSD(aov(litter.chem$`% N` ~ Treatment), data = litter.chem))
tkh

set.seed(1280)
fit <- aovp(`C:N` ~ Treatment, data = litter.chem, perm = "")
summary(fit)
# Degrees of freedom = 2,3
# F = 35.95
# p-value = 0.00802**

tkh <- (TukeyHSD(aov(litter.chem$`C:N` ~ Treatment), data = litter.chem))
tkh

# Condensed tannins (CT) 
# ANOVA for litter CT using aov()
litter.CT
Treatment <- as_factor(litter.CT$`Litter treatment`)
levels(Treatment)

lv <- leveneTest(`% CT` ~ Treatment, data = litter.CT)
lv
sw <- shapiro.test(resid(lm(`% CT` ~ Treatment, data = litter.CT)))
sw
mod.results <- aov(lm(`% CT`~ Treatment, data = litter.CT), type = "II")
mod.results
# Meets assumptions
# Degrees of freedom = 2,10
# F = 64.66
# p-value < 0.00001
summary(mod.results)

set.seed(1456)
fit <- aovp(`% CT` ~ `Litter treatment`, data = litter.CT, perm = "")
summary(fit)
# p-value < 0.00001

tkh <- (TukeyHSD(aov(litter.CT$`% CT` ~ Treatment), data = litter.CT))
tkh

# Chi-squared for litter chemistry
# Subsetting by treatment
litter.chem
hybrid <- filter(litter.chem, `Litter treatment` == "Hybrid")
p.trich <- filter(litter.chem, `Litter treatment` == "P. trich.")
p.max <- filter(litter.chem, `Litter treatment` == "P. max.")

# Observed values
Obs.Hy.C <- hybrid$`% C`
Obs.Hy.N <- hybrid$`% N`
Obs.Hy.CN <- hybrid$`C:N`

# Expected values
Exp.C <- mean(c(p.trich$`% C`, p.max$`% C`))
Exp.N <- mean(c(p.trich$`% N`, p.max$`% N`))
Exp.CN <- mean(c(p.trich$`C:N`, p.max$`C:N`))

# Chi-squared analysis
X2.Hy.C <- map(Obs.Hy.C, ~ (((.x-Exp.C)^2)/Exp.C))
X2.Hy.C 
X2.Hy.C <- unlist(X2.Hy.C)
X2.Hy.C
sum(X2.Hy.C) 
# Sum = 0.5542608
# Degrees of freedom = 1
# Critical X2(4) = 3.841 
# p = 0.4566

X2.Hy.N <- map(Obs.Hy.N, ~ (((.x-Exp.N)^2)/Exp.N))
X2.Hy.N 
X2.Hy.N <- unlist(X2.Hy.N)
X2.Hy.N
sum(X2.Hy.N) 
# Sum = 0.05809609
# Degrees of freedom = 1
# Critical X2(4) = 3.841 
# p = 0.8096

X2.Hy.CN <- map(Obs.Hy.CN, ~ (((.x-Exp.CN)^2)/Exp.CN))
X2.Hy.CN 
X2.Hy.CN <- unlist(X2.Hy.CN)
X2.Hy.CN
sum(X2.Hy.CN) 
# Sum = 12.36685
# Degrees of freedom = 1
# Critical X2(4) = 3.841 
# Significant difference, p = 0.0004***

# Chi-squared test for litter CT
litter.CT
hybrid <- filter(litter.CT, `Litter treatment` == "Hybrid")
p.trich <- filter(litter.CT, `Litter treatment` == "P. trich.")
p.max <- filter(litter.CT, `Litter treatment` == "P. max.")

# Observed values for CT
Obs.Hy.CT <- hybrid$`% CT`

# Expected values for CT
Exp.CT <- mean(c(p.trich$`% CT`, p.max$`% CT`))

# Chi-squared analysis for CT
X2.Hy.CT <- map(Obs.Hy.CT, ~ (((.x-Exp.CT)^2)/Exp.CT))
X2.Hy.CT 
X2.Hy.CT <- unlist(X2.Hy.CT)
X2.Hy.CT
sum(X2.Hy.CT) 
# Sum = 8.615705
# Degrees of freedom = 3
# Critical X2(4) = 7.815 
# p = 0.034862*




# Litter decomposition analysis
# Assigning factors
lnp_afdmr
trt <- as.factor(lnp_afdmr$trt)
day <- as.factor(lnp_afdmr$day)

# ANCOVA for litter decomposition rates
leveneTest(lnpAFDMR ~ trt, lnp_afdmr)
shapiro.test(resid(lm(lnpAFDMR ~ trt*day, lnp_afdmr)))



shapiro.test(resid(lm((lnpAFDMR)^2 ~ trt*day, lnp_afdmr)))
# Fails assumptions

# Permutative ANCOVA using lmPerm
set.seed(1468)
fit <- aovp(lnpAFDMR ~ trt*day, data = lnp_afdmr)
summary(fit) 


# k rates
# fitted.values k rates
hybrid_litter <- filter(lnp_afdmr, trt == "P. max x Chill 1-13")
mix_litter <- filter(lnp_afdmr, trt == "P.max, Chill 1-13 mix")
ptrich_litter <- filter(lnp_afdmr, trt == "Chill 1-13")
pmax_litter <- filter(lnp_afdmr, trt == "P. max")

# k rates and associated values for SE included in Table 1 in manuscript
hybrid_k <- lm(lnpAFDMR ~ day, data = hybrid_litter)
hybrid_k
# k = 0.0216 
summary(hybrid_k)
# SE = 0.004534

mix_k <- lm(lnpAFDMR ~ day, data = mix_litter)
mix_k
# k = 0.02495
summary(mix_k)
# SE = 0.003451

ptrich_k <- lm(lnpAFDMR ~ day, data = ptrich_litter)
ptrich_k


# k = 0.02766
summary(ptrich_k)
# SE = 0.005147

pmax_k <- lm(lnpAFDMR ~ day, data = pmax_litter)
pmax_k
# k = 0.01194
summary(pmax_k)
# SE = 0.002197

# Percent Mass Remaining
names(p_afdmr)
p_afdmr$`Percent mass remaining`
p_afdmr$`Harvest date`
H1 <- filter(p_afdmr, `Harvest date` == "1")
H2 <- filter(p_afdmr, `Harvest date` == "2")
H3 <- filter(p_afdmr, `Harvest date` == "3")

set.seed(1479)
fit <- aovp(`Percent mass remaining` ~ `Litter`, data = H1, perm = "")
summary(fit)

set.seed(1479)
fit <- aovp(`Percent mass remaining` ~ `Litter`, data = H2, perm = "")
summary(fit)

set.seed(1479)
fit <- aovp(`Percent mass remaining` ~ `Litter`, data = H3, perm = "")
summary(fit)


# Chi-square test for percent ash free dry mass remaining;
# p-values from online calculator
# at https://www.socscistatistics.com/pvalues/chidistribution.aspx
# Subsetting by treatment
Hybrid <- filter(p_afdmr, Litter == "P. max x Chill 1-13")
Mixture <- filter(p_afdmr, Litter  == "P. max, Chill 1-13 mix")
P.trich <- filter(p_afdmr, Litter == "Chill 1-13")
P.max <- filter(p_afdmr, Litter == "P. max")

# Chi-square test for percent afdm remaining for each harvest date
# Subsetting by harvest date for each treatment
hybrid.h1 <- filter(Hybrid, `Harvest date` == "1")
mixture.h1 <- filter(Mixture,`Harvest date` == "1")
p.trich.h1 <- filter(P.trich, `Harvest date` == "1")
p.max.h1 <- filter(P.max, `Harvest date` == "1")

hybrid.h2 <- filter(Hybrid, `Harvest date` == "2")
mixture.h2 <- filter(Mixture,`Harvest date` == "2")
p.trich.h2 <- filter(P.trich, `Harvest date` == "2")
p.max.h2 <- filter(P.max, `Harvest date` == "2")

hybrid.h3 <- filter(Hybrid, `Harvest date` == "3")
mixture.h3 <- filter(Mixture,`Harvest date` == "3")
p.trich.h3 <- filter(P.trich, `Harvest date` == "3")
p.max.h3 <- filter(P.max, `Harvest date` == "3")

# For harvest 1
hybrid.h1.pafdmr <- as_tibble(hybrid.h1$`Percent mass remaining`)
hybrid.h1.pafdmr
E <- mean(c(p.max.h1$`Percent mass remaining`, p.trich.h1$`Percent mass remaining`))
E # E = 71.27461

Obs <- hybrid.h1.pafdmr
Exp <- E
Obs
Exp
X2.hybrid.h1 <- map(Obs, ~ (((.x-Exp)^2)/Exp))
X2.hybrid.h1 
X2.hy.h1 <- as_tibble(X2.hybrid.h1)
sum(X2.hy.h1) # Gives X2 value 
# Sum = 0.8818883
# Degrees of freedom (n-1) = 4
# Critical X2(4) = 9.488

mixture.h1.pafdmr <- as_tibble(mixture.h1$`Percent mass remaining`)
mixture.h1.pafdmr
E
Obs <- mixture.h1.pafdmr
Exp <- E
Obs
Exp
X2.mixture.h1 <- map(Obs, ~(((.x-Exp)^2)/Exp))
X2.mixture.h1
X2.mix.h1 <- as_tibble(X2.mixture.h1)
sum(X2.mix.h1) 
# Sum = 1.740161
# Degrees of freedom = 2
# Critical X2(2) = 5.991

# For harvest 2
hybrid.h2.pafdmr <- as_tibble(hybrid.h2$`Percent mass remaining`)
hybrid.h2.pafdmr
E <- mean(c(p.max.h2$`Percent mass remaining`, p.trich.h2$`Percent mass remaining`))
E # E = 56.33918
Obs <- hybrid.h2.pafdmr
Exp <- E
Obs
Exp
X2.hybrid.h2 <- map(Obs, ~ (((.x-Exp)^2)/Exp))
X2.hybrid.h2 
X2.hy.h2 <- as_tibble(X2.hybrid.h2)
sum(X2.hy.h2) 
# Sum = 9.407949
# Degrees of freedom = 4
# Critical X2(4) = 9.488

mixture.h2.pafdmr <- as_tibble(mixture.h2$`Percent mass remaining`)
mixture.h2.pafdmr
E
Obs <- mixture.h2.pafdmr
Exp <- E
Obs
Exp
X2.mixture.h2 <- map(Obs, ~(((.x-Exp)^2)/Exp))
X2.mixture.h2
X2.mix.h2 <- as_tibble(X2.mixture.h2)
sum(X2.mix.h2) #Gives X2 value for Mixture litter for harvest 2
# Sum = 1.364969
# Degrees of freedom = 2
# Critical X2(2) = 5.991

# For harvest 3
hybrid.h3.pafdmr <- as_tibble(hybrid.h3$`Percent mass remaining`)
hybrid.h3.pafdmr
E <- mean(c(p.max.h3$`Percent mass remaining`, p.trich.h3$`Percent mass remaining`))
E # E = 42.95913
Obs <- hybrid.h3.pafdmr
Exp <- E
Obs
Exp
X2.hybrid.h3 <- map(Obs, ~ (((.x-Exp)^2)/Exp))
X2.hybrid.h3 
X2.hy.h3 <- as_tibble(X2.hybrid.h3)
sum(X2.hy.h3) #Gives X2 value for Hybrid litter for harvest 3
# Sum = 35.32059
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# P-value < 0.0001***

mixture.h3.pafdmr <- as_tibble(mixture.h3$`Percent mass remaining`)
mixture.h3.pafdmr
E
Obs <- mixture.h3.pafdmr
Exp <- E
Obs
Exp
X2.mixture.h3 <- map(Obs, ~(((.x-Exp)^2)/Exp))
X2.mixture.h3
X2.mix.h3 <- as_tibble(X2.mixture.h3)
sum(X2.mix.h3) #Gives X2 value for Mixture litter for harvest 3
# Sum = 20.66261
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# P-value = 0.000369***




# Macroinvertebrate analysis
# Calculating diversity indices for: Abundance, S, J, H, D
# For harvest 1 and 3

# H1 diversity indices
h1
names(h1)
dim(h1)

# For Harvest 1 Macroinvertebrates
# Abundance
abundance.data <- data.mod <- h1[,5:23] %>%
  mutate(Abundance = rowSums(.))
abundance.data <- print(data.mod)
Abundance <- abundance.data$Abundance
Abundance

# Richness (S)
# This code for richness should be adjusted for each different data frame
richness.data <- h1[,5:23] %>%
  mutate(S = rowSums(.[2:(ncol(.))] > 0))
S <- richness.data$S
S

# Shannon-Wiener index (H)
H <- diversity(h1[,5:23])  # for the Shannon index
H

# Pielou's evenness (J')
J.prime <- H/log(S)
J.prime

# Simpson's index (D)
D <- diversity(h1[,5:23], "simpson")  # Simpson index
D

# Combining indices into a tibble with the original data frame
names(h1)
diversity.data <-tibble(Abundance, S, J.prime, H, D)
h1.diversity.indexes <- bind_cols(h1[,3], diversity.data)
head(h1.diversity.indexes)
dim(h1.diversity.indexes)


# For harvest 3 macroinvertebrates
# Abundance
abundance.data <- data.mod <- h3[,5:34] %>%
  mutate(Abundance = rowSums(.))
abundance.data <- print(data.mod)
Abundance <- abundance.data$Abundance
Abundance

# Richness (S)
# This code for richness should be adjusted for each different data frame
richness.data <- h3[,5:34] %>%
  mutate(S = rowSums(.[2:(ncol(.))] > 0))

S <- richness.data$S
S

# Shannon-Weiner index (H)
H <- diversity(h3[,5:34])  # for the Shannon index
H

# Pielou's evenness (J')
J.prime <- H/log(S)
J.prime

# Simpson's index (D)
D <- diversity(h3[,5:34], "simpson")  # Simpson index
D


# Combine various indices into a tibble with the original data frame
names(h3)
diversity.data <-tibble(Abundance, S, J.prime, H, D)
h3.diversity.indexes <- bind_cols(h3[,3], diversity.data)
head(h3.diversity.indexes)
dim(h3.diversity.indexes)

# file_path <- "/Users/waltonandrews/Downloads/h3.diversity.indexes.csv"
# write_csv(h3.diversity.indexes, file = file_path)


# One-way ANOVAs for macroinvertebrate data
# Harvest 1 ANOVA for macroinvertebrate diversity indexes
h1.diversity.indexes
h1.div <- h1.diversity.indexes
Treatment <- as_factor(h1.div$`Litter treatment`)
lv <- (~leveneTest(.~Treatment))
map(h1.div[,2:6], lv)
sw <- (~shapiro.test(resid(lm(.~ Treatment))))
map(h1.div[,2:6], sw)

# Analysis for S
hist(h1.div$S)
hist(sqrt(h1.div$S))
shapiro.test(resid(lm(sqrt(h1.div$S) ~ Treatment)))
# Passes assumptions
Anova(lm(h1.div$S ~ Treatment), type = "III")

# Harvest 3 ANOVA for macroinvertebreate diversity indexes
h3.diversity.indexes
h3.div <- h3.diversity.indexes
Treatment <- as_factor(h3.div$`Litter treatment`)
lv <- (~leveneTest(.~Treatment))
map(h3.div[,2:6], lv)
sw <- (~shapiro.test(resid(lm(.~ Treatment))))
map(h3.div[,2:6], sw)
# D failed S-W test
h3.div[,2:6] %>%
  map(~Anova(lm(. ~ Treatment), type = "III")) 
h3.div[,2:6] %>%
  map(~summary(Anova(lm(. ~ Treatment), type = "III"))) 
# These results are included in Table 2.
# Abundance doesn't show a significant difference, p = 0.082326
# F(3,17) = 2.7457

# S doesn't show a significant difference, p =  0.0777
# F(3,17) = 2.8132

# J' doesn't show a significant difference, p = 0.2897 
# F(3,17) = 1.3909

# Significant difference for H, p-value = 0.008037**
# F(3, 17) =  5.9044

# Transformation and analysis for D 
hist(h3.div$D)
hist(h3.div$D^2)
lv <- leveneTest(((h3.div$D)^2) ~ Treatment)
lv
sw <- shapiro.test(resid(lm(((h3.div$D)^2)~ Treatment)))
sw
# Passes assumptions
Anova(lm(h3.div$D^2 ~ Treatment), type = "III") 
# For D, p-value = 0.003849**
# F(3,17) = 7.1303

# Continued, R-squared results included in Table 2 in manuscript
# R-squared values at day 41
h3.div
summary(lm(h3.div$Abundance ~ Treatment)) 
# Abundance R-squared: 0.3704
summary(lm(h3.div$S ~ Treatment)) 
# S R-squared: 0.3761
summary(lm(h3.div$J.prime ~ Treatment)) 
# R-squared: 0.243
summary(lm(h3.div$H ~ Treatment)) 
# R-squared: 0.5585
summary(lm(h3.div$D^2 ~ Treatment))
# Adjusted R-squared: 0.6044

# Permutation tests
set.seed(1456)
fit <- aovp(`Abundance` ~ `Litter treatment`, data = h3.div, perm = "")
summary(fit)

set.seed(1479)
fit <- aovp(`S` ~ `Litter treatment`, data = h3.div, perm = "")
summary(fit)

set.seed(1403)
fit <- aovp(`J.prime` ~ `Litter treatment`, data = h3.div, perm = "")
summary(fit)

set.seed(1448)
fit <- aovp(`H` ~ `Litter treatment`, data = h3.div, perm = "")
summary(fit)
# p-value = 0.00804**

set.seed(1427)
fit <- aovp(`D` ~ `Litter treatment`, data = h3.div, perm = "")
summary(fit)
# p-value = 0.0331*


h3.div
hybrid <- filter(h3.div, `Litter treatment` == "Hybrid")
mixture <- filter(h3.div, `Litter treatment` == "Mixture")
p.trich <- filter(h3.div, `Litter treatment` == "P. trich.")
p.max <- filter(h3.div, `Litter treatment` == "P. max.")

mean(hybrid$H)
# Mean hybrid H = 1.67807
mean(mixture$H)
# Mean mixture H = 1.645688
mean(p.trich$H)
# Mean P. trich. H = 1.114851
mean(p.max$H)
# Mean P. max. H = 0.08551023

mean(hybrid$D)
# Mean hybrid D = 0.7581925 
mean(mixture$D)
# Mean mixture D = 0.736923
mean(p.trich$D)
# Mean P. trich. D = 0.5658635
mean(p.max$D)
# Mean P. max. D = 0.4533086


# Post-hoc pairwise comparisons
tkh <- (TukeyHSD(aov(H~Treatment), data = h3.div))
tkh
# Significant difference for H, P. max.-Hybrid, p = 0.0120191*
# Significant difference for H, P. max.-Mixture, p = 0.0383774*
tkd <- (TukeyHSD(aov(D~Treatment), data = h3.div))
tkd
# Significant difference for D, P. trich.-Hybrid, p = 0.0266548* 
# Significant difference for D, P. max.- Hybrid, p =0.0071556**
# Significant difference for D, P.max.-Mixture, p = 0.0417552*




# Chi-squared calculations for observed hybrid and mixture abundance and 
# macroinvertebrate diversity values compared to expected (parental average),
# including DF and p-values from online calculator
# at https://www.socscistatistics.com/pvalues/chidistribution.aspx
# Analysis for harvest 3
# Response variables: Abundance, Richness (S), Evenness (Pielou's J'), 
# Shannon-Wiener Index (H), and Simpson's Index (D)

# Subsetting by treatment
h3
hybrid <- filter(h3.div, `Litter treatment` == "Hybrid")
mixture <- filter(h3.div, `Litter treatment` == "Mixture")
p.trich <- filter(h3.div, `Litter treatment` == "P. trich.")
p.max <- filter(h3.div, `Litter treatment` == "P. max.")
p.max <- na.omit(p.max)

# Observed values
Obs.Hy.A <- hybrid$Abundance
Obs.Hy.S <- hybrid$S
Obs.Hy.J.prime <- hybrid$J.prime
Obs.Hy.H <- hybrid$H
Obs.Hy.D <- hybrid$D

Obs.Mix.A <- mixture$Abundance
Obs.Mix.S <- mixture$S
Obs.Mix.J.prime <- mixture$J.prime
Obs.Mix.H <- mixture$H
Obs.Mix.D <- mixture$D

# Expected values
Exp.A <- mean(c(p.trich$Abundance, p.max$Abundance))
Exp.S <- mean(c(p.trich$S, p.max$S))
Exp.J.prime <- mean(c(p.trich$J.prime, p.max$J.prime))
Exp.H <- mean(c(p.trich$H, p.max$H))
Exp.D <- mean(c(p.trich$D, p.max$D))

# Chi-squared analysis
X2.Hy.A <- map(Obs.Hy.A, ~ (((.x-Exp.A)^2)/Exp.A))
X2.Hy.A 
X2.Hy.A <- unlist(X2.Hy.A)
X2.Hy.A
sum(X2.Hy.A) 
# Sum = 18.42765
# Degrees of freedom = 4
# Critical X2(4) = 9.488 
# for hybrid abundance, p = 0.0001**

# For hybrid treatment observed compared to expected
X2.Hy.S <- map(Obs.Hy.S, ~ (((.x-Exp.S)^2)/Exp.S))
X2.Hy.S 
X2.Hy.S <- unlist(X2.Hy.S)
X2.Hy.S
sum(X2.Hy.S) 
# Sum = 10.8
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# For hybrid richness, p = 0.028906*

X2.Hy.J.prime <- map(Obs.Hy.J.prime, ~ (((.x-Exp.J.prime)^2)/Exp.J.prime))
X2.Hy.J.prime 
X2.Hy.J.prime <- unlist(X2.Hy.J.prime)
X2.Hy.J.prime
sum(X2.Hy.J.prime) 
# Sum = 0.1271504
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# For hybrid J', p = 0.998064

X2.Hy.H <- map(Obs.Hy.H, ~ (((.x-Exp.H)^2)/Exp.H))
X2.Hy.H 
X2.Hy.H <- unlist(X2.Hy.H)
X2.Hy.H
sum(X2.Hy.H) 
# Sum = 1.922914
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# For hybrid H, p-value = 0.749937

X2.Hy.D <- map(Obs.Hy.D, ~ (((.x-Exp.D)^2)/Exp.D))
X2.Hy.D 
X2.Hy.D <- unlist(X2.Hy.D)
X2.Hy.D
sum(X2.Hy.D) 
# Sum = 0.3793533
# Degrees of freedom = 4
# Critical X2(4) = 9.488
# For hybrid D, p = 0.984136

# For mixture treatment observed compared to expected
X2.Mix.A <- map(Obs.Mix.A, ~ (((.x-Exp.A)^2)/Exp.A))
X2.Mix.A 
X2.Mix.A <- unlist(X2.Mix.A)
X2.Mix.A
sum(X2.Mix.A) 
# Sum = 28.7345
# Degrees of freedom = 2
# Critical X2(2) = 5.991
# For abundance, p-value < 0.00001***

X2.Mix.S <- map(Obs.Mix.S, ~ (((.x-Exp.S)^2)/Exp.S))
X2.Mix.S 
X2.Mix.S <- unlist(X2.Mix.S)
X2.Mix.S
sum(X2.Mix.S) 
# Sum = 6.6 
# Degrees of freedom = 2 
# Critical X2(2) = 5.991
# For mixture richness, p-value = 0.036883*

X2.Mix.J.prime <- map(Obs.Mix.J.prime, ~ (((.x-Exp.J.prime)^2)/Exp.J.prime))
X2.Mix.J.prime 
X2.Mix.J.prime <- unlist(X2.Mix.J.prime)
X2.Mix.J.prime
sum(X2.Mix.J.prime) 
# Sum = 0.02606923
# Degrees of freedom = 2
# Critical X2(2) = 5.991
# For mixture J', p = 0.987084

X2.Mix.H <- map(Obs.Mix.H, ~ (((.x-Exp.H)^2)/Exp.H))
X2.Mix.H 
X2.Mix.H <- unlist(X2.Mix.H)
X2.Mix.H
sum(X2.Mix.H) 
# Sum = 0.8372721
# Degrees of freedom = 2
# Critical X2(2) = 5.991
# For mixture H, p = 0.657967

X2.Mix.D <- map(Obs.Mix.D, ~ (((.x-Exp.D)^2)/Exp.D))
X2.Mix.D 
X2.Mix.D <- unlist(X2.Mix.D)
X2.Mix.D
sum(X2.Mix.D) 
# Sum = 0.1568114
# Degrees of freedom = 2
# Critical X2(2) = 5.991
# For mixture D,  p = 0.924595

ggplot(h3, aes(`Litter treatment`, Abundance)) + geom_boxplot()
ggplot(h3, aes(`Litter treatment`, S)) + geom_boxplot()
ggplot(h3, aes(`Litter treatment`, J.prime)) + geom_boxplot()
ggplot(h3, aes(`Litter treatment`, H)) + geom_boxplot()
ggplot(h3, aes(`Litter treatment`, D)) + geom_boxplot()


# Community analysis
# Ordinations

# NMDS for Harvest 1
h1 
dim(h1)
comm<-h1[,5:23]
trt<-h1[,3]
comm
dis <- vegdist(comm, method = "bray")
scree <- nmds(dis, mindim = 1, maxdim = 5, nits = 1)
stress <- scree$stress
plot(stress, xlab = "Number of Axes")
set.seed(123) 
m.mds <- metaMDS(comm, k = 3, trymax = 1001)
stressplot(m.mds)

m_comm<-as.matrix(comm)
m_comm
mrpp(dat = comm, trt$`Litter treatment`, distance = "bray")


# NMDS for harvest 3
h3
dim(h3)
comm<-h3[,5:34]
trt<-h3[,3]
comm
dis <- vegdist(comm, method = "bray")
scree <- nmds(dis, mindim = 1, maxdim = 5, nits = 1)
stress <- scree$stress
plot(stress, xlab = "Number of Axes")
set.seed(123) 
m.mds <- metaMDS(comm, k = 2, trymax = 1001)
stressplot(m.mds)

m_comm<-as.matrix(comm)
m_comm
mrpp(dat = comm, trt$`Litter treatment`, distance = "bray")
#Non-significant, A = 0.01586, significance of delta = 0.246





# Functional group analysis
# For percent shredders
# Harvest 1
h1.func
Treatment <- as_factor(h1.func$`Litter treatment`)
lv <- leveneTest(sqrt(h1.func$`Percent shredders`) ~Treatment)
lv
sw <- shapiro.test(resid(lm(sqrt(h1.func$`Percent shredders`) ~ Treatment)))
sw
Anova(lm(sqrt(`Percent shredders`) ~ Treatment, data = h1.func), Type = "II")

set.seed(1431)
fit <- aovp(`Percent shredders` ~ `Litter treatment`, data = h1.func, perm = "")
summary(fit)

# Harvest 3
h3.func
Treatment <- as_factor(h3.func$`Litter treatment`)
lv <- leveneTest(h3.func$`Percent shredders` ~Treatment)
lv
sw <- shapiro.test(resid(lm(h3.func$`Percent shredders` ~ Treatment)))
sw
Anova(lm(`Percent shredders` ~ Treatment, data = h3.func), Type = "II")

Treatment <- as_factor(h3.func$`Litter treatment`)
lv <- leveneTest(h3.func$`Percent predators` ~Treatment)
lv
sw <- shapiro.test(resid(lm(h3.func$`Percent predators` ~ Treatment)))
sw
Anova(lm(`Percent predators` ~ Treatment, data = h3.func), Type = "II")

set.seed(1441)
fit <- aovp(`Percent shredders` ~ `Litter treatment`, data = h3.func, perm = "")
summary(fit)

# Harvest 1 for all functional groups
h1.func
Treatment <- as_factor(h1.func$`Litter treatment`)
lv <- (~leveneTest(.~Treatment))
map(h1.func[,5:12], lv)
sw <- (~shapiro.test(resid(lm(.~ Treatment))))
map(h1.func[,5:12], sw)
h1.func[,5:12] %>%
  map(~Anova(lm(. ~ Treatment))) 

# Transforming data for data that failed assumptions
# All other variables were not significant
# Shredders
lv <- leveneTest(sqrt(h1.func$`Shredders`) ~Treatment)
lv
sw <- shapiro.test(resid(lm(sqrt(h1.func$`Shredders`) ~ Treatment)))
sw
Anova(lm(sqrt(`Shredders`) ~ Treatment, data = h1.func), Type = "II")

# Scrapers
lv <- leveneTest(sqrt(h1.func$`Scrapers`) ~Treatment)
lv
sw <- shapiro.test(resid((lm(sqrt(h1.func$`Scrapers`) ~ Treatment))))
sw
hist(-1/(h1.func$`Scrapers`)^2)

# Permutative ANOVA for Scrapers
set.seed(1437)
fit <- aovp(`Scrapers` ~ Treatment, data = h1.func)
summary(fit)

# Collector-filterers
lv <- leveneTest(h1.func$`Collector-filterers` ~Treatment)
lv
sw <- shapiro.test(resid(lm(h1.func$`Collector-filterers` ~ Treatment)))
sw
hist(h1.func$`Collector-filterers`)
sw <- shapiro.test(resid(lm(-1/(h1.func$`Collector-filterers`)^2 ~ Treatment)))
sw

set.seed(1421)
fit <- lmp(h1.func$`Collector-filterers` ~ Treatment, data = h1.func)
summary(aov(fit))

# Predator-shredders
lv <- leveneTest(h1.func$`Predator-shredders` ~Treatment)
lv
sw <- shapiro.test(resid(lm(h1.func$`Predator-shredders` ~ Treatment)))
sw
hist(h1.func$`Predator-shredders`)
sw <- shapiro.test(resid(lm(-1/(h1.func$`Predator-shredders`)^2 ~ Treatment)))
sw

set.seed(1474)
fit <- lmp(h1.func$`Predator-shredders` ~ Treatment, data = h1.func)
summary(aov(fit))


# Harvest 3 for all functional groups
h3.func
Treatment <- as_factor(h3.func$`Litter treatment`)

lv <- (~leveneTest(.~Treatment))
map(h3.func[,5:13], lv)
sw <- (~shapiro.test(resid(lm(.~ Treatment))))
map(h3.func[,5:13], sw)
h3.func[,5:13] %>%
  map(~Anova(lm(. ~ Treatment))) 

# Shredders
lv <- leveneTest(sqrt(h3.func$`Shredders`) ~Treatment)
lv
sw <- shapiro.test(resid(lm(sqrt(h3.func$`Shredders`) ~ Treatment)))
sw
Anova(lm(sqrt(`Shredders`) ~ Treatment, data = h3.func), Type = "II")

# Predators
lv <- leveneTest((h3.func$`Predators`) ~Treatment)
lv
sw <- shapiro.test(resid(lm(-1/(h3.func$`Predators`)^2 ~ Treatment)))
sw

# Permutative ANOVA for Predators
set.seed(1459)
fit <- lmp(h3.func$`Predators` ~ Treatment, data = h3.func, perm = "")
summary(aov(fit))
# Degrees of freedom: 3,14
# F = 4.605
# p = 0.0192*
ggplot(h3.func, aes(Treatment, Predators)) + geom_boxplot()

set.seed(1421)
fit <- aovp(`Predators` ~ Treatment, data = h3.func, perm = "")
summary(fit)

# Permutative ANOVA for Percent Predators
set.seed(1459)
fit <- lmp(h3.func$`Percent predators` ~ Treatment, data = h3.func, perm = "")
summary(aov(fit))
# Degrees of freedom: 3,14
# F = 5.696
# p = 0.000918**
ggplot(h3.func, aes(Treatment, `Percent predators`)) + geom_boxplot()

# Predators-shredders
lv <- leveneTest(sqrt(h3.func$`Predators-shredders`) ~Treatment)
lv
sw <- shapiro.test(resid(lm(sqrt(h3.func$`Predators-shredders`) ~ Treatment)))
sw
Anova(lm(sqrt(`Predators-shredders`) ~ Treatment, data = h3.func), Type = "II")

# Collector-gatherers
lv <- leveneTest(h3.func$`Collector-gatherers` ~Treatment)
lv
hist(h3.func$`Collector-gatherers`)
sw <- shapiro.test(resid(lm(-1/(h3.func$`Collector-gatherers`)^2 ~ Treatment)))
sw

# Permutative ANOVA for collector-gatherers
set.seed(1497)
fit <- lmp(h3.func$`Collector-gatherers` ~ Treatment, data = h3.func)
summary(aov(fit))

# Shredder-collectors
lv <- leveneTest(h3.func$`Shredder-collectors` ~Treatment)
lv
hist(h3.func$`Shredder-collectors`)
sw <- shapiro.test(resid(lm(-1/(h3.func$`Collector-filterers`)^2 ~ Treatment)))
sw

# Permutative ANOVA for shredder-collectors
set.seed(1473)
fit <- lmp(h3.func$`Shredder-collectors` ~ Treatment, data = h3.func)
summary(aov(fit))

# Collector-filterers
lv <- leveneTest((h3.func$`Collector-filterers`)^2 ~Treatment)
lv
sw <- shapiro.test(resid(lm((h3.func$`Collector-filterers`)^2 ~ Treatment)))
sw

# Permutative ANOVA for collector-filterers  
set.seed(1488)
fit <- lmp(h3.func$`Collector-filterers` ~ Treatment, data = h3.func)
summary(aov(fit))




# Percent EPT 
# Harvest 1 EPT
h1.ept
Litter <- as_factor(h1.ept$`Litter treatment`)

lv <- leveneTest(h1.ept$`% EPT` ~ Litter)
lv
sw <- shapiro.test(resid(lm(h1.ept$`% EPT` ~ Litter)))
sw
ggplot(h1.ept, aes(x = `% EPT`)) + geom_histogram()
ggplot(h1.ept, aes(x = sqrt(`% EPT`))) + geom_histogram()

# Using sqrt transformation #Double check that this is the correct transformation to use
lv <- leveneTest(sqrt(h1.ept$`% EPT`) ~ Litter)
lv
sw <- shapiro.test(resid(lm(sqrt(`% EPT`) ~ Litter, data = h1.ept)))
sw
Anova(lm(sqrt(`% EPT`) ~ Litter, data = h3.ept), Type = "II")

# Harvest 3 EPT
h3.ept
Litter <- as_factor(h3.ept$`Litter treatment`)

lv <- leveneTest(`% EPT` ~ Litter, data = h3.ept)
lv
sw <- shapiro.test(resid(lm(`% EPT` ~ Litter, data = h3.ept)))
sw
Anova(lm(`% EPT` ~ Litter, data = h3.ept), Type = "II")

set.seed(1478)
fit <- aovp(`% EPT` ~ Litter, data = h3.ept, perm = "")
summary(fit)


# Indicator species analysis using the indval() function in the labdsv package 
# as specified in Dufrene and Legendre 1997
# Harvest 1 indicator species analysis
h1
dim(h1)
spec <- h1[,5:23]
head(spec)
Treatment <- h1$`Litter treatment`
h1.ind <- indval(spec, Treatment, numitr=1001)
summary(h1.ind)

# Harvest 3 indicator species analysis
h3
head(h3)
dim(h3)
spec1 <- h3[,5:34]
head(spec1)
Treatment <- h3$`Litter treatment`
h3.ind <- indval(spec1, Treatment, numitr=1001)
summary(h3.ind)
# For Tipulidae A, p = 0.02497502*







