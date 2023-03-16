####### community cage experiment #######
####### Jinlin Chen #######
####### Last update: 2023 Feb #######

### Table of Content
### 1. dependencies
### 2. data formating
### 3. establishment of PANDORA
### 4. community structure
### 5. impact on individual species (census)
### 6. temperal trends of individual species (sampling)
### 7. infering indirect effects
### 8. climate record

############## dependencies ##############
# libraries
library(tidyverse)
library(cowplot) # plot multiple ggplot into one graph.
library(emmeans) # multiple conparisons
library(lme4) # mix effect model
library(brms) # for fitting bayesian linear model 
library(vegan) # conducting NMDS analysis 
library(multcompView) # to create the levels from multiple comparison

# color code
# color-blind friendly: https://davidmathlogic.com/colorblind/#%23E1BE6A-%2340B0A6
temp_cc <- c("upland" = "#003333", "heatwave" = "#990000", "warming" = "#CC6600")
#temp_cc <- c("upland" = "#003333", "heatwave" = "#A0438A", "warming" = "#949436")
treatment_cc <- c("HW_CLOSED" = "#FF0033", "HW_INVAS" = "#990000",
                  "CTRL_CLOSED" = "#009999", "CTRL_INVAS" = "#003333",
                  "WARM_CLOSED" = "#FFCC99", "WARM_INVAS" = "CC6600")
spSex_cc <- c("PAL.M" = "dark blue", "PAL.F" = "blue", "RUB.M" = "dark green", "RUB.F" = "green",
              "PST.M" = "purple", "PAN.M" = "red", "BIR.M" = "yellow", "unID.F" = "orange",
              "Asobara" = "black", "Diapriid" = "gray")
sp_shape <- c("PAL" = "dotted", "PST" = "dashed", "PAN" = "solid")
sp_cc <- c("PAL" = "#076F53", "PST" = "#565252", "PAN" = "#772553")

############## data formating ##############
# census data
dat_census <- read.csv("commCageCensus.csv")
dat_census <- dat_census %>% 
  mutate(treatment = paste0(temp, "_", invasion))
# mutate(PAL.prop = (PAL.M + PAL.F)/(PAL.M+PAL.F+RUB.M+RUB.F+BIR.M+PST.M+PAN.M+unID.F), 
#        RUB.prop = (RUB.M + RUB.F)/(PAL.M+PAL.F+RUB.M+RUB.F+BIR.M+PST.M+PAN.M+unID.F),
#        BIR.prop = BIR.M/(PAL.M+RUB.M+BIR.M+PST.M+PAN.M),
#        PST.prop = PST.M/(PAL.M+RUB.M+BIR.M+PST.M+PAN.M),
#        PAN.prop = PAN.M/(PAL.M+RUB.M+BIR.M+PST.M+PAN.M)
#        )
swfun.cc <- function(x){switch (x,
                                "HW_CLOSED" = "#FF0033", "HW_INVAS" = "#990000",
                                "CTRL_CLOSED" = "#009999", "CTRL_INVAS" = "#003333",
                                "WARM_CLOSED" = "#FFCC99", "WARM_INVAS" = "CC6600")}
dat_census$cc <- sapply(as.character(dat_census$treatment), swfun.cc)
dat_census <- dat_census %>%
  mutate(temp = ifelse(temp == "CTRL", "upland", 
                       ifelse(temp == "HW", "heatwave", "warming"))) %>% 
  mutate(invasion = ifelse(invasion == "INVAS", "introduction", "without"))
dat_census$invasion <- factor(dat_census$invasion, levels = c("without", "introduction"))
dat_census$temp <- factor(dat_census$temp, levels = c("upland", "warming", "heatwave"))
# wide data format is useful in analyzing community composition. long format is useful to plotting 
# transform to long data format
dat_census_long <- dat_census %>% gather(spSex, value, -c(vID, sID, week, cageID, block, temp, invasion, rep, treatment, cc))
dat_census_long$spSex <- factor(dat_census_long$spSex, levels = c("PAL.F", "PAL.M", "RUB.M", "RUB.F", "PAN.M", "PST.M", "BIR.M",  "unID.F", "Asobara", "Diapriid"))

# regular samples
dat_sample <- read.csv("commCageSample.csv")
dat_sample <- dat_sample %>%
  mutate(treatment = paste0(temp, "_", invasion)) %>%
  mutate(temp = ifelse(temp == "CTRL", "upland", 
                       ifelse(temp == "HW", "heatwave", "warming"))) %>% 
  mutate(invasion = ifelse(invasion == "INVAS", "introduction", "without"))
dat_sample$invasion <- factor(dat_sample$invasion, levels = c("without", "introduction"))
dat_sample$temp <- factor(dat_sample$temp, levels = c("upland", "warming", "heatwave"))

dat_sample_perCage <- dat_sample %>%
  dplyr::group_by(sID, temp, invasion, block, cageID, treatment) %>%
  dplyr::summarize(PAL.M = sum(PAL.M, na.rm = TRUE), RUB.M = sum(RUB.M, na.rm = TRUE),  # both plyr and dplyr have functions called summarize, need to use the dplyr
            BIR.M = sum(BIR.M, na.rm = TRUE), PST.M = sum(PST.M, na.rm = TRUE), # use the sum of three vials from the same cage --- that is the new input each week into the community
            PAN.M = sum(PAN.M, na.rm = TRUE),
            Asobara = sum(Asobara, na.rm = TRUE), Diapriid = sum(Diapriid, na.rm = TRUE))
dat_sample_perCage <- dat_sample_perCage %>% # no data of PAN is collcted in the first 3 samplings
  mutate(totalFly = ifelse(sID=="S0.3"|sID=="S0.1"|sID=="S0.2", PAL.M+RUB.M+BIR.M+PST.M, PAL.M+RUB.M+BIR.M+PST.M+PAN.M)) %>% 
  mutate(PAL.prop = PAL.M/totalFly,
         RUB.prop = RUB.M/totalFly,
         BIR.prop = BIR.M/totalFly,
         PST.prop = PST.M/totalFly,
         PAN.prop = PAN.M/totalFly)

# heatwave sample
hwSample <- read.csv("hwSample.csv")
hwSample <- hwSample %>% 
  mutate(PAL = PAL.M + PAL.F, PAN = PAN.F + PAN.M, PST = PST.F + PST.M) %>% 
  dplyr::select(cageID, block, timing, temp, invasion, rep, PAL, PAN, PST) %>% 
  gather(species, value, -c(cageID, block, timing, temp, invasion, rep))
hwSample$species <- factor(hwSample$species, levels = c("PAL", "PST", "PAN"))
hwSample <- hwSample %>% mutate(timing = ifelse(timing == "preH", "before HW", 
                                    ifelse(timing == "heat5d", "end of HW", "5d after HW")))
hwSample$timing <- factor(hwSample$timing, levels = c("before HW", "end of HW", "5d after HW"))
hwSample_perCage <-  hwSample %>%
  group_by(cageID, block, timing, temp, invasion, species) %>% 
  summarize(totalOffspring = sum(value)) %>% # add up the offspring of the three vials from the same cage
  ungroup()

############## establishment of PANDORA  ##############
## plots
p1PAN <- dat_census %>%
  filter(invasion == "introduction") %>%
  ggplot(aes(x = temp, y = (PAN.M))) + 
  geom_point(aes(color = temp, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = temp, color = temp), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = temp, color = temp), size = 1, width = 0.1, alpha = 0.6) + 
  #geom_boxplot(aes(group = treatment)) + 
  ylab("Number of male D. pandora") + xlab("Temperature treatment") + 
  scale_color_manual(values = temp_cc, name = "temperature") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = c(0.2, 0.85)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))

p2PAN <- dat_sample_perCage %>%
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  filter(invasion == "introduction") %>%
  ggplot(aes(x = sID, y = (PAN.M), group = temp)) + 
  geom_point(aes(color = temp, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = temp, color = temp), size = 2) +
  stat_summary(fun=mean, geom="line", aes(group = temp, color = temp), size = 1, alpha = 0.6) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = temp, color = temp), size = 1, width = 0.1, alpha = 0.6) + 
  xlab("Week") + ylab("Number of male offspring") + 
  scale_color_manual(values = temp_cc, name = "temprature") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))

p3PAN <- hwSample_perCage %>%
  filter(invasion != "CLOSED") %>%
  ggplot(aes(x = timing, y = (totalOffspring))) +
  geom_point(aes(color = species, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = species, color = species), size = 2) +
  stat_summary(fun=mean, geom="line", aes(group = species, linetype = species, color = species), size = 1) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = species, color = species), size = 1, width = 0.1, alpha = 0.6) + 
  scale_color_manual(values = sp_cc, name = "species") + 
  scale_linetype_manual(values = sp_shape, name = "species") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  ylab("Number of offspring") + xlab("Sample time relative to heatwave (HW)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = c(0.45, 0.85)) +
  theme(text = element_text(family = "Times"))

# Figure 2: heatwave facillitate establishment of D.pandora at upland
plot_grid(p1PAN, p2PAN, p3PAN, nrow = 1, ncol = 3)
ggsave("figure3_pandora_establishment.png", width = 12, height = 7)

# stats
# census
# info on how to inteprete glm result: https://www.statology.org/interpret-glm-output-in-r/
mod.pan.census <- dat_census %>%
  filter(invasion == "introduction") %>%
  glm(formula = PAN.M ~ temp * block, family = poisson)
plot(mod.pan.census)
summary(mod.pan.census)
# get the estimated mean for each treatment
emmeans(mod.pan.census, specs = pairwise ~ temp, type = "response")$emmeans

# sampling data (reproduction)
# the rule to specify nested random effect: https://www.muscardinus.be/2017/07/lme4-random-effects/
# because in my data, cageID is a unique lable that doesn't overlap between block, (1|block/cageID) is equivalent to (1|block) + (1|cageID)
mod.pan.sample <- dat_sample_perCage %>%
  filter(invasion == "introduction") %>%
  glmer(formula = PAN.M ~ temp + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pan.sample)
summary(mod.pan.sample)
emmeans(mod.pan.sample, specs = pairwise ~ temp, type = "response")$emmeans

# heatwave sampling - effect of heatwave on fly reproduction
mod.hwSample0 <- hwSample_perCage %>%
  filter(species != "PAN" | invasion != "CLOSED") %>% ## in "CLOSED" treatment, number of PAN is def zero, this is not real data point
  glm(formula = totalOffspring ~ timing*species + invasion + block, family = poisson)
plot(mod.hwSample0)
mod.hwSample <- hwSample_perCage %>%
  filter(species != "PAN" | invasion != "CLOSED") %>% ## in "CLOSED" treatment, number of PAN is def zero, this is not real data point
  glm(formula = log(totalOffspring+1) ~ timing*species + invasion + block, family = gaussian)
plot(mod.hwSample)
summary(mod.hwSample)
emm.hwSample <- emmeans(mod.hwSample, specs = pairwise ~ timing|species, type = "response")
emm.hwSample$emmeans
emm.hwSample$contrasts %>% summary(infer = TRUE)

sink("PAN_results.txt")
print("========= PAN CENSUS ==========")
summary(mod.pan.census)
emmeans(mod.pan.census, specs = pairwise ~ temp, type = "response")$emmeans
print("========= PAN SAMPLING ==========")
summary(mod.pan.sample)
summary(mod.pan.sample_indirect)
print("========= PAN HW ==========")
emm.hwSample$emmeans
emm.hwSample$contrasts %>% summary(infer = TRUE)
sink()


############## community structure ###############
## plot host and wasps (different trophic levels) separately
## all upland fly species 
# extract relevant data
composition <- dat_census %>%
  select(cageID, treatment, cc, PAL.M, RUB.M, BIR.M, PST.M)

# calculate NMDS
resident_NMDS=metaMDS(composition[, 4:7],k=2,trymax=100)
stressplot(resident_NMDS) # Run 20 stress 0.1110647 
residentNMDS.spp.fit <- envfit(resident_NMDS, composition[, 4:7], permutations = 999)
head(residentNMDS.spp.fit) # again, their abs values are much larger than those of each cage - not good looking if plotting them in reaal scale

# plot in ggplot
# extract NMDS axis values
site.scrs <- as.data.frame(scores(resident_NMDS, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, treatment = composition$treatment) #add grouping variable "Management" to dataframe
# finding the minimal convex hull that contains all community within each treatment
find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]
hulls <- plyr::ddply(site.scrs, "treatment", find_hull)
# Figure 3a: NMDs of upland flies 
NMDS_plot_flies <- site.scrs %>% ggplot(aes(x=NMDS1, y=NMDS2, color = treatment, fill = treatment))+ #sets up the plot
  geom_point(size = 2, shape = 7)+ #adds site points to plot
  geom_polygon(data = hulls, alpha = 0.5) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  scale_fill_manual(values = treatment_cc, name = "Treatment") + 
  ggtitle("resident flies") + 
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "right")
# ggsave("NMDS_uplandFlies.png", width = 8, height = 6)

# statistic difference
# The ANOSIM statistic compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups.Significance of the R statistic is determined by permuting group membership a large number of times to obtain the null distribution of the R statistic. 
# ANOSIM statistic R: An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups”
# ANOSIM significance: similar to p value
# note: Running ANOSIM on groups with very different dispersions (size of the polygon) can lead to unreliable results. Groups with very different dispersions may produce high R values, even if there's no real difference in their centroids. If differences in group dispersion are as meaningful to your analysis as differences in group centre, this may not be an issue.
anosim(composition[, 4:7], composition$treatment, distance = "bray", permutations = 9999) # Bray-Curtis dissimilarity is an asymmetrical measure often used for raw count data.
# NPMANOVA (non-parametric MANOVA, as community composition data is usually not normal distributed)
# NPMANOVA uses permutation to assess the significance of the pseudo F-statistic described above.
# It is generally accepted that any separation between groups is not significant if more than ~ 5% of the permuted F-statistics have values greater than that of the observed statistic (i.e. a P-value > 0.05).  
resident_adonis <- adonis(composition[, 4:7] ~ composition$treatment, method = "bray", permutations = 9999)
resident_adonis # treatment has significant effects on community composition: p = 1e-04 *** (df = 5)

# posthoc analysis
# A posteriori testing, using NPMANOVA, of each pair of groups can be performed after a significant result to determine this. 
treatments <- unique(composition$treatment)
treatments <- treatments[order(treatments, decreasing = FALSE)]
result <- data.frame(pair = NA, Pvalue = NA)
for (i in 1:5) {  
  for (j in (i+1):6){ # all possible pairs of comparison
    pair <- paste0(treatments[i], "-", treatments[j])
    sub <- composition %>% filter(treatment == treatments[i] | treatment == treatments[j])
    t <- adonis(sub[,4:7] ~ sub$treatment, method = "bray", permutations = 9999)
    Pvalue <- t$aov.tab$`Pr(>F)`[1]
    result <- rbind(result, cbind(pair, Pvalue))
  }
}
posthoc_flies <- result[-1,]
t <- posthoc_flies$Pvalue
names(t) <- posthoc_flies$pair
NPMANOVA.labels_flies <- data.frame(multcompLetters(t)['Letters'])
NPMANOVA.labels_flies
# flies results:
# Letters
# CTRL_CLOSED       a
# CTRL_INVAS       ab
# HW_CLOSED         c
# HW_INVAS          d
# WARM_CLOSED      ab
# WARM_INVAS        b


## two wasp species 
# no need to run NMDS (also unable) as there are only 2 dimensions in the original data. Just plot Asobara and Diapriid as two coordinates
# extract relevant data and scale the number of wasps (centered, then divided by sd)
composition <- dat_census %>%
  select(cageID, treatment, cc, Asobara, Diapriid) %>%
  mutate(Asobara.s = scale(Asobara), Diapriid.s = scale(Diapriid))

# plot - use scaled Asobara and Diapriid number to plot
find_hull <- function(df) df[chull(df$Asobara.s, df$Diapriid.s), ]
hulls <- plyr::ddply(composition, "treatment", find_hull)
NMDS_plot_wasps <- composition %>% 
  ggplot(aes(x = Asobara.s, y = Diapriid.s, group = treatment, fill = treatment, color = treatment)) + 
  geom_point(shape = 7) + geom_polygon(data = hulls, alpha = 0.5) +
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_fill_manual(values = treatment_cc, name = "treatment") + 
  xlab("Asobara") + ylab("Trichopria") + 
  ggtitle("both wasps") + 
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "right")
#ggsave("NMDS_uplandWasps.png", width = 8, height = 6)

# statistic difference (use direct nnumber to calculate community dissimilarity) - significant
# ANOSIM statistic R: 0.2735 ; Significance: 2e-04
anosim(composition[, 4:5], composition$treatment, distance = "bray", permutations = 9999) # Bray-Curtis dissimilarity is an asymmetrical measure often used for raw count data.
# NPMANOVA: p = 1e-04 *** (df = 5)
wasp_NOMANOVA <- adonis(composition[, 4:5] ~ composition$treatment, method = "bray", permutations = 9999)
wasp_NOMANOVA
# posthoc analysis
# A posteriori testing, using NPMANOVA, of each pair of groups can be performed after a significant result to determine this. 
treatments <- unique(composition$treatment)
treatments <- treatments[order(treatments, decreasing = FALSE)]
result <- data.frame(pair = NA, Pvalue = NA)
for (i in 1:5) {  
  for (j in (i+1):6){ # all possible pairs of comparison
    pair <- paste0(treatments[i], "-", treatments[j])
    sub <- composition %>% filter(treatment == treatments[i] | treatment == treatments[j])
    t <- adonis(sub[,4:5] ~ sub$treatment, method = "bray", permutations = 9999)
    Pvalue <- t$aov.tab$`Pr(>F)`[1]
    result <- rbind(result, cbind(pair, Pvalue))
  }
}
posthoc_wasp <- result[-1,]
t <- posthoc_wasp$Pvalue
names(t) <- posthoc_wasp$pair
NPMANOVA.labels_wasp <- data.frame(multcompLetters(t, threshold = 0.05)['Letters'])
NPMANOVA.labels_wasp
# wasp result
# Letters
# CTRL_CLOSED       a
# CTRL_INVAS       bc
# HW_CLOSED         b
# HW_INVAS          c
# WARM_CLOSED       b
# WARM_INVAS       ab

# Figure 3: response of host and parasitoid trophic level to temperature treatments and invasion
plot_grid(NMDS_plot_flies, NMDS_plot_wasps, nrow = 1, ncol = 2)
ggsave("figure4_community_responses.png", width = 16, height = 5)

# stat result output
sink("NMDS_results.txt")
print("========= upland flies ==========")
resident_adonis
NPMANOVA.labels_flies
print("========= wasps ==========")
wasp_NOMANOVA
NPMANOVA.labels_wasp
sink()


############## impact on individual resident species (census) ##############
# Asobara
# plot
dat_census %>%
  ggplot(aes(x = invasion, y = Asobara, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  facet_grid(. ~ temp) + xlab("Novel species (D.pandora)") + ylab("Number of Asobara wasp") + ylim(c(0,700)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5a_Asobara census.png", width = 8, height = 6)
# stat - emmeans
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
# https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#links
# dt = Inf: This is simply the way that emmeans labels asymptotic results (z test rather than t test). https://cran.r-project.org/web/packages/emmeans/vignettes/FAQs.html#asymp
mod.asobara <- glm(data = dat_census, formula = Asobara ~ temp * invasion + block, family = poisson(link = "log"))
plot(mod.asobara)
summary(mod.asobara)
# emm.asobara <- emmeans(mod.asobara, specs = pairwise~temp:invasion, type = "response")  # comparing all possible pairs - not our aims
# the effect of temperature with or without invasion
emm1.asobara <- emmeans(mod.asobara, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm1.asobara$emmeans
emm1.asobara$contrasts %>% summary(infer = TRUE)
# the effect of invasion in each temperature scenario
emm2.asobara <- emmeans(mod.asobara, specs = pairwise ~ invasion|temp, type = "response")
emm2.asobara$emmeans
emm2.asobara$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.asobara, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# Diapriid
dat_census %>%
  ggplot(aes(x = invasion, y = Diapriid, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  facet_grid(. ~ temp) + xlab("Novel species (D.pandora)") + ylab("Number of Diapriid wasp") + ylim(c(0,40)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5b_Diapriid census.png", width = 8, height = 6)

mod.diapriid <- glm(data = dat_census, formula = Diapriid ~ temp * invasion + block, family = poisson(link = "log")) 
plot(mod.diapriid)
summary(mod.diapriid)
emm1.diapriid <- emmeans(mod.diapriid, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.diapriid <- emmeans(mod.diapriid, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.diapriid, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# PST
dat_census %>%
  ggplot(aes(x = invasion, y = PST.M, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  facet_grid(. ~ temp) + xlab("Novel species (D. pandora)") + ylab("Number of D. pseudotakahashii male") + ylim(c(0,90)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5c_PST census.png", width = 8, height = 6)

# the overall effects treatments
mod.pst0 <- glm(data = dat_census, formula = PST.M ~ temp * invasion + block, family = poisson(link = "log")) 
plot(mod.pst0)
summary(mod.pst0)
emm1.pst0 <- emmeans(mod.pst0, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.pst0 <- emmeans(mod.pst0, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.pst0, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# test if the effects are mediated by other species
mod.pst <- glm(data = dat_census, formula = PST.M ~ temp * invasion + Asobara + Diapriid + block, family = poisson(link = "log")) 
plot(mod.pst)
summary(mod.pst)
emm1.pst <- emmeans(mod.pst, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.pst <- emmeans(mod.pst, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.pst, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# PAL
dat_census %>%
  ggplot(aes(x = invasion, y = PAL.M+PAL.F, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + facet_grid(. ~ temp) + xlab("Novel species (D. pandora)") + ylab("Number of D. pallidifrons") + ylim(c(0,400)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5d_PAL census.png", width = 8, height = 6)

# the overall effects treatments
mod.pal0 <- dat_census %>%
  mutate(PAL.T = PAL.M + PAL.F) %>%
  glm(formula = PAL.T ~ temp * invasion + block, family = poisson(link = "log")) 
summary(mod.pal0)
emm1.pal0 <- emmeans(mod.pal0, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.pal0 <- emmeans(mod.pal0, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.pal0, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# test indirect effects: invasion is significant when wasps numbers are not conditioned on.
mod.pal <- dat_census %>%
  mutate(PAL.T = PAL.M + PAL.F) %>%
  glm(formula = PAL.T ~ temp * invasion + Asobara + Diapriid + block, family = poisson(link = "log")) 
plot(mod.pal) # diagnosis accepted
summary(mod.pal) # invasion is no longer significant when wasps numbers are conditioned on.
emm1.pal <- emmeans(mod.pal, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.pal <- emmeans(mod.pal, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.pal, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# RUB
dat_census %>%
  ggplot(aes(x = invasion, y = RUB.M+RUB.F, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  facet_grid(. ~ temp) + xlab("Novel species (D. pandora)") + ylab("Number of D. rubida") + ylim(c(0,60)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5e_RUB census.png", width = 8, height = 6)
mod.rub <- dat_census %>%
  mutate(RUB.T = RUB.M + RUB.F) %>%
  glm(formula = RUB.T ~ temp * invasion + block, family = poisson(link = "log")) 
plot(mod.rub)
summary(mod.rub)
emm1.rub <- emmeans(mod.rub, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.rub <- emmeans(mod.rub, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.rub, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# BIR
dat_census %>%
  ggplot(aes(x = invasion, y = BIR.M, group = block)) + 
  geom_point(aes(color = treatment, shape = block), size = 1.5) + 
  stat_summary(fun=mean, geom="point", aes(group = treatment, color = treatment), size = 2) +
  stat_summary(fun.data=mean_se, geom = "errorbar", aes(group = treatment, color = treatment), size = 1, width = 0.2, alpha = 0.6) + 
  scale_color_manual(values = treatment_cc, name = "treatment") + 
  scale_shape_manual(values = c("B1" = 1, "B2" = 2), name = "block") + 
  facet_grid(. ~ temp) + xlab("Novel species (D. pandora)") + ylab("Number of D. birchii male") + ylim(c(0,20)) +
  theme(legend.position = "right") + 
  theme(strip.text.x = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times New Roman"))
ggsave("figure5f_BIR census.png", width = 8, height = 6)
mod.bir <- glm(data = dat_census, formula = BIR.M ~ temp * invasion + block, family = poisson(link = "log")) 
plot(mod.bir)
emm1.bir <- emmeans(mod.bir, specs = pairwise ~ temp|invasion, type = "response")  # back-transform the data (type = "unlink" and type = "response" are the same here. Unlink only undo the transform by the link function of the glm)
emm2.bir <- emmeans(mod.bir, specs = pairwise ~ invasion|temp, type = "response")
multcomp::cld(emmeans(mod.bir, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)

# bonferroni correction for multiple comparisons
mod1 <- as.data.frame(cbind(as.numeric(coef(summary(mod.asobara))[,'Estimate']), as.numeric(coef(summary(mod.asobara))[,'Pr(>|z|)'])))
mod2 <- as.data.frame(cbind(as.numeric(coef(summary(mod.diapriid))[,'Estimate']), as.numeric(coef(summary(mod.diapriid))[,'Pr(>|z|)'])))
mod3 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pst0))[,'Estimate']), as.numeric(coef(summary(mod.pst0))[,'Pr(>|z|)'])))
mod4 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pal0))[,'Estimate']), as.numeric(coef(summary(mod.pal0))[,'Pr(>|z|)'])))
mod5 <- as.data.frame(cbind(as.numeric(coef(summary(mod.bir))[,'Estimate']), as.numeric(coef(summary(mod.bir))[,'Pr(>|z|)'])))
mod6 <- as.data.frame(cbind(as.numeric(coef(summary(mod.rub))[,'Estimate']), as.numeric(coef(summary(mod.rub))[,'Pr(>|z|)'])))
mod.adjust <- rbind(mod1, mod2, mod3, mod4, mod5, mod6)
names(mod.adjust) <- c("coef", "pValue")
mod.adjust <- mod.adjust %>%
  mutate(pAdjust = p.adjust(pValue, method = "bonferroni")) %>%
  mutate(pAdjust = round(pAdjust, digits = 5))

# Effect of temp*invasion on ending population size - glm output and emmeans pairwise comparison
sink("census_regression and emmeans_results.txt")
print("========= Asobara ==========")
summary(mod.asobara)
emm1.asobara$contrasts %>% summary(infer = TRUE)
emm2.asobara$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.asobara, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= Diapriid ==========")
summary(mod.diapriid)
emm1.diapriid$contrasts %>% summary(infer = TRUE)
emm2.diapriid$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.diapriid, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= PST ==========")
summary(mod.pst0)
emm1.pst0$contrasts %>% summary(infer = TRUE)
emm2.pst0$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.pst0, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= PAL ==========")
summary(mod.pal0)
emm1.pal0$contrasts %>% summary(infer = TRUE)
emm2.pal0$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.pal0, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= BIR ==========")
summary(mod.bir)
emm1.bir$contrasts %>% summary(infer = TRUE)
emm2.bir$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.bir, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= RUB ==========")
summary(mod.rub)
emm1.rub$contrasts %>% summary(infer = TRUE)
emm2.rub$contrasts %>% summary(infer = TRUE)
multcomp::cld(emmeans(mod.rub, specs = pairwise ~ temp:invasion, type = "response"), alpha = 0.05, Letters = LETTERS)
print("")
print("========= adjusted p values for all ==========")
mod.adjust
sink()


############## temperal trends of individual species across treatments ##############
## plot trends of a species in different treatments on one plot
# Asobara
p2Asobara <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = Asobara)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Offspring") + ggtitle("a (Asobara)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))

# Diapriid
p2Diapriid <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = Diapriid)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Offspring") + ggtitle("b (Trichopria)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))

# PST
p2pst <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = PST.M)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Male offspring") + ggtitle("c (D. pseudotakahashii)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))


# PAL
p2pal <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = PAL.M)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Male offspring") + ggtitle("d (D. pallidifrons)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))


# RUB
p2rub <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = RUB.M)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Male offspring") + ggtitle("e (D. rubida)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))

# BIR
p2bir <- dat_sample_perCage %>% 
  filter(sID!="S0.3"&sID!="S0.1"&sID!="S0.2") %>%
  ggplot(aes(x = sID, y = BIR.M)) + geom_path(aes(group = cageID, color = treatment, linetype = invasion), alpha = 0.5) +
  stat_summary(fun=mean, geom="line", aes(group = treatment, color = treatment, linetype = invasion), size = 1) +
  xlab("Sample ID (every three weeks)") + ylab("Male offspring") + ggtitle("f (D. birchii)") + 
  facet_grid(. ~ temp) + 
  scale_color_manual(values = treatment_cc, name = "Treatment") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) + 
  theme(legend.text=element_text(size=12), legend.title = element_text(size=12), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(family = "Times"))


# Supp Figure 1. Temporal trends of individual species in different treatments
plot_grid(p2Asobara, p2Diapriid, p2pst, p2pal, p2rub, p2bir, nrow = 6, ncol = 1) 
ggsave("suppFigure1_temporal trends.png", width = 12, height = 15)


############## Inferring indirect effects ##############
## first of all, check treatments' direct effects to wasps ---- we assume that there were not indirect effects through hosts because there were always abundant hosts
mod.asobara.sample <- dat_sample_perCage %>%
  glmer(formula = Asobara ~ temp * invasion + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.asobara.sample) # residuals get smaller for bigger estimate
summary(mod.asobara.sample)
emmeans(mod.asobara.sample, specs = pairwise ~ invasion|temp, type = "response") # introduction significantly decrease pst number (overall effect)
emmeans(mod.asobara.sample, specs = pairwise ~ temp|invasion, type = "response")

mod.asobara.sample2 <- dat_sample_perCage %>%
  mutate(hostN = PAL.M + PAN.M + RUB.M + BIR.M) %>%
  glmer(formula = Asobara ~ temp * invasion + scale(hostN) + block + (1|cageID) + (1|sID), family = poisson)
summary(mod.asobara.sample2) # hostN has a negative association with Asobara number.
# the negative association can not be explained by the bottom-up effect of host, instead, it is more likely due to the top-down effect of the parasitoid!
# host number is therefore not included as preditor of wasp number

mod.diapriid.sample <- dat_sample_perCage %>%
  mutate(hostN = PAL.M + PAN.M + RUB.M + BIR.M) %>%
  glmer(formula = Diapriid ~ temp * invasion + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.diapriid.sample) # residuals get smaller for bigger estimate
summary(mod.diapriid.sample)

mod7 <- as.data.frame(cbind(as.numeric(coef(summary(mod.asobara.sample))[,'Estimate']), as.numeric(coef(summary(mod.asobara.sample))[,'Pr(>|z|)'])))
mod8 <- as.data.frame(cbind(as.numeric(coef(summary(mod.diapriid.sample))[,'Estimate']), as.numeric(coef(summary(mod.diapriid.sample))[,'Pr(>|z|)'])))
mod.sample.w.adjust <- rbind(mod7, mod8)
names(mod.sample.w.adjust) <- c("coef", "pValue")
mod.sample.w.adjust <- mod.sample.w.adjust %>%
  mutate(pAdjust = p.adjust(pValue, method = "bonferroni")) %>%
  mutate(pAdjust = round(pAdjust, digits = 4))

## focal species
## D. pandora
mod.pan.sample <- dat_sample_perCage %>%
  filter(invasion == "introduction") %>%
  glmer(formula = PAN.M ~ temp + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pan.sample)
summary(mod.pan.sample)

# warning message: Model is nearly unidentifiable: very large eigenvalue - Rescale variables?
# Therefore, I rescaled the number of the species in the updated model - Then there is no warning
mod.pan.sample_indirect <- dat_sample_perCage %>%
  filter(invasion == "introduction") %>%
  glmer(formula = PAN.M ~ temp + block + scale(PAL.M) + scale(Asobara) + scale(Diapriid)  + (1|cageID) + (1|sID), family = poisson)
plot(mod.pan.sample_indirect)
summary(mod.pan.sample_indirect)

## D. pallidifrons 
## (hypotheses: indirect effects could be mediated by wasps, not flies (we know that pallidifrons is very good competitor))
# sampling data (reproduction)
mod.pal.sample <- dat_sample_perCage %>%
  glmer(formula = PAL.M ~ temp * invasion + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pal.sample)
summary(mod.pal.sample)

# test if the effect of warming is mediated by pallidifrons and/or parasitoid
# rescale the abundance of the wasps in the predictor due to the same reason as the analysis of D. pandora
mod.pal.sample_indirect <- dat_sample_perCage %>%
  glmer(formula = PAL.M ~ temp * invasion + scale(Asobara) + scale(Diapriid) + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pal.sample_indirect)
summary(mod.pal.sample_indirect)

## D. pseudotakahashii 
## (hypotheses: indirect effects could be mediated by wasps and D. pallidifrons flies (PST is inferior competitor))
# sampling data (reproduction)
mod.pst.sample <- dat_sample_perCage %>%
  glmer(formula = PST.M ~ temp * invasion + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pst.sample) # acceptable
summary(mod.pst.sample)
emmeans(mod.pst.sample, specs = pairwise ~ invasion|temp, type = "response") # introduction significantly decrease pst number (overall effect)

# test if the effect of warming is mediated by pallidifrons and/or parasitoid
# rescale the abundance of the wasps in the predictor due to the same reason as the analysis of D. pandora
mod.pst.sample_indirect <- dat_sample_perCage %>%
  glmer(formula = PST.M ~ temp * invasion + scale(PAL.M)  + scale(Diapriid) + block + (1|cageID) + (1|sID), family = poisson)
plot(mod.pst.sample_indirect)
summary(mod.pst.sample_indirect)
emmeans(mod.pst.sample_indirect, specs = pairwise ~ invasion|temp, type = "response") # introduction has no direct effect in control and heatwave, but remained significantly negative in warming


## correct for multiple comparisons
mod9 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pan.sample))[,'Estimate']), as.numeric(coef(summary(mod.pan.sample))[,'Pr(>|z|)'])))
mod10 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pan.sample_indirect))[,'Estimate']), as.numeric(coef(summary(mod.pan.sample_indirect))[,'Pr(>|z|)'])))
mod.sample.pan.adjust <- rbind(mod9, mod10)
names(mod.sample.pan.adjust) <- c("coef", "pValue")
mod.sample.pan.adjust <- mod.sample.pan.adjust %>%
  mutate(pAdjust = p.adjust(pValue, method = "bonferroni")) %>%
  mutate(pAdjust = round(pAdjust, digits = 4))

mod11 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pal.sample))[,'Estimate']), as.numeric(coef(summary(mod.pal.sample))[,'Pr(>|z|)'])))
mod12 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pal.sample_indirect))[,'Estimate']), as.numeric(coef(summary(mod.pal.sample_indirect))[,'Pr(>|z|)'])))
mod.sample.pal.adjust <- rbind(mod11, mod12)
names(mod.sample.pal.adjust) <- c("coef", "pValue")
mod.sample.pal.adjust <- mod.sample.pal.adjust %>%
  mutate(pAdjust = p.adjust(pValue, method = "bonferroni")) %>%
  mutate(pAdjust = round(pAdjust, digits = 4))

mod13 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pst.sample))[,'Estimate']), as.numeric(coef(summary(mod.pst.sample))[,'Pr(>|z|)'])))
mod14 <- as.data.frame(cbind(as.numeric(coef(summary(mod.pst.sample_indirect))[,'Estimate']), as.numeric(coef(summary(mod.pst.sample_indirect))[,'Pr(>|z|)'])))
mod.sample.pst.adjust <- rbind(mod13, mod14)
names(mod.sample.pst.adjust) <- c("coef", "pValue")
mod.sample.pst.adjust <- mod.sample.pst.adjust %>%
  mutate(pAdjust = p.adjust(pValue, method = "bonferroni")) %>%
  mutate(pAdjust = round(pAdjust, digits = 4))

# result output - Table 1
sink("reproduction_regression and emmeans_results.txt")
print("========= Asobara ==========")
summary(mod.asobara.sample)
emmeans(mod.asobara.sample, specs = pairwise ~ invasion|temp, type = "response")
emmeans(mod.asobara.sample, specs = pairwise ~ temp|invasion, type = "response")
print("")
print("========= Diapriid ==========")
summary(mod.diapriid.sample)
emmeans(mod.diapriid.sample, specs = pairwise ~ invasion|temp, type = "response")
emmeans(mod.diapriid.sample, specs = pairwise ~ temp|invasion, type = "response")
print("")
print("========= PAN ==========")
summary(mod.pan.sample)
summary(mod.pan.sample_indirect)
mod.sample.pan.adjust
print("")
print("========= PAL ==========")
summary(mod.pal.sample)
summary(mod.pal.sample_indirect)
mod.sample.pal.adjust
print("")
print("========= PST ==========")
summary(mod.pst.sample)
emmeans(mod.pst.sample, specs = pairwise ~ invasion|temp, type = "response")
summary(mod.pst.sample_indirect)
emmeans(mod.pst.sample_indirect, specs = pairwise ~ invasion|temp, type = "response")
mod.sample.pst.adjust
sink()


############## climate record ##############
## monthly temperature summary
climate.dat <- read.csv("Data/longtermClimate_formated.csv")
climate.feb <- climate.dat %>%
  filter(Site2016 == "P780") %>% 
  filter(year_month == "2017-02") %>%
  mutate(heatwave = ifelse(Date == "13/02/2017", "YES", 
                           ifelse(Date == "14/02/2017", "YES", "NO"))) 
## summplementary figure 1a
climate.feb %>%
  ggplot(aes(x = justTime.c, y = Celsius)) +  geom_line(aes(group = day, color = heatwave), alpha = 0.4) + 
  scale_color_manual(values = c("YES" = "red", "NO" = "blue")) + 
  scale_y_continuous(limits = c(17, 32), breaks = seq(17, 32, by = 1)) + ylab("Temperature (degree Celsius)") +
  scale_x_discrete(breaks = c("00:15", "04:15", "08:15", "12:15", "16:15", "20:15", "23:15")) + xlab("Time of the day") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("upland feb heatwave.png", width = 4, height = 3)



