
####################################

# Mixed Assemblage Artifical Reefs
# Stats & Figures 
# Author: Róisín Hayden 
# Date: January 2024

##################################

#---------- ADMIN ------------------

library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(gridExtra)
library(ggfortify)
library(car)
library(performance)
library(qqplotr)
library(see)
library(tidyverse)
library(broom)
library(plotrix)
library(rcompanion)
library(glmmTMB)
library(DHARMa)
library(Matrix)
library(emmeans)
library(effects)
library(fitdistrplus) 
library(corrplot)
library(viridis)

mytheme <- theme(axis.text=element_text(size=14), #levels
                  axis.title=element_text(size=18), #titles
                  strip.text.x = element_text(size = 16),
                  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
                  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())

mytheme2 <- theme(legend.position = "right", 
                  legend.background = element_rect(fill = "white", colour = "black"),
                  legend.title = element_text(size=14),
                  legend.text = element_text(size=14),
                  axis.text=element_text(size=14), #levels
                  axis.title=element_text(size=18), #titles
                  strip.text.x = element_text(size = 16),
                  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
                  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#---------- IMPORT DATA ------------------

recruits = read.csv("recruitment.csv", header=TRUE)
recruits = na.omit(recruits)
recruits$Assemblage <- factor(recruits$Assemblage, levels = c("P", "PAT", "PT", "PA", "T", "AT", "A", "C"))
recruits$Family <- as.factor(recruits$Famille)
levels(recruits$Family) <- c("Other", "Pocilliporidae")
str(recruits)

recruits2 = read.csv("percentages.csv", header=TRUE)
str(recruits2)

#---------- STATS ------------------

# Response variable: Number of recruits 
# Explanatory variable(s): Assemblage 
# Random effect: Site 

# Calculate total recruits per assemblage, site and structure
recruits_count <- recruits %>% 
  group_by(Assemblage, Site, Structure) %>% 
  summarise(total_recruits = sum(Nombre_recrues))

# Calculate total Pocillo recruits per assemblage, site and structure 

poc = subset(recruits, Family=="Pocilliporidae")

poc_count <- poc %>% 
  group_by(Assemblage, Site, Structure, Nombre_Polype) %>% 
  summarise(total_recruits = sum(Nombre_recrues))

poc_sum <- poc_count %>% 
  group_by(Assemblage) %>% 
  summarise(avg = mean(total_recruits),
            ste = std.error(total_recruits))

# Checking distribution of response variable 

plotdist(recruits_count$total_recruits, hist= TRUE, demp = TRUE) 
descdist(recruits_count$total_recruits, boot = 1000) # Strong positive skew 

# Testing common positive-skew distributions 

fit.gamma <- fitdist(recruits_count$total_recruits, "gamma")
plot(fit.gamma)
fit.lnorm <- fitdist(recruits_count$total_recruits, "lnorm")
plot(fit.lnorm)
fit.poisson <- fitdist(recruits_count$total_recruits, "pois")
plot(fit.poisson)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot.legend <- c("gamma", "lognormal", "poisson")
denscomp(list(fit.gamma, fit.lnorm, fit.poisson), legendtext = plot.legend)
qqcomp(list(fit.gamma, fit.lnorm, fit.poisson), legendtext = plot.legend)
cdfcomp(list(fit.gamma, fit.lnorm, fit.poisson), legendtext = plot.legend)
ppcomp(list(fit.gamma, fit.lnorm, fit.poisson), legendtext = plot.legend)

gofstat(list(fit.gamma, fit.lnorm, fit.poisson), 
        fitnames = c("gamma", "lnorm", "poisson")) # gamma distribution best fits response variable

# Model 1 - All Recruits

model1 = glmmTMB(total_recruits ~ Assemblage + (1|Site), data=recruits_count, family=Gamma(link = "log"))

# Testing residuals 
model1_res = simulateResiduals(fittedModel = model1, plot = T)
plot(allEffects(model1))
check_model(model1)

# Model summary 
summary(model1)
emmeans(model1, list(pairwise ~ Assemblage), adjust = "tukey")


# Model 2 - Pocillo Recruits *Only* 

model2 = glmmTMB(total_recruits ~ Assemblage + (1|Site), data=poc_count, family=Gamma(link=log))
model2_res = simulateResiduals(fittedModel = model2, plot = T)
plot(allEffects(model2))
check_model(model2)

# Model summery 
summary(model2)
emmeans(model2, list(pairwise ~ Assemblage), adjust = "tukey")


#---------- FIGURES ------------------

# Revised Figure 2 - All Recruits 

# Box plot
fig2_box = ggplot(recruits_count,aes(x=Assemblage, y=total_recruits, fill=Assemblage)) + 
  stat_boxplot(geom='errorbar', linetype=1, width=0.2) + geom_boxplot() +
  theme_bw() + ylab("Number of Recruits") +
  scale_fill_viridis(option="D", discrete = TRUE) + 
  theme(legend.position = "none") + mytheme

# Dot plot
# Calculate average # of recruits per assemblage (+ standard error)
recruit_avg <- recruits_count %>% 
  group_by(Assemblage) %>% 
  summarise(avg = mean(total_recruits),
            ste = std.error(total_recruits))

fig2_dot <- ggplot(recruit_avg)+ 
  geom_point(aes(x=Assemblage, y=avg, colour=Assemblage), position=position_dodge(.8), size=5) + 
  geom_errorbar(aes(x=Assemblage, y=avg, ymin=avg-ste, ymax=avg+ste, colour=Assemblage), width=.05,
                position=position_dodge(.8)) + theme_bw() + ylab("Number of Recruits") +
  scale_colour_viridis(option="D", discrete=TRUE) + theme(legend.position = "none") + mytheme


fig2_box | fig2_dot 


# Revised Figure 2 - Pocillo only

fig2_box_poc = ggplot(poc_count,aes(x=Assemblage, y=total_recruits, fill=Assemblage)) + 
  stat_boxplot(geom='errorbar', linetype=1, width=0.2) + geom_boxplot() +
  theme_bw() + ylab("Number of Recruits") +
  scale_fill_viridis(option="D", discrete = TRUE) + 
  theme(legend.position = "none") + mytheme

fig2_box_poc


# Revised Figure 3 
# Calculate total recruits per Assemblage and Family
recruits_count2 <- recruits %>% 
  group_by(Assemblage, Family) %>% 
  summarise(total_recruits = sum(Nombre_recrues))

# Calculated % of each outside R so import data 
recruits2 = read.csv("percentages.csv", header=TRUE)
str(recruits2)

# Arrange them in descending order 
recruits2$Assemblage <- factor(recruits2$Assemblage, levels = c("PA", "PT", "P", "PAT", "C", "AT", "T", "A"))


fig3 <- ggplot(data=recruits2, aes(x=Assemblage, y=percetage , fill=Family)) +
  geom_bar(stat="identity") + ylab("Total Number of Recruits (%)") +
  scale_fill_manual(values = c(Other = "#98B5CF", Pocilliporidae = "#2F5C84")) +
  theme_bw() + mytheme2

fig3


# Figure 4 Revised 
fig4 <- ggplot(poc, aes(poc$Nombre_Polype)) +
  geom_histogram(binwidth=2, colour="black", fill="#2F5C84") +
  ylim(0, 20) + xlim(0,80) +
  ylab("Count") + theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_blank())

# Calculate number of recruits per size (# corallites)
poc_count <- poc %>% 
  group_by(Nombre_Polype) %>% 
  summarise(total_recruits = sum(Nombre_recrues))

fig4.1= ggplot() +
  geom_point(data=poc_count, aes(x=Nombre_Polype, y=total_recruits), size = 2) +
  geom_smooth(data=poc_count, aes(x=Nombre_Polype, y=total_recruits),
              method = lm, colour="#2F5C84", size=1.5) +
  ylab("Number of Recruits") +
  xlab("Size of Recruit (No. of Corallites)") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

fig4 | fig4.1

require(gridExtra)
grid.arrange(fig4, fig4.1, nrow=2)

#---------- EXPORT FIGURES ------------------


png("figure2_box.png", width = 8, height = 8, units = 'in', res = 300)
fig2_box
dev.off()

png("figure2_dot.png", width = 8, height = 8, units = 'in', res = 300)
fig2_dot 
dev.off()

png("figure3_revised.png", width = 10, height = 8, units = 'in', res = 300)
fig3
dev.off()

png("figure4_revised.png", width = 10, height = 10, units = 'in', res = 300)
grid.arrange(fig4, fig4.1, nrow=2)
dev.off()

