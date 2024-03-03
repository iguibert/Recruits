
####################################

# Effect of Coral-Giant Clam Artificial Reefs on Coral Recruitment: 
# Insights for Restoration and Conservation Efforts
# Figures & Statistics (Revised)  
# R Script Author: Róisín Hayden 
# Date: January 2024

##################################

#---------- ADMIN ------------------

library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(gridExtra)
library(ggfortify)
library(performance)
library(qqplotr)
library(see)
library(tidyverse)
library(broom)
library(plotrix)
library(rcompanion)
library(viridis)
library(readxl)
library(MCMCglmm)

mytheme <- theme(axis.text=element_text(size=14), #levels
                 axis.title=element_text(size=18), #titles
                 strip.text.x = element_text(size = 16),
                 axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
                 axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())

mytheme2 <- theme(legend.position = c(0.15, 0.875), 
                  legend.background = element_rect(fill = "white", colour = "black"),
                  legend.title = element_text(size=14),
                  legend.text = element_text(size=14),
                  axis.text=element_text(size=14), #levels
                  axis.title=element_text(size=18), #titles
                  strip.text.x = element_text(size = 16),
                  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
                  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#---------- IMPORT & ORGANISE DATA ------------------

data <- read_excel("Supplementary_S2_Guibert_et_al_recruits.xlsx",
                   na = "NA")
data$Site <- as.factor(data$Site)
str(data)

###### Variables

recruit <- data[,c(1,4,5)]
recruit$Assemblage[grep("C",recruit$Assemblage)] <- "C"

recruit_tile <- recruit %>% group_by(Site,Assemblage,Tile_name) %>%
  summarise(Nbr=n())



recruit2 <- data[,c(1,4,5,6)]
recruit2$Assemblage[grep("C",recruit2$Assemblage)] <- "C"

recruit_tile2 <- recruit2 %>% group_by(Site,Assemblage,Tile_name, Group) %>%
  summarise(Nbr=n())

recruit3 <- data[,c(1,4,5,6, 10)]
recruit3 = na.omit(recrut3)
recruit3$Assemblage[grep("C",recruit3$Assemblage)] <- "C"

recruit_tile3 <- recruit3 %>% group_by(Site,Assemblage,Tile_name, Family) %>%
  summarise(Nbr=n())

recruit_count <- recruit_tile3 %>% group_by(Assemblage, Family) %>%
  summarise(avg = mean(Nbr), 
            sd = sd(Nbr))

recruit_count=na.omit(recruit_count)
write.csv(recruit_count, file = "recruits.csv", quote = FALSE, row.names = F)


recruit_count2 <- recruit_tile3 %>% group_by(Assemblage) %>%
  summarise(avg = mean(Nbr), 
            sd = sd(Nbr))


### Proportion of Recruit Types
percent <- read.csv("Supplementary_S3_Guibert_et_al_recruits.csv",
                    na = "NA")
percent$Site <- as.factor(data2$Site)

### Sup. Data: Light & Temperature
env_data = read_excel("Supplementary_S4_Guibert_et_al_recruits.xlsx")
env_data <- env_data %>% filter(intensity_lux > 0)



#----------  STATISTICS ------------------


### ANOVA: (Sup. Data) Light & Temp 

light_aov = aov(log(intensity_lux) ~ site, data=env_data)
par(mfrow=c(2,2))
plot(light_aov)
summary(light_aov)
TukeyHSD(light_aov)

temp_aov = aov(temp ~ site, data=env_data)
par(mfrow=c(2,2))
plot(temp_aov)
summary(temp_aov)
TukeyHSD(temp_aov)


### Total Number of Recruits 
#### SITE EFFECT

#ON THE NUMBER OF RECRUITS
Site_Recruit <- data %>% filter(!is.na(Family)) %>% count(Site,.drop=F)

x<-chisq.test(Site_Recruit$n,simulate.p.value = F)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-paste("Site",Site_Recruit$Site)
Tab


#ON THE NUMBER OF Pocilloporidae RECRUITS
Site_Poc <- data %>% filter(Family=="PD") %>% count(Site,.drop=F)
x<-chisq.test(Site_Poc$n,simulate.p.value = F)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-paste("Site",Site_Poc$Site)
Tab

#ON THE NUMBER OF Other RECRUITS
Site_Other <- data %>% filter(Family=="OTHER") %>% count(Site,.drop=F)
x<-chisq.test(Site_Other$n,simulate.p.value = F)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-paste("Site",Site_Other$Site)
Tab

#ON THE size and shape of structures
Site_Structure <- data %>% filter(!is.na(Family)) %>%
  count(Site,Artificial_reef_shape,.drop=F) %>%
  pivot_wider(names_from = Site, values_from = n) %>%
  column_to_rownames(var= "Artificial_reef_shape")

x<-chisq.test(Site_Structure,simulate.p.value = F)
x

Tab<-round(x$observed,2)
row.names(Tab)<-paste("Structure",row.names(Tab),sep="-")
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Obs"
Tab

Tab<-round(x$expected,2)
row.names(Tab)<-paste("Structure",row.names(Tab),sep="-")
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Expected"
Tab

Tab<-round(x$stdres,2)
row.names(Tab)<-paste("Structure",row.names(Tab),sep="-")
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Stdres"
Tab

#ON THE Assemblages
Site_Assemblage <- data %>% filter(!is.na(Family)) %>%
  count(Site,Assemblage,.drop=F) %>%
  pivot_wider(names_from = Site, values_from = n) %>%
  column_to_rownames(var= "Assemblage")

x<-chisq.test(Site_Assemblage,simulate.p.value = T, B=5000)
x

Tab<-round(x$observed,2)
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Obs"
Tab

Tab<-round(x$expected,2)
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Expected"
Tab

Tab<-round(x$stdres,2)
colnames(Tab)<-paste("Site",colnames(Tab),sep="-")
"Stdres"
Tab

###############################
####
#### ASSEMBLAGE EFFECT

#All recruits
Assemblage_Recruit <- data %>% filter(!is.na(Family)) %>%
  count(Assemblage,.drop=F) %>%
  column_to_rownames(var= "Assemblage")

x<-chisq.test(Assemblage_Recruit$n,simulate.p.value = F)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-rownames(Assemblage_Recruit)
Tab

# Assemblage recruit Family
Assemblage_family <- data %>% filter(!is.na(Family)) %>%
  group_by(Site,Assemblage,Family) %>%
  summarise(nb=n(), .groups="keep") %>%
  pivot_wider(names_from = Family, values_from = nb) %>%
  replace_na(list(PD=0,OTHER=0))


#PD recruits
Assemblage_PD <- data %>% filter(Family=="PD") %>%
  count(Assemblage,.drop=F) %>%
  column_to_rownames(var= "Assemblage")

x<-chisq.test(Assemblage_PD$n,simulate.p.value = F)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-rownames(Assemblage_PD)
Tab

#Other recruits
Assemblage_other <- data %>% filter(Family=="OTHER") %>%
  count(Site,Assemblage,.drop=F) %>%
  column_to_rownames(var= "Assemblage")

x<-chisq.test(Assemblage_other$n,simulate.p.value = T, B=5000)
x
Tab<-round(rbind(x$observed,x$expected,x$stdres),2)
row.names(Tab)<-c('Obs','The','Std')
colnames(Tab)<-rownames(Assemblage_other)
Tab

# Family - Assemblage
res <- chisq.test(table(data$Assemblage,Data$Family),simulate.p.value = T)
res

### GLMM: Number of Recruits per Group (Containing P. acuta or Not)
# random: site

grp <- MCMCglmm(Nbr~Group-1, random=~Site,data=recruit_tile2,
                nitt=2000, burnin=100, thin = 2, family="poisson")
summary(grp)


#----------  MAIN FIGURES ------------------

### FIGURE 2 (Panel)
### Number of Recruits per Tile by: 

### 1) ASSEMBLAGE

recruit_avg <- recruit_tile %>% group_by(Assemblage) %>%
  summarise(avg = mean(Nbr), 
            std = std.error(Nbr))

# Re-order points in descending order 
recruit_avg$Assemblage <- factor(recruit_avg$Assemblage, levels = c("PAT", "PT","P", "PA", "T", "C","AT", "A"))

fig2.1 <- ggplot(recruit_avg)+ 
  geom_point(aes(x=Assemblage, y=avg, colour=Assemblage), position=position_dodge(.8), size=5) + 
  geom_errorbar(aes(x=Assemblage, y=avg, ymin=avg-std, ymax=avg+std, colour=Assemblage), width=.05,
                position=position_dodge(.8)) + theme_bw() + ylim(0, 6) +
  scale_colour_viridis(option="D", discrete = TRUE) + 
  ylab("Number of Recruits per Tile") + theme(legend.position = "none") + mytheme

fig2.1

## 2) GROUP - Assemblages Containing P. acuta (P) or Not Containing P. acuta (NP)

recruit_tile2$Group <- factor(recruit_tile2$Group, levels = c("P", "NP"))

recruit_avg2 <- recruit_tile2 %>% group_by(Group) %>%
  summarise(avg = mean(Nbr), 
            std = std.error(Nbr))

fig2.2 <- ggplot(recruit_avg2)+ 
  geom_point(aes(x=Group, y=avg, colour=Group), position=position_dodge(.8), size=5) + 
  geom_errorbar(aes(x=Group, y=avg, ymin=avg-std, ymax=avg+std, colour=Group), width=.05,
                position=position_dodge(.8)) + theme_bw() + ylim(0, 6) +
  scale_colour_manual(values = c(NP = "#98B5CF", P = "#2F5C84")) +
  scale_x_discrete(name="Assemblage", labels =c("P. acuta present", "P. acuta absent")) +
  mytheme2 + theme(legend.position = "none", axis.title.y = element_blank())

fig2.2

## Create panel figure

require(gridExtra)
grid.arrange(fig2.1, fig2.2, ncol=2)

### FIGURE 3
### Proportion of Recruit Type - Pocilliporidae or Other - by Assemblage: 

# Re-order points in descending order 
percent$Assemblage <- factor(percent$Assemblage, levels = c("PA", "PT","P", "C", "PAT", "AT","T", "A"))

fig3 <- ggplot(data=percent, aes(x=Assemblage, y=percetage, fill=Family)) +
  geom_bar(stat="identity") +
  ylab("Relative Abundance of Recruit Type (%)") + 
  scale_fill_manual(values = c(Other = "#98B5CF", Pocilliporidae = "#2F5C84")) +
  theme_bw() + mytheme2 + theme(legend.position = "right") 

fig3

### FIGURE 4 
### Distribution of Recruit Sizes 

# 1) Histogram 

PD_live = subset(data, Family=="PD" & Health=="Alive")

fig4.1 <- ggplot(PD_live, aes(PD_live$Corallite_number)) +
  geom_histogram(binwidth=2, colour="black", fill="#2F5C84") +
  ylim(0, 20) + xlim(0,80) +
  ylab("Count") + theme_bw() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x = element_blank())

fig4.1

# 2) Dot plot with exponential curve 

surv <- data %>% dplyr::select(Family, Corallite_number, Health) %>%
  filter(Health=="Alive" & Family=="PD") %>% group_by(Corallite_number) %>%
  summarise(Nbr=n())

model <- lm(log(surv$Nbr)~surv$Corallite_number)

plot(surv$Corallite_number, surv$Nbr,
     ylab = "Number of Recruits", xlab = "Size of Recruit", ylim = c(0,10), xlim=c(0,80), cex.lab = 1.45)
lines(surv$Corallite_number,
      exp(surv$Corallite_number*model$coefficients[2])*exp(model$coefficients[1]),lwd = 3, type="l",col="#2F5C84")

fig4.2 <- recordPlot()
plot.new()
fig4.2


#----------  SUPPLEMENTARY FIGURES ------------------

### Sup. Figure 1: Number of Recruits per Tile by site  

recruit_avg3 <- recruit_tile3 %>% group_by(Site) %>%
  summarise(avg = mean(Nbr), 
            std = std.error(Nbr))

supfig1 <- ggplot(recruit_avg3)+ 
  geom_point(aes(x=Site, y=avg, colour=Site), position=position_dodge(.8), size=5) + 
  geom_errorbar(aes(x=Site, y=avg, ymin=avg-std, ymax=avg+std, colour=Site), width=.05,
                position=position_dodge(.8)) + theme_bw() + ylim(0, 6) +
  scale_colour_viridis(option="D", discrete = TRUE) +
  ylab("Number of Recruits per Tile") + theme(legend.position = "none") + mytheme

supfig1

### Sup. Figure 2: Light intensity and temperature per Site  

env_sum <- env_data %>%
  group_by(site) %>%
  summarise(avg_light = mean(intensity_lux),
            std_light = std.error(intensity_lux), 
            avg_temp = mean(temp), 
            std_temp = std.error(temp))

supfig2.1 <- ggplot(env_sum)+ 
  geom_errorbar(aes(x=site, y=avg_light, ymin=avg_light-std_light, ymax=avg_light+std_light), width=.05,
                position=position_dodge(.8)) + 
  geom_point(aes(x=site, y=avg_light), position=position_dodge(.8), size=5, colour="goldenrod1")+
  ylab("Light Intensity (lx)") + xlab("Site") +
  theme_bw() + theme(legend.position = "none") + mytheme

supfig2.1

supfig2.2 <- ggplot(env_sum)+ 
  geom_errorbar(aes(x=site, y=avg_temp, ymin=avg_temp-std_temp, ymax=avg_temp+std_temp), width=.2,
                position=position_dodge(.8)) + ylim(28, 30) + 
  geom_point(aes(x=site, y=avg_temp, colour=site), position=position_dodge(.8), size=5, colour="darkorange1") + 
  ylab("Average Temperature (°C)") + xlab("Site") +
  theme_bw() + theme(legend.position = "none") + mytheme 

supfig2.2 

#----------  EXPORT HIGH-QUALITY FIGURES ------------------

png("FIG2_final_panel.png", width = 16, height = 10, units = 'in', res = 300)
grid.arrange(fig2.1, fig2.2, ncol=2)
dev.off()

png("FIG3_final.png", width = 12, height = 10, units = 'in', res = 300)
type_bar
dev.off()

png("FIG4_dot.png", width = 12, height = 10, units = 'in', res = 300)
fig4.1
dev.off()

png("FIG4_hist.png", width = 12, height = 10, units = 'in', res = 300)
fig4.2
dev.off()

png("SUPPFIG_site.png", width = 16, height = 10, units = 'in', res = 300)
supfig1
dev.off()

png("SUPPFIG_light_temp.png", width = 16, height = 10, units = 'in', res = 300)
grid.arrange(supfig2.1, supfig2.2, ncol=2)
dev.off()

