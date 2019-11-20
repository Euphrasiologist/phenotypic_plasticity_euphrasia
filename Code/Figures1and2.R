# Figures 1 and 2 for the manuscript #
# updated 19.11.19 #
# Max Brown #

library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(lme4)

# Figure 1 #

# need four traits, julian days to flower, corolla length, nodes to flower, height
manyh <- fread("./Data/Manyhosts.csv")[, .(Experiment = as.factor("Phenotypic plasticity experiment"),
                                            `Host or Species` = Host, 
                                           `Height (mm)` = Height, 
                                           `Nodes to flower` = Nodes_to_flower_inc_coty, 
                                           `Corolla length (mm)` = Standard_corolla_l,
                                           `Julian days to flower` = Julian_days_to_flower)]
manysp <- fread("./Data/Manyspecies.csv")[, .(Experiment = as.factor("Species differences experiment"),
                                              `Host or Species` = Expert_ID_living,
                                              `Height (mm)` = Height, 
                                              `Nodes to flower` = Nodes_to_flower_inc_coty, 
                                              `Corolla length (mm)` = Standard_corolla_l,
                                              `Julian days to flower` = Julian_days_to_flower)]

# need a data table with columns: experiment, Host or species with levels, variable, value

com <- rbind(manyh, manysp)

meltcom <- melt.data.table(data = com, id.vars = c("Experiment", "Host or Species"))

# reorder factor levels for experiment
meltcom$Experiment <- factor(x = meltcom$Experiment, levels = c("Species differences experiment",
                                          "Phenotypic plasticity experiment"))
meltcom$`Host or Species` <- as.factor(meltcom$`Host or Species`)
levels(meltcom$`Host or Species`)[17]<- "No host"
levels(meltcom$`Host or Species`)[16]<- "M. polymorpha"
levels(meltcom$`Host or Species`)[13]<- "E. arvense"
levels(meltcom$`Host or Species`)[18]<- "P. sylvestris"
levels(meltcom$`Host or Species`)[14]<- "F. rubra"
levels(meltcom$`Host or Species`)[15]<- "H. lanatus"
levels(meltcom$`Host or Species`)[1]<- "A. thaliana"
levels(meltcom$`Host or Species`)[19]<- "P. lanceolata"
levels(meltcom$`Host or Species`)[20]<- "T. repens"
levels(meltcom$`Host or Species`)[10] <- "E. confusa x nemorosa"
levels(meltcom$`Host or Species`)[12] <- "E. confusa x tetraquetra"


meltcom$`Host or Species` <- factor(meltcom$`Host or Species`, levels = c( "No host", "M. polymorpha", "E. arvense", "P. sylvestris", "F. rubra",
                                                                       "H. lanatus", "A. thaliana", "P. lanceolata", "T. repens",
                                                                       "E. arctica", "E. confusa", 
                                                                       "E. micrantha", "E. nemorosa", "E. pseudokerneri", 
                                                                       "E. anglica x nemorosa","E. anglica x rostkoviana",
                                                                       "E. arctica x confusa", 
                                                                       "E. arctica x nemorosa",
                                                                       "E. confusa x nemorosa",
                                                                       "E. confusa x tetraquetra"))


figure1 <- ggplot(meltcom, aes(x=`Host or Species`, y=value))+
  geom_point(alpha=0.4, position = position_jitter(width = 0.3),
             size = 1, aes(colour=Experiment))+
  geom_boxplot(alpha=0, col="black")+
  facet_grid(variable ~ Experiment, scales = "free", switch = "y")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, angle = 60, face = "italic"),
        legend.position = "none",
        strip.background = element_rect(fill = "white", colour="white"),
        strip.text.y = element_text(),
        strip.placement = "outside",
        panel.grid = element_blank())+
  ylab(label = "")+
  xlab(label = "")

ggsave(filename = "Figure1_updated.pdf", plot = figure1, device = "pdf", 
       path = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2016/Manuscript/AJB/V4",
       width = 7.25, height = 7.25, units = "in")

# Figure 2 #

# we plot the common vs wild data in an intuitive way
# The variables to consider are nodes to flower, corolla length, number of teeth on the lower floral leaf and internode ratio
comwild <- fread("./Data/Wildcommon.csv") # = pope4ecom

tidy.comwild <- data.table::melt.data.table(data = comwild, id.vars = c("Experiment", "E4E", "Species")) # = tidy.pop

# Nodes to flower - updated 19.11.19
tidy.comwild2 <- tidy.comwild[variable == "Nodes to flower inc coty", .(variable, value), by = .(E4E, Experiment, Species)]

sum.tidy.pop <- tidy.comwild2[, .(N=.N,
                  mean = mean(value),
                  sem = sd(value)/sqrt(.N)), by = .(E4E, Experiment, Species, variable)]


common.1 <- filter(sum.tidy.pop, sum.tidy.pop$Experiment == "Common")
wild.1 <- filter(sum.tidy.pop, sum.tidy.pop$Experiment == "Wild")

commonandwild.1<-data.frame(common.1, wild.1)

plot_colours <- data.frame(colour = c("#AEA200", "#00C1A7",
                                      "#00A6FF", "#EF67EB",
                                      "#F8766D", "#DB8E00",
                                      "#64B200", "#00BD5C",
                                      "#B385FF", "#FF63B6"),
                           Species = c("E. arctica", 
                                      "E. confusa",
                                      "E. nemorosa", 
                                      "E. pseudokerneri", 
                                      "E. anglica x nemorosa",
                                      "E. anglica x rostkoviana",
                                      "E. arctica x confusa", 
                                      "E. arctica x nemorosa",
                                      "E. confusa x nemorosa",
                                      "E. confusa x tetraquetra"))
commonandwild.1 <- merge(plot_colours, commonandwild.1, by = c("Species"))


####### Get the legend first

get.legend <- ggplot(commonandwild.1, aes(x=mean,
                                          y=mean.1, 
                                          group = colour))+
  geom_errorbar(aes(ymin = mean.1-2*sem.1, ymax = mean.1+2*sem.1), width = 0.15)+
  geom_errorbarh(aes(xmin = mean-2*sem, xmax = mean+2*sem), width=0.5)+
  geom_point(aes(colour = colour), size=3.5) +
  geom_abline(intercept = coef(lm(mean.1 ~ mean, data = commonandwild.1[-1,]))[1], slope = coef(lm(mean.1 ~ mean, data = commonandwild.1[-1,]))[2], col="red", lty = 2)+
  theme_bw()+
  labs(title = "Nodes to flower")+
  xlab(label = "")+
  ylab(label = "")+
  scale_colour_discrete(name = "Taxa",
                        breaks = c("#AEA200", "#00C1A7",
                                   "#00A6FF", "#EF67EB",
                                   "#F8766D", "#DB8E00",
                                   "#64B200", "#00BD5C",
                                   "#B385FF", "#FF63B6"),
                        labels = c("E. arctica", 
                                   "E. confusa",
                                   "E. nemorosa", 
                                   "E. pseudokerneri", 
                                   "E. anglica x nemorosa",
                                   "E. anglica x rostkoviana",
                                   "E. arctica x confusa", 
                                   "E. arctica x nemorosa",
                                   "E. nemorosa x confusa",
                                   "E. tetraquetra x confusa"))+
  theme(legend.text = element_text(size = 8, face = "italic"))

legendforerrorbars<-get_legend(get.legend)

####### Nodes to flower

plote1<-ggplot(commonandwild.1, aes(x=mean,
                                    y=mean.1))+
  geom_errorbar(aes(ymin = mean.1-2*sem.1, ymax = mean.1+2*sem.1), width = 0.15)+
  geom_errorbarh(aes(xmin = mean-2*sem, xmax = mean+2*sem), width=0.5)+
  geom_point(aes(colour = Species), size=3.5) +
  geom_abline(intercept = coef(lm(mean.1 ~ mean, data = commonandwild.1[-1,]))[1], slope = coef(lm(mean.1 ~ mean, data = commonandwild.1[-1,]))[2], col="red", lty = 2)+
  theme_bw()+
  labs(title = "")+
  xlab(label = "Nodes to flower (wild collected)")+
  ylab(label = "Nodes to flower (common garden)")+
  scale_colour_discrete(name = "Taxa",
                        breaks = c("#AEA200", "#00C1A7",
                                   "#00A6FF", "#EF67EB",
                                   "#F8766D", "#DB8E00",
                                   "#64B200", "#00BD5C",
                                   "#B385FF", "#FF63B6"),
                        labels = c("E. arctica", 
                                   "E. confusa",
                                   "E. nemorosa", 
                                   "E. pseudokerneri", 
                                   "E. anglica x nemorosa",
                                   "E. anglica x rostkoviana",
                                   "E. arctica x confusa", 
                                   "E. arctica x nemorosa",
                                   "E. nemorosa x confusa",
                                   "E. tetraquetra x confusa"))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))

# Corolla length

sum.tidy.pop2 <- tidy.comwild[variable == "Standard corolla l"][, .(N = .N,
                                                   mean = mean(value),
                                                   sem = sd(value)/sqrt(.N)), by = .(E4E, Experiment, Species, variable)]

common.2 <- filter(sum.tidy.pop2, sum.tidy.pop2$Experiment == "Common")
wild.2 <- filter(sum.tidy.pop2, sum.tidy.pop2$Experiment == "Wild")

commonandwild.2<-data.frame(common.2, wild.2)

plote2<-ggplot(commonandwild.2, aes(x=mean,
                                    y=mean.1))+
  geom_errorbar(aes(ymin = mean.1-2*sem.1, ymax = mean.1+2*sem.1), width = 0.15, size=0.5)+
  geom_errorbarh(aes(xmin = mean-2*sem, xmax = mean+2*sem), width=0.5, size=0.5)+
  geom_point(aes(colour = Species), size=3.5) +
  geom_abline(intercept = coef(lm(mean.1 ~ mean, data = commonandwild.2))[1], 
              slope = coef(lm(mean.1 ~ mean, data = commonandwild.2))[2], col="red", lty = 2)+
  theme_bw()+
  labs(title = "")+
  xlab(label = "Corolla length (wild collected)")+
  ylab(label = "Corolla length (common garden)")+
  scale_colour_discrete(name = "Taxa",
                        breaks = c("#AEA200", "#00C1A7",
                                   "#00A6FF", "#EF67EB",
                                   "#F8766D", "#DB8E00",
                                   "#64B200", "#00BD5C",
                                   "#B385FF", "#FF63B6"),
                        labels = c("E. arctica", 
                                   "E. confusa",
                                   "E. nemorosa", 
                                   "E. pseudokerneri", 
                                   "E. anglica x nemorosa",
                                   "E. anglica x rostkoviana",
                                   "E. arctica x confusa", 
                                   "E. arctica x nemorosa",
                                   "E. nemorosa x confusa",
                                   "E. tetraquetra x confusa"))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))

# Number of leaf teeth on the lower floral leaf

sum.tidy.pop3 <- tidy.comwild[variable == "No. teeth lower floral leaf"][-c(12,23,144,195)][, .(N = .N,
                                                                    mean = mean(value),
                                                                    sem = sd(value)/sqrt(.N)), by = .(E4E, Experiment, Species, variable)]

common.3 <- filter(sum.tidy.pop3, sum.tidy.pop3$Experiment == "Common")
wild.3 <- filter(sum.tidy.pop3, sum.tidy.pop3$Experiment == "Wild")

commonandwild.3 <- data.frame(common.3, wild.3)


plote3<-ggplot(commonandwild.3, aes(x=mean,
                                    y=mean.1))+
  geom_errorbar(aes(ymin = mean.1-2*sem.1, ymax = mean.1+2*sem.1), width = 0.15, size=0.5)+
  geom_errorbarh(aes(xmin = mean-2*sem, xmax = mean+2*sem), width=0.5, size=0.5)+
  geom_point(aes(colour = Species), size=3.5) +
  geom_abline(intercept = coef(lm(mean.1 ~ mean, data = na.omit(commonandwild.3)))[1], 
              slope = coef(lm(mean.1 ~ mean, data = na.omit(commonandwild.3)))[2], col="red", lty = 2)+
  theme_bw()+
  labs(title = "")+
  xlab(label = "Number of leaf teeth (wild collected)")+
  ylab(label = "Number of leaf teeth (common garden)")+
  scale_colour_discrete(name = "Taxa",
                        breaks = c("#AEA200", "#00C1A7",
                                   "#00A6FF", "#EF67EB",
                                   "#F8766D", "#DB8E00",
                                   "#64B200", "#00BD5C",
                                   "#B385FF", "#FF63B6"),
                        labels = c("E. arctica", 
                                   "E. confusa",
                                   "E. nemorosa", 
                                   "E. pseudokerneri", 
                                   "E. anglica x nemorosa",
                                   "E. anglica x rostkoviana",
                                   "E. arctica x confusa", 
                                   "E. arctica x nemorosa",
                                   "E. nemorosa x confusa",
                                   "E. tetraquetra x confusa"))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))

# Leaf ratio

sum.tidy.pop4 <- tidy.comwild[variable == "Cauline internode leaf ratio"][-c(83, 144, 23, 12)][, .(N = .N,
                                                                              mean = mean(value),
                                                                              sem = sd(value)/sqrt(.N)), by = .(E4E, Experiment, Species, variable)]

common.4 <- filter(sum.tidy.pop4, sum.tidy.pop4$Experiment == "Common")
wild.4 <- filter(sum.tidy.pop4, sum.tidy.pop4$Experiment == "Wild")

commonandwild.4 <- data.frame(common.4, wild.4)


plote4<-ggplot(commonandwild.4, aes(x=mean,
                                    y=mean.1))+
  geom_errorbar(aes(ymin = mean.1-2*sem.1, ymax = mean.1+2*sem.1), width = 0.15, size=0.5)+
  geom_errorbarh(aes(xmin = mean-2*sem, xmax = mean+2*sem), width=0.5, size=0.5)+
  geom_point(aes(colour = Species), size=3.5) +
  geom_abline(intercept = coef(lm(mean.1 ~ mean, data = na.omit(commonandwild.4)))[1], 
              slope = coef(lm(mean.1 ~ mean, data = na.omit(commonandwild.4)))[2], col="red", lty = 2)+
  theme_bw()+
  labs(title = "")+
  xlab(label = "Internode ratio (wild collected)")+
  ylab(label = "Internode ratio (common garden)")+
  scale_colour_discrete(name = "Taxa",
                        breaks = c("#AEA200", "#00C1A7",
                                   "#00A6FF", "#EF67EB",
                                   "#F8766D", "#DB8E00",
                                   "#64B200", "#00BD5C",
                                   "#B385FF", "#FF63B6"),
                        labels = c("E. arctica", 
                                   "E. confusa",
                                   "E. nemorosa", 
                                   "E. pseudokerneri", 
                                   "E. anglica x nemorosa",
                                   "E. anglica x rostkoviana",
                                   "E. arctica x confusa", 
                                   "E. arctica x nemorosa",
                                   "E. nemorosa x confusa",
                                   "E. tetraquetra x confusa"))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))

figure2<-plot_grid(plote1, plote2, plote3, plote4, labels = "auto")
# needs small amount of tweaking but basically there
figure2 <- plot_grid(figure2, legendforerrorbars, rel_widths = c(1, 0.3), scale = 0.9)

ggsave(filename = "Figure2_updated.pdf", plot = figure2, device = "pdf", 
       path = "/Users/mbrown/Dropbox/Euphrasia 2016 common garden hosts vs pops/Experiments 2016/Manuscript/AJB/V4",
       width = 7.25, units = "in")
