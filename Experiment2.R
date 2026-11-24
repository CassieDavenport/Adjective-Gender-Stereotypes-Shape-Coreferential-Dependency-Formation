##########
# Set-Up #
##########

# Remove pre-loaded variables/references for a clean slate
rm(list=ls())
# Load necessary packages 
library(plyr)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(ggsignif)
library(lattice)
library(tidyverse)
library(grid)

################
# SSE Function #
################
# summarySE for mean, sd, se, and ci
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

######################
# Data Pre-Processing#
######################

# read in data

data = read.csv('Exp2Results.csv')

# Format Critical parts of sub-experiment
data$correct = as.character(ifelse(data$correct == "yes", 1,0))
data$item <- as.numeric(data$item)
# Extract data for trials by subexpt
exp2.data = subset(data, substr(cond, 0,4) == "exp.") 

# Convert hash to subject number 
exp2.data <- exp2.data %>%
  mutate(subj = paste(subj, time),
         subj = factor(subj, levels = unique(subj),
                       labels = 1:length(unique(subj))))

# change classifications of columns so we can more easily work with them
exp2.data$time <- as.numeric(as.character(exp2.data$time))
exp2.data$subj <- as.character(exp2.data$subj)
exp2.data$rt <- as.numeric(as.character(exp2.data$rt))
exp2.data$cond <- as.factor(exp2.data$cond)
exp2.data$wordnumber <- as.numeric(as.character(exp2.data$wordnumber))
#take everything off the label except condition letter
exp2.data$cond = substring(exp2.data$cond, 5) 

#calculate error rate:
# proportion incorrect
error_rate <- mean(exp2.data$correct == "0")

#keep only correct responses
exp2.data <- exp2.data %>%
  filter(correct == "1")

exp2.data <- exp2.data %>% 
  # just to make sure we are only looking at the four conditions
  filter(cond %in% c("a","b","c","d","e","f","g","h")) %>% 
  # assign contrasts (scaled sum coded)
  mutate(PronounMatch = ifelse(cond %in% c("a","c","e","g"), "Match", "Mismatch"), 
         StructureType = dplyr::case_when(
           cond %in% c("a","b") ~ "Attri",
           cond %in% c("c","d") ~ "ARC",
           cond %in% c("e","f") ~ "RRC",
           cond %in% c("g","h") ~ "Pred"
         ),
  )

# remove outliers using standard outlier removal for maze

exp2.fastslow <- exp2.data %>%
  filter((rt<100 | rt>10000)) # drops unusually fast (<100 ms) or unusually slow (>10,000 ms) RTs

exp2.data <- exp2.data %>%
  filter(!(rt<100 | rt>10000)) #this resulted in 1 observation removed at the critical region


# log transform RTs
exp2.data <- exp2.data %>%
  mutate(LogRT = log(rt))

# assign the critical region

Region = "Reflexive"

exp2.data <- exp2.data %>%
  mutate(
    Region = case_when(
      word %in% c("herself", "himself") ~ "Reflexive",
      TRUE ~ "Other"
    )
  )


##########################
# Descriptive Statistics #
##########################

# sse summary tables
reflexive_sse <- summarySE(data = subset(exp2.data, Region == "Reflexive"),
                     measurevar = "rt", 
                     groupvars = c("PronounMatch","StructureType"),
                     na.rm = TRUE)

reflexive_sse_log <- summarySE(data = subset(exp2.data, Region == "Reflexive"),
                     measurevar = "LogRT", 
                     groupvars = c("PronounMatch","StructureType"),
                     na.rm = TRUE)


all_sse <- summarySE(data = exp2.data,
                     measurevar = "rt", 
                     groupvars = c("PronounMatch","StructureType","wordnumber"),
                     na.rm = TRUE)


# plots

reflexive_sse <- reflexive_sse %>%
  mutate(StructureType = factor(StructureType, 
                                levels = c("Attri","ARC","RRC","Pred")))

reflexive_sse_log <- reflexive_sse_log %>%
  mutate(StructureType = factor(StructureType, 
                                levels = c("Attri","ARC","RRC","Pred")))

plots_theme <- theme_classic(base_size = 30) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.5,
      linetype = "solid"
    ),
    
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

plot_reflexive <- function(data, ylab = "Mean RT (ms)") {
  ggplot(data, aes(x = StructureType, y = mean, color = PronounMatch)) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.12,
      linewidth = 1.5,
      position = position_dodge(width = 0.35)
    ) +
    geom_point(
      size = 2.6,
      position = position_dodge(width = 0.35)
    ) +
    scale_color_manual(values = c("Match" = "#E76F51", 
                                  "Mismatch" = "#6A4C93")) +
    labs(x = "Adjective Information Structure", y = ylab) +
    plots_theme
}

plot_reflexive_log <- function(data, ylab = "Mean log RT") {
  ggplot(data, aes(x = StructureType, y = mean, color = PronounMatch)) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.12,
      linewidth = 0.6,
      position = position_dodge(width = 0.35)
    ) +
    geom_point(
      size = 2.6,
      position = position_dodge(width = 0.35)
    ) +
    scale_color_manual(values = c("Match" = "#E76F51", 
                                  "Mismatch" = "#6A4C93")) +
    labs(x = "Adjective Information Structure", y = ylab) +
    plots_theme
}

p <- plot_reflexive(reflexive_sse)
p_log <- plot_reflexive_log(reflexive_sse_log)

p
p_log

#full sentence plots

RRC_plot <- ggplot(subset(all_sse, StructureType=='RRC'), aes(x=wordnumber, y=mean, color=PronounMatch)) + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
    geom_point()+
    geom_line()+
    scale_color_manual(values = c("Match" = "#E76F51", 
                                  "Mismatch" = "#6A4C93")) +
    labs(title="Restrictive Relative Clause",x="Region", y = "Mean RT (ms)")+
    scale_y_continuous(breaks = c(0:1500*100)) +
    scale_x_continuous(breaks = c(0:14*1),
                       labels = c("The","cashier","who","was","handsome",
                                  "earned","a","promotion","after", "committing","himself/herself","to","excellent","customer","service.")) +
    theme(axis.text.x = element_text(size=9, angle =90),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_line(
            color = "grey90",
            linewidth = 0.5,
            linetype = "solid"
          ),
          
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()
    )

ARC_plot <- ggplot(subset(all_sse, StructureType=='ARC'), aes(x=wordnumber, y=mean, color=PronounMatch)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("Match" = "#E76F51", 
                                "Mismatch" = "#6A4C93")) +
  labs(title="Appositive Relative Clause",x="Region", y = "Mean RT (ms)")+
  scale_y_continuous(breaks = c(0:1500*100)) +
  scale_x_continuous(breaks = c(0:14*1),
                     labels = c("The","cashier,","who","was","handsome,",
                                "earned","a","promotion","after", "committing","himself/herself","to","excellent","customer","service."))+
  theme(axis.text.x = element_text(size=9, angle =90),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(
          color = "grey90",
          linewidth = 0.5,
          linetype = "solid"
        ),
        
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
  )

Attri_plot <- ggplot(subset(all_sse, StructureType=='Attri'), aes(x=wordnumber, y=mean, color=PronounMatch)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("Match" = "#E76F51", 
                                "Mismatch" = "#6A4C93")) +
  labs(title="Attributive Position",x="Region", y = "Mean RT (ms)")+
  scale_y_continuous(breaks = c(0:1500*100)) +
  scale_x_continuous(breaks = c(0:12*1),
                     labels = c("The","handsome","cashier",
                                "earned","a","promotion","after", "committing","himself/herself","to","excellent","customer","service.")) +
  theme(axis.text.x = element_text(size=9, angle =90),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(
          color = "grey90",
          linewidth = 0.5,
          linetype = "solid"
        ),
        
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
  )

Pred_plot <- ggplot(subset(all_sse, StructureType=='Pred'), aes(x=wordnumber, y=mean, color=PronounMatch)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("Match" = "#E76F51", 
                                "Mismatch" = "#6A4C93")) +
  labs(title="Predicative Position",x="Region", y = "Mean RT (ms)")+
  scale_y_continuous(breaks = c(0:1500*100)) +
  scale_x_continuous(breaks = c(0:14*1),
                     labels = c("The","cashier",
                                "was","handsome","and", "earned", "a", "promotion","after", "committing","himself/herself","to","excellent","customer","service.")) +
  theme(axis.text.x = element_text(size=9, angle =90),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(
          color = "grey90",
          linewidth = 0.5,
          linetype = "solid"
        ),
        
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()
  )

Attri_plot
ARC_plot
RRC_plot
Pred_plot

########################
# Statistical Modeling #
########################

# Pronoun Match: sum-coded
exp2.data$PronounMatch <- factor(exp2.data$PronounMatch, levels = c("Match", "Mismatch"))
contrasts(exp2.data$PronounMatch) <- c(-0.5, 0.5)

# Infor Structure: forward Helmter, Attri < ARC < RRC < Pred
# make sure level order is correct before assigning contrasts
exp2.data$StructureType <- factor(exp2.data$StructureType, levels = c("Attri", "ARC", "RRC", "Pred"))

contrasts(exp2.data$StructureType) <- cbind(
  Attri_vs_rest   = c(-3/4,  1/4,  1/4, 1/4),
  ARC_vs_RRCPred = c(   0, -2/3,  1/3, 1/3),
  RRC_vs_Pred    = c(   0,    0, -1/2, 1/2)
)

exp2.mdl <- lmer(LogRT ~ PronounMatch * StructureType + (1|subj) + (1|item),
                 data = subset(exp2.data, Region=="Reflexive"), 
                 REML=F,control = lmerControl(optimizer = "bobyqa"))

summary(exp2.mdl)

##############################################################
# Exploratory: Comparing aesthetic and behavioral adjectives #
##############################################################

# 1. Create a small dataframe for adjective types
adj_types <- data.frame(
  item = 1:32,
  adj = c("petite","fashionable","dainty","graceful","girly","stunning","chic","artistic",
          "gossipy","bewitching","curvy","catty","sassy","adorable","beautiful","pretty",
          "chivalrous","boastful","strapping","lanky","violent","dashing","burly","charming",
          "studly","powerful","tough","commanding","chiseled","hunky","handsome","buff"),
  AdjType = c("aes","aes","aes","beh","beh","aes","beh","beh","beh","beh","aes","beh",
              "beh","aes","aes","aes","beh","beh","aes","aes","beh","beh","aes","beh",
              "aes","beh","beh","beh","aes","aes","aes","aes")
)

# 2. Join with your main dataframe
exp2.data <- exp2.data %>%
  left_join(adj_types %>% select(item, adj, AdjType), by = "item")

adj_type_sse <- summarySE(data = subset(exp2.data, Region == "Reflexive"),
                          measurevar = "rt", 
                          groupvars = c("PronounMatch","StructureType","AdjType"),
                          na.rm = TRUE)

plot_reflexive_adj <- function(data, ylab = "Mean RT (ms)") {
  ggplot(data, aes(x = StructureType, y = mean, color = PronounMatch)) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.12,
      linewidth = 1.5,
      position = position_dodge(width = 0.35)
    ) +
    geom_point(
      size = 2.6,
      position = position_dodge(width = 0.35)
    ) +
    facet_wrap(~AdjType) +
    scale_color_manual(values = c("Match" = "#E76F51", 
                                  "Mismatch" = "#6A4C93")) +
    labs(x = "Adjective Information Structure", y = ylab) +
    plots_theme
}

p_adj <- plot_reflexive_adj(adj_type_sse)

p_adj

# statistical analysis

exp2.data$AdjType <- factor(exp2.data$AdjType, levels = c("aes", "beh"))
contrasts(exp2.data$AdjType) <- c(-0.5, 0.5)

exp2adj.mdl <- lmer(LogRT ~ PronounMatch * StructureType * AdjType + (1|subj) + (1|item),
                 data = subset(exp2.data, Region=="Reflexive"), 
                 REML=F,control = lmerControl(optimizer = "bobyqa"))

summary(exp2adj.mdl)





