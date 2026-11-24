###############################################################
# Experiment 1
# Maze
#
# This script reproduces all preprocessing, figures,
# descriptive statistics, and mixed-effects analyses
# reported for Experiment 1 in the manuscript.
###############################################################

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

# Specify columns and read in the data, and the list of critical regions
cols = c("time","subj","type","square","null","cond","item","wordnumber","word", "foil","Value","correct","rt","fulltext")
data = read.csv('Exp1Results.csv', header = 1, sep = ',', comment.char = '#', col.names = cols)

# Format Critical parts of sub-experiment
data$correct = as.character(ifelse(data$correct == "yes", 1,0))
data$item <- as.numeric(data$item)
# Extract data for trials by subexpt
exp1.data = subset(data, substr(cond, 0,4) == "exp.") 

#calculate error rate:
# proportion incorrect
error_rate <- mean(exp1.data$correct == "0")

#keep only correct responses
exp1.data <- exp1.data %>%
  filter(correct == "1")

# Convert hash to subject number 
exp1.data <- exp1.data %>%
  mutate(subj = paste(subj, time),
         subj = factor(subj, levels = unique(subj),
                       labels = 1:length(unique(subj))))

# change classifications of columns so we can more easily work with them
exp1.data$time <- as.numeric(as.character(exp1.data$time))
exp1.data$subj <- as.character(exp1.data$subj)
exp1.data$rt <- as.numeric(as.character(exp1.data$rt))
exp1.data$cond <- as.factor(exp1.data$cond)
exp1.data$wordnumber <- as.numeric(as.character(exp1.data$wordnumber))
#take everything off the label except condition letter
exp1.data$cond = substring(exp1.data$cond, 5) 

exp1.data <- exp1.data %>% 
  # just to make sure we are only looking at the four conditions
  filter(cond %in% c("a","b","c","d")) %>% 
  # assign contrasts (scaled sum coded)
  mutate(PronounMatch_c = ifelse(cond %in% c("a","d"), -1/2, 1/2), 
         PronounMatch = ifelse(cond %in% c("a","d"), "Match", "Mismatch"))

# remove outliers using standard outlier removal for maze
  
exp1.fastslow <- exp1.data %>%
  filter((rt<100 | rt>10000)) # drops unusually fast (<100 ms) or unusually slow (>10,000 ms) RTs
exp1.data <- exp1.data %>%
  filter(!(rt<100 | rt>10000)) #this resulted in no observations removed at the critical region

# log transform RTs
exp1.data <- exp1.data %>%
  mutate(LogRT = log(rt))

##########################
# Descriptive Statistics #
##########################

# sse summary tables
reflexive_sse <- summarySE(data = subset(exp1.data, wordnumber==9), # has reflexive region only
                           measurevar = "rt", 
                           groupvars = c("PronounMatch"),
                           na.rm = TRUE)

reflexive_sse_log <- summarySE(data = subset(exp1.data, wordnumber==9), # has reflexive region only
                     measurevar = "LogRT", 
                     groupvars = c("PronounMatch"),
                     na.rm = TRUE)



all_sse <- summarySE(data = exp1.data, # has every word region
                         measurevar = "rt", 
                         groupvars = c("PronounMatch",
                                       "wordnumber"),
                         na.rm = TRUE)

###### Plots #########

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
  ggplot(data, aes(x = PronounMatch, y = mean, color = PronounMatch)) +
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
    labs(x = "Adjective Reflexive Agreement", y = ylab) +
    plots_theme
}

plot_reflexive_log <- function(data, ylab = "Mean log RT") {
  ggplot(data, aes(x = PronounMatch, y = mean, color = PronounMatch)) +
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
    labs(x = "Adjective Reflexive Agreement", y = ylab) +
    plots_theme
}

p <- plot_reflexive(reflexive_sse)
p_log <- plot_reflexive_log(reflexive_sse_log)

p
p_log

# Full sentence
p_sentence <- ggplot(all_sse, aes(x=wordnumber, y=mean, color=PronounMatch)) + 
              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
              geom_point()+
              geom_line()+
              scale_color_manual(values = c("Match" = "#E76F51", 
                                            "Mismatch" = "#6A4C93")) +
        labs(title="Experiment 1: Full Sentence",x="Region", y = "Mean RT (ms)")+
        scale_y_continuous(breaks = c(0:1500*100)) +
        scale_x_continuous(breaks = c(0:13*1),
                     labels = c("Apparently","the,","pretty/handsome","artisan","gained,",
                                "many","clients","after", "dedicating","herself/himself","to","fixing","products","efficiently.")) +
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

p_sentence

########################
# Statistical Modeling #
########################

exp1.mdl <- lmer(LogRT ~ PronounMatch_c + (1 | subj) + (1 | item), data = subset(exp1.data, wordnumber==9), 
                       REML=F,control = lmerControl(optimizer = "bobyqa"))
summary(exp1.mdl) 


