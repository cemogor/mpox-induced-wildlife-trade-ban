rm(list=ls())


# LOAD PACKAGES -----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(performance)                    
library(ggpubr)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(ggsci)
library(tidyr)
                           

################################################################################
################################################################################
################################################################################
###############                                                #################
###############        SENSITIVITY MODEL                       #################
###############                                                #################
################################################################################
################################################################################
################################################################################


sensitivity.model <- read.csv("01. Data/sensitivity.model.csv")


str(sensitivity.model)

sensitivity.model <- sensitivity.model %>% 
  dplyr::mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
                week = factor(week),
                location = factor(location),
                season = factor(season),
                vendor.id = factor(vendor.id),
                mass = as.numeric(summed.mass),
                lagged.mass = as.numeric(lagged.mass),
                transformed.mass = log(mass)) 

str(sensitivity.model)

# VISUALISE THE DATA ------------------------------------------------------

boxplot(sensitivity.model$transformed.mass ~ sensitivity.model$week)

boxplot(sensitivity.model$transformed.mass ~ sensitivity.model$season)
boxplot(sensitivity.model$transformed.mass ~ sensitivity.model$season*sensitivity.model$period)
boxplot(sensitivity.model$transformed.mass ~ sensitivity.model$period)
boxplot(sensitivity.model$transformed.mass ~ sensitivity.model$week * sensitivity.model$period)

hist(sensitivity.model$mass)
hist(sensitivity.model$transformed.mass)

#Range of response variable
range(sensitivity.model$transformed.mass)


# UNIVARIATE PLOT ---------------------------------------------------------

sensitivity.season <- ggplot(sensitivity.model, aes(period, mass, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.season

sensitivity.week<- ggplot(sensitivity.model, aes(week, mass, fill = period)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Julian week (2022)", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        axis.text.x = element_text(angle=90, vjust = 0.4, size =8),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.week

sensitivity.market <- ggplot(sensitivity.model, aes(location, mass, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.market

sensitivity.lagg <- ggplot(sensitivity.model %>% filter(!is.na(mass)), aes(lagged.mass, mass)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged mass (kg)", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.lagg

sensitivity.order <- ggplot(sensitivity.model %>% filter(!is.na(mass)), aes(order, mass, fill = order)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Taxonomic order", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.order


sensitivity.fig.a <- ggarrange(sensitivity.season, sensitivity.week, sensitivity.market, sensitivity.lagg,
                             labels = "AUTO",
                             ncol = 2,
                             nrow = 2,
                             common.legend = F,
                             font.label = list(size = 16))


sensitivity.fig.b <- ggarrange(sensitivity.fig.a, sensitivity.order,
                               labels = c("", "E"),
                               ncol = 1,
                               nrow = 2,
                               heights = c(1,0.5),
                               font.label = list(size = 16))

ggsave("Figure S5.png",
       sensitivity.fig.b, 
       path = "..//figure/final-figures",
       height = 24,
       width = 21,
       units = "cm",
       dpi = 300)


# RUN MODEL ---------------------------------------------------------------

sensitivity.model.out <- lmer(transformed.mass ~ period + season + location + order + scale(lagged.mass) + (1|week) + (1|vendor.id), 
                              na.action = na.exclude,
                              data = sensitivity.model)

summary(sensitivity.model.out)

# Extract fixed effects coefficients
coefs <- summary(sensitivity.model.out)$coefficients

# Convert to data frame
coefs_df <- as.data.frame(coefs)

# Save to CSV
write.csv(coefs_df, "TABLE_S7.csv", row.names = TRUE)

# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S11.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(sensitivity.model.out) 

dev.off()

r2(sensitivity.model.out) 

check_convergence(sensitivity.model.out)
check_singularity(sensitivity.model.out)
check_outliers(sensitivity.model.out)


# MODEL FIGURES -----------------------------------------------------------


# PERIOD ------------------------------------------------------------------

model.sensitivity.prediction <- ggpredict(sensitivity.model.out, c("period"), back.transform = TRUE, type = "fixed")


model.sensitivity.prediction <- model.sensitivity.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>% 
  drop_na()


# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

sensitivity.model.plot <- ggplot()+
  geom_point(data=sensitivity.model, aes(x=period, y=mass, fill = season), alpha = 0.2, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.sensitivity.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.sensitivity.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                       ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Standardised vendor sensitivity", col = "") +
  scale_color_manual(values = "#B82132") +
  scale_fill_manual(values = c("#2973B2", "#2973B2"))

sensitivity.model.plot




################################################################################
################################################################################
################################################################################
###############                                                #################
###############        SENSITIVITY MODEL (EXTENDED)            #################
###############                                                #################
################################################################################
################################################################################
################################################################################


sensitivity.model.extended <- read.csv("01. Data/sensitivity.model.extended.csv")


str(sensitivity.model.extended)

sensitivity.model.extended <- sensitivity.model.extended %>% 
  dplyr::mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
                week = factor(week),
                location = factor(location),
                season = factor(season),
                vendor.id = factor(vendor.id),
                mass = as.numeric(summed.mass),
                lagged.mass = as.numeric(lagged.mass),
                transformed.mass = log(mass)) 

str(sensitivity.model.extended)


# VISUALISE THE DATA ------------------------------------------------------

boxplot(sensitivity.model.extended$transformed.mass ~ sensitivity.model.extended$week)

boxplot(sensitivity.model.extended$transformed.mass ~ sensitivity.model.extended$season)
boxplot(sensitivity.model.extended$transformed.mass ~ sensitivity.model.extended$season*sensitivity.model.extended$period)
boxplot(sensitivity.model.extended$transformed.mass ~ sensitivity.model.extended$period)
boxplot(sensitivity.model.extended$transformed.mass ~ sensitivity.model.extended$week * sensitivity.model.extended$period)

hist(sensitivity.model.extended$mass)
hist(sensitivity.model.extended$transformed.mass)

#Range of response variable
range(sensitivity.model.extended$transformed.mass)


# UNIVARIATE PLOT ---------------------------------------------------------

sensitivity.season.ext <- ggplot(sensitivity.model.extended, aes(period, mass, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.season.ext

sensitivity.week.ext <- ggplot(sensitivity.model.extended, aes(week, mass, fill = period)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Julian week (2021-22)", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(angle=90, vjust = 0.4, size =6),
        strip.text = element_text(color = "black", size = 6)) + 
  scale_fill_aaas()

sensitivity.week.ext

sensitivity.market.ext <- ggplot(sensitivity.model.extended, aes(location, mass, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.market.ext

sensitivity.lagg.ext <- ggplot(sensitivity.model.extended %>% filter(!is.na(mass)), aes(lagged.mass, mass)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged mass (kg)", y = "Mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.lagg.ext


sensitivity.order.ext <- ggplot(sensitivity.model.extended %>% filter(!is.na(mass)), aes(order, mass, fill = order)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Taxonomic order", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

sensitivity.order.ext


sensitivity.fig.ext.a <- ggarrange(sensitivity.season.ext, sensitivity.week.ext, sensitivity.market.ext, sensitivity.lagg.ext,
                                 labels = "AUTO",
                                 ncol = 2,
                                 nrow = 2,
                                 common.legend = F,
                                 font.label = list(size = 16))

sensitivity.fig.ext.b <- ggarrange(sensitivity.fig.ext.a, sensitivity.order.ext,
                               labels = c("", "E"),
                               ncol = 1,
                               nrow = 2,
                               heights = c(1,0.5),
                               font.label = list(size = 16))

ggsave("Figure S6.png",
       sensitivity.fig.ext.b, 
       path = "..//figure/final-figures",
       height = 24,
       width = 21,
       units = "cm",
       dpi = 300)


# RUN MODEL ---------------------------------------------------------------

sensitivity.model.extended.out <- lmer(transformed.mass ~ period + season + location + order + scale(lagged.mass) + (1|week) + (1|vendor.id), 
                              na.action = na.exclude,
                              data = sensitivity.model.extended)

summary(sensitivity.model.extended.out)


coefS8 <- summary(sensitivity.model.extended.out)$coefficients

# Convert to data frame
coefS8.df <- as.data.frame(coefS8)

# Save to CSV
write.csv(coefS8.df, "TABLE_S9.csv", row.names = TRUE)

# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S12.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(sensitivity.model.extended.out) 

dev.off()

r2(sensitivity.model.extended.out) 

check_convergence(sensitivity.model.extended.out)
check_singularity(sensitivity.model.extended.out)
check_outliers(sensitivity.model.extended.out)


# MODEL FIGURES -----------------------------------------------------------


# PERIOD ------------------------------------------------------------------

model.sensitivity.prediction <- ggpredict(sensitivity.model.extended.out, c("period"), back.transform = TRUE, type = "fixed")


model.sensitivity.prediction <- model.sensitivity.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>%
  drop_na()


# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

sensitivity.model.extended.plot <- ggplot()+
  geom_point(data=sensitivity.model.extended, aes(x=period, y=mass, fill = season), alpha = 0.2, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.sensitivity.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.sensitivity.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                       ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Mass (kg)", col = "") +
  scale_color_manual(values = "#B82132") +
  scale_fill_manual(values = c("#2973B2", "#2973B2"))

sensitivity.model.extended.plot




################################################################################
################################################################################
################################################################################
###############                                                #################
###############        COUNT MODEL                             #################
###############                                                #################
################################################################################
################################################################################
################################################################################


count.model <- read.csv("01. Data/number.vendors.csv")


str(count.model)

count.model <- count.model %>% 
  dplyr::mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
                week = factor(week),
                year = factor(year),
                location = factor(location),
                season = factor(season),
                standardized.n.vendors = as.numeric(standardized.n.vendors),
                lagged.n.vendors = as.numeric(lagged.n.vendors),
                transformed.vendor.count = log(standardized.n.vendors))

str(count.model)


# VISUALISE THE DATA ------------------------------------------------------

boxplot(count.model$transformed.vendor.count ~ count.model$year)

boxplot(count.model$transformed.vendor.count ~ count.model$season)
boxplot(count.model$transformed.vendor.count ~ count.model$season*count.model$period)
boxplot(count.model$transformed.vendor.count ~ count.model$period)
boxplot(count.model$transformed.vendor.count ~ count.model$week * count.model$period)

hist(count.model$standardized.n.vendors)
hist(count.model$transformed.vendor.count)

#Range of response variable
range(count.model$transformed.vendor.count)


# UNIVARIATE PLOT ---------------------------------------------------------


vc.season <- ggplot(count.model, aes(period, standardized.n.vendors, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Standardised vendor count", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

vc.season

vc.year <- ggplot(count.model, aes(year, standardized.n.vendors, fill = year)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Year", y = "Standardised vendor count", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

vc.year

vc.market <- ggplot(count.model, aes(location, standardized.n.vendors, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Standardised vendor count", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

vc.market

vc.lagg <- ggplot(count.model %>% filter(!is.na(lagged.n.vendors)), aes(lagged.n.vendors, standardized.n.vendors)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged vendor count", y = "Standardised vendor count", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()


vc.lagg

vendorcount <- ggarrange(vc.season, vc.year, vc.market, vc.lagg,
                     labels = "AUTO",
                     ncol = 2,
                     nrow = 2,
                     common.legend = F,
                     font.label = list(size = 16))



ggsave("Figure S1.png",
       vendorcount, 
       path = "..//figure/final-figures",
       height = 18,
       width = 20,
       units = "cm",
       dpi = 300)


# RUN MODEL ---------------------------------------------------------------

count.model.out <- lmer(transformed.vendor.count ~ period + season  + year + location + scale(lagged.n.vendors) + (1|week), 
                                 na.action = na.exclude,
                                 data = count.model)

summary(count.model.out)

coefS3a <- summary(count.model.out)$coefficients

# Convert to data frame
coefS3a.df <- as.data.frame(coefS3a)

# Save to CSV
write.csv(coefS3a.df, "TABLE_S3A.csv", row.names = TRUE)

# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S7.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(count.model.out) 

dev.off()

r2(count.model.out) 

check_convergence(count.model.out)
check_singularity(count.model.out)
check_outliers(count.model.out)

# POSTHOC TEST ------------------------------------------------------------

emm_options(pbkrtest.limit = 419)

posthoc.count.model.out<- emmeans(count.model.out, ~ year) #Compare the effect of year

pairs(posthoc.count.model.out, adjust="BH") #BH is best for controlling for false discovery rate (given multiple testing)

year.count <- as.data.frame(pairs(posthoc.count.model.out, adjust="BH") )

write.csv(year.count, "TABLE_S3B.csv")


# MODEL FIGURES -----------------------------------------------------------


# PERIOD ------------------------------------------------------------------

model.count.prediction <- ggpredict(count.model.out, c("period"), back.transform = TRUE, type = "fixed")


model.count.prediction <- model.count.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>% 
  drop_na()


# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

count.model.plot <- ggplot()+
  geom_point(data=count.model, aes(x=period, y=standardized.n.vendors, fill = season), alpha = 0.2, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.count.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.count.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Standardised vendor count", col = "") +
  scale_color_manual(values = "#690B22") +
  scale_fill_manual(values = c("#A4B465", "#A4B465")) +
  geom_signif(
    data = count.model,
    mapping = aes(x=period, y=standardized.n.vendors),
    comparisons = list(c("Pre-ban", "Post-ban")), 
    map_signif_level = TRUE,
    annotation = "NS", 
    textsize = 2.5,
    vjust = -0.2,
    y_position = 8,
    tip_length = 0.02)

count.model.plot



################################################################################
################################################################################
################################################################################
###############                                                #################
###############              MASS MODEL                        #################
###############                                                #################
################################################################################
################################################################################
################################################################################

mass.model <- read.csv("01. Data/mass.model.csv")

length(unique(mass.model$vendor.id))

str(mass.model)

mass.model <- mass.model %>% 
  dplyr::mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
                week = factor(week),
                year = factor(year),
                order = factor(order),
                location = factor(location),
                season = factor(season),
                vendor.id = factor(vendor.id),
                standardised.mass = as.numeric(standardized.mass),
                lagged.effect = as.numeric(lagged.effect),
                transformed.mass = log(standardised.mass)) 


table(mass.model$location)

str(mass.model)
    


# VISUALISE THE DATA ------------------------------------------------------

boxplot(mass.model$standardised.mass ~ mass.model$year)
boxplot(mass.model$standardised.mass ~ mass.model$order)

boxplot(mass.model$standardised.mass ~ mass.model$season)
boxplot(mass.model$standardised.mass ~ mass.model$location)
boxplot(mass.model$standardised.mass ~ mass.model$season*mass.model$period)
boxplot(mass.model$standardised.mass ~ mass.model$period)
boxplot(mass.model$standardised.mass ~ mass.model$week * mass.model$period)
boxplot(mass.model$standardised.mass ~ mass.model$vendor.id * mass.model$period)


hist(mass.model$standardised.mass)
hist(mass.model$transformed.mass)

#Range of response variable
range(mass.model$standardised.mass)
range(mass.model$transformed.mass)


# UNIVARIATE PLOT ---------------------------------------------------------


mass.period <-ggplot(mass.model %>%  filter(standardised.mass < 98), aes(period, standardised.mass, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

mass.period


mass.year <- ggplot(mass.model %>% filter(standardised.mass < 98), aes(year, standardised.mass, fill = year)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Year", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

mass.year

mass.market <- ggplot(mass.model %>% filter(standardised.mass < 98), aes(location, standardised.mass, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()

mass.market

mass.lagg <- ggplot(mass.model %>% filter(!is.na(lagged.effect)) %>% filter(standardised.mass < 98), aes(lagged.effect, standardised.mass)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged standardised mass (kg)", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()


mass.lagg


mass.order <- ggplot(mass.model %>% filter(standardised.mass < 98), aes(order, standardised.mass, fill = order)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Taxonomic order", y = "Standardised mass (kg)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas()


mass.order


mass.distribution.a <- ggarrange(mass.period, mass.year, mass.market, mass.lagg, 
                     labels = "AUTO",
                     ncol = 2,
                     nrow = 2,
                     font.label = list(size = 16))


mass.distribution <- ggarrange(mass.distribution.a, mass.order,
                                 labels = c("", "E"),
                                 ncol = 1,
                                 nrow = 2,
                                 heights = c(1,0.5),
                                 font.label = list(size = 16))


mass.distribution

ggsave("Figure S2.jpg",
       mass.distribution, 
       path = "..//figure/final-figures",
       height = 24,
       width = 21,
       units = "cm",
       dpi = 300)


# RUN MODEL ---------------------------------------------------------------

mass.model.out <- lmer(transformed.mass ~ period + year + location + season + order + scale(lagged.effect) + (1|vendor.id) + (1|week), 
                       na.action = na.exclude,
                       data = mass.model)

summary(mass.model.out)


coefS4a <- summary(mass.model.out)$coefficients

# Convert to data frame
coefS4a.df <- as.data.frame(coefS4a)

# Save to CSV
write.csv(coefS4a.df, "TABLE_S4A.csv", row.names = TRUE)

# POSTHOC TEST ------------------------------------------------------------

emm_options(pbkrtest.limit = 3961)

posthoc.mass.model.model <- emmeans(mass.model.out, ~ year) #Compare the effect of year
posthoc.mass.model.model2 <- emmeans(mass.model.out, ~ order, pbkrtest.limit = 3939)


pairs(posthoc.mass.model.model, adjust="BH") #BH is best for controlling for false discovery rate (given multiple testing)
pairs(posthoc.mass.model.model2, adjust="BH") 

coef4b.df <- as.data.frame(pairs(posthoc.mass.model.model, adjust="BH") )
coef4c.df <- as.data.frame(pairs(posthoc.mass.model.model2, adjust="BH") )


write.csv(coef4b.df, "TABLE_S4B.csv")
write.csv(coef4c.df, "TABLE_S4C.csv")



# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S8.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(mass.model.out) #Model fit

dev.off()

r2(mass.model.out) #Variance explained 

check_convergence(mass.model.out)
check_singularity(mass.model.out)
check_outliers(mass.model.out)

# MODEL FIGURES -----------------------------------------------------------


# PERIOD ------------------------------------------------------------------

model.mass.prediction <- ggpredict(mass.model.out, c("period"), back.transform = TRUE, type = "fixed")


model.mass.prediction <- model.mass.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>% 
  drop_na()


# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

mass.model.plot <- ggplot()+
  geom_point(data=mass.model, aes(x=period, y=standardized.mass, fill = season), alpha = 0.2, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.mass.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.mass.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                 ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Standardised mass (kg)", col = "") +
  scale_color_manual(values = "#690B22") +
  scale_fill_manual(values = c("#A4B465", "#A4B465")) +
  geom_signif(
    data = mass.model,
    mapping = aes(x=period, y=standardized.mass),
    comparisons = list(c("Pre-ban", "Post-ban")), 
    map_signif_level = TRUE,
    annotation = "NS", 
    textsize = 2.5,
    vjust = -0.2,
    y_position = 180,
    tip_length = 0.02)
  
mass.model.plot



################################################################################
################################################################################
################################################################################
###############                                                #################
###############             PRICE MODEL                        #################
###############                                                #################
################################################################################
################################################################################
################################################################################

income.model <- read.csv("01. Data/income.model.csv")

length(unique(income.model$vendor.id))


str(income.model)

income.model <- income.model %>% 
  mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
         week = factor(week),
         year = factor(year),
         location = factor(location),
         season = factor(season),
         vendor.id = factor (vendor.id),
         lagged.effect = as.numeric(lagged.effect),
         standardized.income = as.numeric(standardized.income),
         transformed.income = log(standardized.income)) 

str(income.model)



# VISUALISE THE DATA ------------------------------------------------------
boxplot(income.model$standardized.income ~ income.model$year)
boxplot(income.model$standardized.income ~ income.model$season)
boxplot(income.model$standardized.income ~ income.model$period)
boxplot(income.model$standardized.income ~ income.model$week * income.model$period)
boxplot(income.model$standardized.income ~ income.model$season * income.model$period)
boxplot(income.model$standardized.income ~ income.model$vendor.id * income.model$period)


hist(income.model$standardized.income)
hist(income.model$transformed.income)

#Range of response variable
range(income.model$standardized.income)
range(income.model$transformed.income)


# UNIVARIATE PLOT ---------------------------------------------------------


income.period <-ggplot(income.model %>%  filter(standardized.income < 300000), aes(period, standardized.income, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Standardised price (₦)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


income.period


income.year <- ggplot(income.model %>% filter(standardized.income < 300000), aes(year, standardized.income, fill = year)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Year", y = "Standardised price (₦)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


income.year

income.market <- ggplot(income.model %>% filter(standardized.income < 300000), aes(location, standardized.income, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Standardised price (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


income.market

income.lagg <- ggplot(income.model %>% filter(!is.na(lagged.effect)) %>% 
                        filter(standardized.income < 300000) %>% 
                        filter(lagged.effect < 300000), aes(lagged.effect, standardized.income)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged standardised price (₦)", y = "Standardised price (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k")) +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  



income.lagg


income.order <- ggplot(income.model %>% filter(standardized.income < 300000), aes(order, standardized.income, fill = order)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Taxonomic order", y = "Standardised price (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  



income.order


income.distribution.a <- ggarrange(income.period, income.year, income.market, income.lagg, 
                                 labels = "AUTO",
                                 ncol = 2,
                                 nrow = 2,
                                 font.label = list(size = 16))


income.distribution <- ggarrange(income.distribution.a, income.order,
                               labels = c("", "E"),
                               ncol = 1,
                               nrow = 2,
                               heights = c(1,0.5),
                               font.label = list(size = 16))


income.distribution

ggsave("Figure S3.jpg",
       income.distribution, 
       path = "..//figure/final-figures",
       height = 24,
       width = 21,
       units = "cm",
       dpi = 300)


# RUN MODEL ---------------------------------------------------------------

income.model.out <- lmer(transformed.income ~ period + year + location + season + order + scale(lagged.effect) + (1|vendor.id) + (1|week), 
                         na.action = na.exclude,
                         data = income.model)

summary(income.model.out)


coefS5a <- summary(income.model.out)$coefficients

# Convert to data frame
coefS5a.df <- as.data.frame(coefS5a)

# Save to CSV
write.csv(coefS5a.df, "TABLE_S5A.csv", row.names = TRUE)

# POSTHOC TEST ------------------------------------------------------------


emm_options(pbkrtest.limit = 3961)

posthoc.income.model.model <- emmeans(income.model.out, ~ year) #Compare the effect of year
posthoc.income.model.model2 <- emmeans(income.model.out, ~ order)

pairs(posthoc.income.model.model, adjust="BH") #BH is best for controlling for false discovery rate (given multiple testing)
pairs(posthoc.income.model.model2, adjust="BH") 

coefS5b.df <- as.data.frame(pairs(posthoc.income.model.model, adjust="BH") )
coefS5c.df <- as.data.frame(pairs(posthoc.income.model.model2, adjust="BH") )

write.csv(coefS5b.df, "TABLE_S5B.csv")
write.csv(coefS5c.df, "TABLE_S5C.csv")


# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S9.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(income.model.out) #Model fit

dev.off()

r2(income.model.out) #Variance explained 

check_convergence(income.model.out)
check_singularity(income.model.out)
check_outliers(income.model.out)


# PERIOD ------------------------------------------------------------------

model.income.prediction <- ggpredict(income.model.out, c("period"), back.transform = TRUE, type = "fixed")


model.income.prediction <- model.income.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>% 
  drop_na()



# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

income.model.plot <- ggplot()+
  geom_point(data=income.model, aes(x=period, y=standardized.income, fill = season), alpha = 0.2, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.income.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.income.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Standardised price (₦)", col = "") +
  scale_color_manual(values = "#690B22") +
  scale_fill_manual(values = c("#A4B465", "#A4B465")) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-5, suffix = "M")) +
  geom_signif(
    data = income.model,
    mapping = aes(x=period, y=standardized.income),
    comparisons = list(c("Pre-ban", "Post-ban")), 
    map_signif_level = TRUE,
    annotation = "*", 
    textsize = 5,
    vjust = 0.5,
    y_position = 820000,
    tip_length = 0.02)
  

income.model.plot



################################################################################
################################################################################
################################################################################
###############                                                #################
###############              PPKG MODEL                        #################
###############                                                #################
################################################################################
################################################################################
################################################################################

ppkg.model <- read.csv("01. Data/ppkg.model.csv")

length(unique(ppkg.model$vendor.id))


str(ppkg.model)

ppkg.model <- ppkg.model %>% 
  mutate(period = factor(period, levels = c("Pre-ban", "Post-ban")),
         week = factor(week),
         year = factor(year),
         location = factor(location),
         season = factor(season),
         vendor.id = factor (vendor.id),
         lagged.effect = as.numeric(lagged.effect),
         standardized.ppkg = as.numeric(standardized.ppkg),
         transformed.ppkg = log(standardized.ppkg)) 
str(ppkg.model)


# VISUALISE THE DATA ------------------------------------------------------
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$year)
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$season)
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$period)
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$week * ppkg.model$period)
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$season * ppkg.model$period)
boxplot(ppkg.model$standardized.ppkg ~ ppkg.model$vendor.id * ppkg.model$period)


hist(ppkg.model$standardized.ppkg)
hist(ppkg.model$transformed.ppkg)

#Range of response variable
range(ppkg.model$standardized.ppkg)
range(ppkg.model$transformed.ppkg)


# UNIVARIATE PLOT ---------------------------------------------------------

ppkg.period <-ggplot(ppkg.model %>%  filter(standardized.ppkg < 60000), aes(period, standardized.ppkg, fill = season)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Period", y = "Standardised price/kg (₦)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


ppkg.period


ppkg.year <- ggplot(ppkg.model %>% filter(standardized.ppkg < 60000), aes(year, standardized.ppkg, fill = year)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8",size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Year", y = "Standardised price/kg (₦)", fill = " ") +
  theme_pubr(base_size = 11) +
  theme(axis.title = element_text(size = 14),
        legend.position = "NULL",
        legend.text = element_text(size = 10), 
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


ppkg.year

ppkg.market <- ggplot(ppkg.model %>% filter(standardized.ppkg < 60000), aes(location, standardized.ppkg, fill = location)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Market", y = "Standardised price/kg (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  


ppkg.market

ppkg.lagg <- ggplot(ppkg.model %>% filter(!is.na(lagged.effect)) %>% 
                        filter(standardized.ppkg < 60000) %>% 
                        filter(lagged.effect < 60000), aes(lagged.effect, standardized.ppkg)) +
  geom_point(color = "#626da8", shape = 21, size=2.5, fill = "red")+
  labs(x = "Lagged standardised price/kg (₦)", y = "Standardised price/kg (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_x_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  



ppkg.lagg


ppkg.order <- ggplot(ppkg.model %>% filter(standardized.ppkg < 60000), aes(order, standardized.ppkg, fill = order)) +
  geom_boxplot(position = position_dodge(width = 0.8),color="black",  fill = "#626da8", size = .5, alpha = .8, outlier.size = 3, outlier.fill = "grey", outlier.alpha = .4) +
  labs(x = "Taxonomic order", y = "Standardised price/kg (₦)", fill = " ") +
  theme_pubr(base_size = 10) +
  theme(axis.title = element_text(size = 14),
        legend.text = element_text(size = 10), 
        legend.position = "NULL",
        strip.text = element_text(color = "black", size = 11)) + 
  scale_fill_aaas() +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k"))  



ppkg.order


ppkg.distribution.a <- ggarrange(ppkg.period, ppkg.year, ppkg.market, ppkg.lagg, 
                                   labels = "AUTO",
                                   ncol = 2,
                                   nrow = 2,
                                   font.label = list(size = 16))


ppkg.distribution <- ggarrange(ppkg.distribution.a, ppkg.order,
                                 labels = c("", "E"),
                                 ncol = 1,
                                 nrow = 2,
                                 heights = c(1,0.5),
                                 font.label = list(size = 16))


ppkg.distribution

ggsave("Figure S4.jpg",
       ppkg.distribution, 
       path = "..//figure/final-figures",
       height = 25,
       width = 21,
       units = "cm",
       dpi = 300)

# RUN MODEL ---------------------------------------------------------------

ppkg.model.out <- lmer(transformed.ppkg ~ period + year + location + season + order + scale(lagged.effect) + (1|vendor.id) + (1|week), 
                         na.action = na.exclude,
                         data = ppkg.model)

summary(ppkg.model.out)

coefS6a <- summary(ppkg.model.out)$coefficients

# Convert to data frame
coefS6a.df <- as.data.frame(coefS6a)

# Save to CSV
write.csv(coefS6a.df, "TABLE_S6A.csv", row.names = TRUE)

# POSTHOC TEST ------------------------------------------------------------


emm_options(pbkrtest.limit = 3961)

posthoc.ppkg.model.model <- emmeans(ppkg.model.out, ~ year) #Compare the effect of year
posthoc.ppkg.model.model2 <- emmeans(ppkg.model.out, ~ order)

pairs(posthoc.ppkg.model.model, adjust="BH") #BH is best for controlling for false discovery rate (given multiple testing)
pairs(posthoc.ppkg.model.model2, adjust="BH") 

coefS6b.df <- as.data.frame(pairs(posthoc.ppkg.model.model, adjust="BH") )
coefS6c.df <- as.data.frame(pairs(posthoc.ppkg.model.model2, adjust="BH") )

write.csv(coefS6b.df, "TABLE_S6B.csv")
write.csv(coefS6c.df, "TABLE_S6C.csv")


# MODEL DIAGNOSTICS -------------------------------------------------------

tiff(filename = "..//figure/final-figures/Figure S10.tiff", width = 30, height = 37, units = "cm", res = 600)

check_model(ppkg.model.out) #Model fit

dev.off()

r2(ppkg.model.out) #Variance explained 

check_convergence(ppkg.model.out)
check_singularity(ppkg.model.out)
check_outliers(ppkg.model.out)


# MODEL FIGURES -----------------------------------------------------------


# PERIOD ------------------------------------------------------------------

model.ppkg.prediction <- ggpredict(ppkg.model.out, c("period"), back.transform = TRUE, type = "fixed")


model.ppkg.prediction <- model.ppkg.prediction %>% 
  mutate(
    conf.high = exp(conf.high),
    conf.low = exp(conf.low),
    predicted = exp(predicted),
    std.error = exp(std.error)) %>% 
  drop_na()


# Define the dodge position adjustment
dodge <- position_dodge(width = 0.2)

ppkg.model.plot <- ggplot()+
  geom_point(data=ppkg.model, aes(x=period, y=standardized.ppkg, fill = season), alpha = 0.01, shape = 21, size=4, show.legend = FALSE) + 
  geom_point(data = model.ppkg.prediction, aes(x=x,y=predicted, col = group), size = 2, position = dodge) +
  geom_errorbar(data=model.ppkg.prediction, aes(x=x,y=predicted, ymin = conf.low, 
                                                ymax = conf.high, col = factor(group)), width =.1, size =1, position = dodge) +
  theme_pubr(base_size = 11) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.position = "NULL") +
  labs(x ="Period", y= "Standardised price/kg (₦)", col = "") +
  scale_color_manual(values = "#690B22") +
  scale_fill_manual(values = c("#A4B465", "#A4B465")) +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "k")) +
  geom_signif(
    data = ppkg.model,
    mapping = aes(x=period, y=standardized.ppkg),
    comparisons = list(c("Pre-ban", "Post-ban")), 
    map_signif_level = TRUE,
    annotation = "NS", 
    textsize = 2.5,
    vjust = -0.2,
    y_position = 65000,
    tip_length = 0.02)


ppkg.model.plot




# COMBINE PLOTS -----------------------------------------------------------


figure4 <- ggarrange(count.model.plot, mass.model.plot, income.model.plot, ppkg.model.plot,
                         labels = "AUTO",
                         nrow = 2, 
                         widths = c(1, 1.1, 1, 1.1),
                         ncol = 2,
                         legend = "none",
                         align = "hv")

figure4


ggsave("Figure 4.png",
       figure4, 
       path = "..//figure/final-figures",
       height = 16,
       width = 16,
       units = "cm",
       dpi = 300)
