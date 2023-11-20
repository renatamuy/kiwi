# Code for evaluating optimal condictions for Emeiria growth in kiwi poo

#-------------------------

library(lme4)
library(boot)
require(xlsx)
library(glmmTMB)
library(ggplot2)
require(bbmle)

set.seed(123)

df <- xlsx::read.xlsx('D://OneDrive - Massey University//kiwi/kiwi_eimeria.xlsx', sheetIndex = 2)

head(df)

# Fit the model using the beta family + adding 0.001 to all values so the interval falls
# within probabilities interval

df$Total_Sporulation_transformed <- df$Total.sporulated.oocysts /100 + 0.001

df$Total_Sporulation_transformed <- ifelse(df$Total_Sporulation_transformed > 1, 0.99, df$Total_Sporulation_transformed) 

# Model fresh versus old 

df$sample_age <- factor(df$Rep, levels = c(1, 2, 3), labels = c("Old", "Fresh", "Fresh"))

# Rep instead of age

model_sample_rep <- glmmTMB( Total_Sporulation_transformed ~ Temperature + Light + Rep +      (1 + Day),
                                     data = df,
                                     family = beta_family(link = "logit"))   

# Model with sample and age

model_sample_age <- glmmTMB( Total_Sporulation_transformed ~ Temperature + Light + sample_age +      (1 + Day),
                                     data = df,
                                     family = beta_family(link = "logit"))        
  


model_beta_fixed_sample_age <- glmmTMB(Total_Sporulation_transformed ~ Day +
                                         Temperature + 
                                         Light +
                              sample_age, data = df, family = beta_family(link = "logit"))


# Full model 

model_full <- glmmTMB( Total_Sporulation_transformed ~ Temperature + Light + Rep +  sample_age + (1 + Day),
                             data = df,
                             family = beta_family(link = "logit"))  
# Interaction

model_beta_fixed_sample_age_int <- glmmTMB(Total_Sporulation_transformed ~ Day*Temperature + 
                                             Light +
                                         sample_age, data = df, 
                                         family = beta_family(link = "logit"))


summary(model_beta_fixed_sample_age_int)

#

model_beta_fixed <- glmmTMB(Total_Sporulation_transformed ~ Day + Temperature + Light +
                         Rep, data = df, family = beta_family(link = "logit"))


summary(model_beta_fixed)

model_beta_random <- glmmTMB(Total_Sporulation_transformed ~ Day + Temperature + Light +
                                      (1 | Rep), data = df, family = beta_family(link = "logit"))
  
# Calculate confidence intervals
conf_intervals <- round(confint(model_beta_random)[2:4,], 3)

conf_intervals
confint(model_beta_fixed_sample_age)
confint(model_beta_fixed_sample_age_int)

#Just a few things, I tested a some models for the % Total Sporulation as response variable after taking a better look:
#The good news is that the full models are meaningful and superior to a no effect model (null hypothesis) :) 
#After testing models with fixed and random effects, the best model in terms of fit is a fixed effects model that considers all predictors but I included a new variable replacing Reps: instead of Reps (1,2,3) the model considers the age of the sample (so Rep1 becomes ‘Older’, and Reps 2 and 3 becomes ‘Fresh’). This model was slightly better than a model considering Reps as fixed variable.

marginal_effects <- effects::effect("Day:Temperature:Light", model_beta_random)

marginal_effects

plot(marginal_effects, intervals = TRUE)

#We allowed for random intercepts for each level of the "Rep" variable as there were differences between ages of the sample especially for the Rep = 1 data, which came from older samples while 2 and 3 were fresh samples. This way, the model managed to capture unobserved variability or correlation within each group.

# Get variance component estimates
var_components <- VarCorr(model_beta_random)
print(var_components)

# Residuals are equally spread all around groups

dotchart(resid(model_beta_random, type = "re"), groups = df$Rep, main = "Residuals by Group")

simulated_residuals <- simulateResiduals(model_beta_random)

plot(simulated_residuals) # Significant deviations for residuals

# Intercept-only model

no_effect_model <- glmmTMB(Total_Sporulation_transformed ~ 1 ,
                           data = df, family = beta_family(link = "logit"))

# Model selection -exploring

AICctab(no_effect_model,
        model_full,
        model_sample_age, # Total_Sporulation_transformed ~ Day + Temperature + Light + sample_age +      (1 + Day)
        model_sample_rep,
         model_beta_fixed_sample_age_int, # Total_Sporulation_transformed ~ Day * Temperature + Light + sample_age
        #model_beta_fixed, #Total_Sporulation_transformed ~ Day + Temperature + Light + Rep
        #model_beta_fixed_sample_age, 
        logLik=FALSE, 
        weights = TRUE )

AIC(no_effect_model) - AIC(model_beta)

## FINAL Models with the time autocorrelation --------------------
# Model with time structure for day

df <- df[order(df$Day), ]

# Note: Check time dependencies for first data point for reps

df$Total_Sporulation_lag1 <- c(NA, head(df$Total_Sporulation_transformed, -1))

# Remove the first NA lag so model with time is comparable to NO effect model

dfc <- df[2:nrow(df),]

dfc$Total_Sporulation_lag1

# Interaction model
# Does UV have an effect on sporulation rate and is there a correlation between UV influence and temperature? I think it might be the higher the temperature, the more influence UV has.

model_autocor_UVtemp_interaction <- glmmTMB(Total_Sporulation_transformed ~ Temperature*Light + 
                                       Rep + Total_Sporulation_lag1 + (1 | Day),
                                       data = dfc, family = beta_family(link = "logit"))

# Fit the model with lagged values to capture autocorrelation
model_autocor <- glmmTMB(Total_Sporulation_transformed ~ Temperature + Light + Rep + Total_Sporulation_lag1 + (1 | Day),
  data = dfc, family = beta_family(link = "logit") )

model_autocor_age <- glmmTMB(Total_Sporulation_transformed ~ Temperature + Light + sample_age + Total_Sporulation_lag1 + (1 | Day),
                         data = dfc, family = beta_family(link = "logit") )

summary(model_autocor_age)

# Print the summary of the model
summary(model_autocor)

# No effect model

model_autocor_no_effect <-  glmmTMB(
  Total_Sporulation_transformed ~ 1 + (1 | Day),
  data = dfc,
  family = beta_family(link = "logit")
)


#  Model comparison

summary(model_autocor)

summary(model_autocor_UVtemp_interaction)

AICctab(model_autocor_no_effect, 
        model_autocor_age,
        model_autocor_UVtemp_interaction,
        model_autocor, logLik=FALSE, 
        weights = TRUE )

summary(model_autocor_no_effect)

# Confint best model

round(confint(model_autocor_age), 3)

#glht(my.glm, mcp(Temperature="Tukey"))

plot(repeated_measures_model)
anova(repeated_measures_model)

sjPlot::plot_model(repeated_measures_model, type = "re", data = df)

sjPlot::plot_model(repeated_measures_model, show.values = TRUE, sort.est = FALSE)

model <- aov(df$Circular.sporulation.. ~ Rep * Temperature * Light * Day, data = df)

model_beta_fixed_sample_age

hist(residuals(model_beta_fixed_sample_age))

# ANOVA DIAGNOSIS ------------------

summary(model)

hist(residuals(model))

shapiro.test(residuals(model)) # not normal 

# References

shapiro.test(rnorm(100, mean = 5, sd = 3)) # normal

shapiro.test(runif(100, min = 2, max = 4)) # not normal 


# So , since it breaks assumption of normality, we can use proper glmms

## Exploring
data <- df

interaction.plot(x.factor = data$Day, 
                 trace.factor = data$Temperature, response = data$Circular.sporulation.., 
                 type = "b", fixed = TRUE, legend = TRUE, 
                 col = c("red", "blue", 'pink', 'orange', 'black'))

# Display the summary of the model

summary(repeated_measures_model)

interaction.plot(x.factor = data$Day, trace.factor =data$Light ,
                 response = data$Circular.sporulation.., type = "b", 
                 fixed = TRUE, legend = TRUE, col = c( "green", 
                                                       "purple"))

# Morphotypes

ggplot(data, aes(x = Day, y = Circular.sporulation.., 
                 color = Light)) +
  facet_grid(~Rep + Temperature) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(aes(group = Light), position = position_dodge(width = 0.3), 
            linetype = "dashed") +
  labs(title = "Interaction Plot of Day, Temperature, and Light on Circular Sporulation",
       x = "Day",
       y = "Circular Sporulation") +
  scale_color_manual(values = c("Dark" = "purple", "Light" = "orange")) +
  theme_minimal()


# Exploring 

ggplot(data, aes(x = Day, y = Teardrops.sporulation.., 
                 color = Light)) +
  facet_grid(~Rep + Temperature) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(aes(group = Light), position = position_dodge(width = 0.3), 
            linetype = "dashed") +
  labs(title = "Interaction Plot of Day, Temperature, and Light on Teardrop Sporulation",
       x = "Day",
       y = "Teardrop Sporulation") +
  scale_color_manual(values = c("Dark" = "purple", "Light" = "orange")) +
  theme_minimal()

# Total sporulation is the new response variable
head(data)
df$Total.sporulated.oocysts
splot <- ggplot(df, aes(x = Day, y = Total.sporulated.oocysts, 
                 color = Light)) +
  facet_grid(Rep  ~ Temperature, labeller = label_both) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(aes(group = Light), position = position_dodge(width = 0.3), 
            linetype = "dashed") +
  labs(title = "Interaction Plot of Day, Temperature, and Light on Total Sporulation",
       x = "Day",
       y = "Sporulation (%)") +
  scale_color_manual(values = c("Dark" = "purple", "Light" = "orange")) +
  theme_minimal()


setwd('D:/OneDrive - Massey University/kiwi/')

png(filename = 'Total_sporulation.png', res = 400, units = 'cm', width = 20, height = 15      )
splot
dev.off()



#-----
# Create a data frame for predictions
new_data <- expand.grid(
  Day = seq(min(df$Day), max(df$Day), length.out = 100),
  Temperature = median(df$Temperature),  # Adjust this based on your data
  Light = levels(df$Light)  # Assuming Light is a factor
)

# Convert Light to a factor if it's not already
new_data$Light <- as.factor(new_data$Light)

# Add random effects to the predictions
predictions <- predict(repeated_measures_model, newdata = new_data, re.form = NA)
new_data$predicted <- predictions

# Plot the results
ggplot(df, aes(x = Day, y = Circular.sporulation.., color = Light)) +
  geom_point() +
  geom_line(data = new_data, aes(x = Day,
                                 y = predicted, group = Light, 
                                 linetype = "dashed")) +
  labs(title = "Mixed-Effects Model Predictions",
       x = "Day",
       y = "Circular Sporulation") +
  theme_minimal()

# UV Hypothesis
#Hypothesis two
#H¬0 = Ultraviolet light and temperature do not have a significant effect on the percentage damaged oocysts of kiwi Eimeria spp..
#H1 = Ultraviolet light and temperature do have a significant effect on the percentage damaged oocysts of kiwi Eimeria spp..

df$total_damaged  <-  df$X..damaged.oocysts..all.oocysts. /100 + 0.001

df$total_damaged  <- ifelse(df$total_damaged > 1, 0.99, df$total_damaged ) 

hist(df$total_damaged)


## FINAL DAMAGED Models with the time autocorrelation --------------------
# Model with time structure for day

df <- df[order(df$Day), ]

df$total_damaged_lag1 <- c(NA, head(df$total_damaged, -1))

# Remove the first NA lag so model with time is comparable to NO effect model

dfc <- df[2:nrow(df),]

dfc$total_damaged_lag1

# Interaction model
#3.	Does UV have an effect on sporulation rate and is there a correlation between UV influence and temperature? I think it might be the higher the temperature, the more influence UV has.

model_autocor_UVtemp_interaction_damaged <- glmmTMB(total_damaged ~ Temperature*Light + 
                                              Rep + total_damaged_lag1 + (1 | Day),
                                            data = dfc, family = beta_family(link = "logit"))

# Fit the model with lagged values to capture autocorrelation
model_autocor_damaged <- glmmTMB(total_damaged ~ Temperature + Light + Rep + total_damaged_lag1 + (1 | Day),
                         data = dfc, family = beta_family(link = "logit") )


#
model_autocor_age_damaged <- glmmTMB(total_damaged ~ Temperature + Light + sample_age + 
                                       total_damaged_lag1 + (1 | Day),
                             data = dfc, family = beta_family(link = "logit") )



# No effect model

model_autocor_no_effect_damaged <-  glmmTMB(
  total_damaged ~ 1 + (1 | Day),
  data = dfc,
  family = beta_family(link = "logit")
)


#  Model comparison

summary(model_autocor)

summary(model_autocor_UVtemp_interaction)

AICctab(model_autocor_no_effect_damaged, 
        model_autocor_age_damaged,
        model_autocor_UVtemp_interaction_damaged,
        model_autocor_damaged, logLik=FALSE, 
        weights = TRUE )

summary(model_autocor_no_effect)

# Confint best model

round(confint(model_autocor_age_damaged), 3)


#--- Damaged
head(data)

sdamage <- ggplot(df, aes(x = Day, y = total_damaged, 
                            color = Light)) +
  facet_grid(Rep  ~ Temperature, labeller = label_both) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(aes(group = Light), position = position_dodge(width = 0.3), 
            linetype = "dashed") +
  labs(title = "Interaction Plot of Day, Temperature, and Light on Damaged Oocysts",
       x = "Day",
       y = "Sporulation (%)") +
  scale_color_manual(values = c("Dark" = "purple", "Light" = "orange")) +
  theme_minimal()

setwd('D:/OneDrive - Massey University/kiwi/')

png(filename = 'Total_Damaged.png', res = 400, units = 'cm', width = 20, height = 15      )
sdamage
dev.off()

# Interaction plot damaged and sporulated

plot(df$total_damaged ~ df$Total.sporulated.oocysts, pch=19)

sall <- ggplot(df, aes(x = Total_Sporulation_transformed, y = total_damaged, 
                          color = Light)) +
  facet_grid(Day  ~ Temperature, labeller = label_both) + 
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(aes(group = Light), position = position_dodge(width = 0.3), 
            linetype = "dashed") +
  labs(title = "Correlation between % Damaged and % Sporulated Oocysts",
       x = "Sporulated",
       y = "Damaged (%)") +
  scale_color_manual(values = c("Dark" = "purple", "Light" = "orange")) +
  theme_minimal()

png(filename = 'Damaged_Sporulated_correlation.png', res = 400, units = 'cm', width = 20, height = 15      )
sall
dev.off()

#-----------------------------------