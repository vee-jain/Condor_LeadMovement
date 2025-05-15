#' Tracking solutions to a persistent threat: 
#' Spatial movement patterns reflect lead exposure in critically endangered 
#' California Condors

#' Script by: Varalika Jain

#' The final processed dataset 'metrics_final' is available on the Movebank data
#' repository associated with this study
#' This code relies on custom functions currently not provided with the code - 
#' please contact me for more information

####----LIBRARIES----####
library(dplyr)
library(ggplot2)
library(kyotil)
library(lme4)
library(car)
library(viridis)
library(ggpubr)

#### (1) Prepare working space----
rm(list = ls())
dev.off()

#### (2) Load data ----
xdata <- read.csv("metrics_final.csv", stringsAsFactors = T)

#### (3) check NA and str ----
any(is.na(xdata))
str(xdata)

xx <- xdata %>% group_by(raw_id, test_date) %>% tally() %>%
  mutate(test_date = as.Date(test_date),
         diff = test_date - lag(test_date)) #all tests at least 30d apart

#### (4) look at predictors ----
#' indv > tests > days

#' multiple obs per individual
plot(table(table(xdata$raw_id)))

#' multiple testing dates per individual
xx= xdata %>% group_by(raw_id, test_date) %>% tally() %>%
  group_by(raw_id) %>% tally()
hist(xx$n)
mean(xx$n>1)

#' range of testing days per individual and test
range(table(xdata$id_test))

#' lead exposure levels in individuals (30 days of data max)
table(xdata$id_test, xdata$lead_level)

hist(xdata$age) #in days
plot(xdata$lead_level)
hist(xdata$day_prior)

#' obs per ID_test combo
plot(table(table(xdata$id_test)))
plot(table(table(xdata$raw_id)))

#### (5) look at the response ----
hist(xdata$mcp)
range(xdata$mcp)
#' data has a positive skew and zeroes

#' log1p tranformation is possible but seems like some distortions in data
hist(log1p(xdata$mcp))

#' transforming by half the minimun non-zero value
xx = signif(0.5*sort(unique(xdata$mcp))[2], 2)
hist(log(xdata$mcp+xx))

#' this is equivalent to adding a small constant and log transforming
xdata$mcp.2 <- xdata$mcp + min(xdata$mcp[xdata$mcp>0])/2
hist(log(xdata$mcp.2))

#### (6) grouping factors----
#' observations per subject and testing date
xx=aggregate(x=1:nrow(xdata), by=xdata[, c("raw_id", "test_date")],
             FUN=length)
sum(xx$x > 1)
mean(xx$x > 1)
#' for 154 (99%) combinations of indv and test date, there are multiple samples
#' include random effect 

xx=aggregate(x=1:nrow(xdata), by=xdata[, c("raw_id", "date")],
             FUN=length)
sum(xx$x > 1)
mean(xx$x > 1) 
#' for 0 combinations of indv and date, there are multiple samples
#' no random effect 

#### (6) theoretically identifiable random mcpopes ----
xx.fe.re = fe.re.tab(fe.model = "mcp~lead_level*day_prior + age",
                     re = "(1|id_test)+(1|raw_id)", 
                     other.vars="n_points",
                     data = xdata)

xx=xx.fe.re$detailed[1]
xx
xx.fe.re$summary[1] 
#' in 154 id_tests, we have only one level of the factor lead_level with at least two observations
#' don't include random nsdope

xx=xx.fe.re$detailed[2]
xx
xx.fe.re$summary[2]
#' in 21 raw_id, we have only one level of the factor lead_level with at least two observations
#' 31 raw_id have 1 or more unique levels of lead level with at least 2 obs
#' include random nsdope

xx=xx.fe.re$detailed[3]
xx
xx.fe.re$summary[3]
#' for all but 3 individuals we have at least 3 unique values of day_prior in id_test
#' include random nsdope
xx.fe.re$summary[4]
#' each id_test has only one day prior with less than two levels

xx=xx.fe.re$detailed[4]
xx
xx.fe.re$summary[5]
#' for all but 1 raw_id we have at least 3 unique values of day_prior
#' include random nsdope
xx.fe.re$summary[6]
#' for 16 raw_id we have no single unique value of day prior with at least two observation
#' in 36 raw_id, we have 5 or more single unique values of day prior with at least two observations

xx=xx.fe.re$detailed[5]
xx
xx.fe.re$summary[7]
#' for all but 3 id_tests we have at least 3 unique values of age
#' include random nsdope 
xx.fe.re$summary[8]

xx=xx.fe.re$detailed[6]
xx
xx.fe.re$summary[9]
#' for all but 1 raw_id we have at least 3 unique values of age
#' include age as random nsdope in raw_id
xx.fe.re$summary[10]

xx=xx.fe.re$detailed[7]
xx
xx.fe.re$summary[11]
#' in the body of the table this shows whether (>2) or not (!>2) 
#' there are more than two unique values of day_prior per level 
#' of lead_level. We need >2 unique values of the covariate in at least two levels 
#' of the fixed effects factor
#' we can’t include the random nsdope of the interaction between 
#' lead_level and day_prior (since lead_level doesn’t vary within id_test)
#' no random nsdope

xx=xx.fe.re$detailed[8]
xx
xx.fe.re$summary[12]

xx=xx.fe.re$detailed[9]
xx
xx.fe.re$summary[13]
#' there are more than two unique values of day_prior per level 
#' of lead_level. we need >2 unique values of day_prior in at least two levels of the lead_level
#' include random nsdope
xx.fe.re$summary[14]
#' there are 28 raw_id in which the levels of the factor lead_level does not occur in combination 
#' with at least two unique values of day_prior, each associated with at least two observations

#### (7) setting up model ----
str(xx.fe.re$data)
t.data=xx.fe.re$data

#' center the dummy variables
t.data$lead_level.2.Exposure=t.data$lead_level.2.Exposure-mean(t.data$lead_level.2.Exposure)
t.data$lead_level.3.High.Exposure =t.data$lead_level.3.High.Exposure -mean(t.data$lead_level.3.High.Exposure)
t.data$lead_level.4.Clinically.affected=t.data$lead_level.4.Clinically.affected-mean(t.data$lead_level.4.Clinically.affected)

#' z-transform the covariates
range(t.data$age)
t.data$z.age=as.vector(scale(t.data$age))
mean(t.data$age)
sd(t.data$age)
range(t.data$day_prior)
t.data$z.day_prior=as.vector(scale(t.data$day_prior))
mean(t.data$day_prior)
sd(t.data$day_prior)

#' weight term allows to give more weight to observations that are more 
#' accurate (i.e., based on a larger sample)
plot(t.data$mcp, t.data$n_points)

#' response
t.data$mcp.2 <- t.data$mcp
t.data$log1p.mcp.2 = log1p(t.data$mcp.2)
hist(t.data$log1p.mcp.2) #not the best distribution

#### (8) running the model----
full.wac = lme4::lmer(log1p.mcp.2~
                        lead_level*z.day_prior + 
                        z.age+
                        (1 + z.day_prior +
                           z.age | id_test)+
                        (1 + lead_level.2.Exposure + 
                           lead_level.3.High.Exposure + 
                           lead_level.4.Clinically.affected +
                           z.day_prior+
                           z.age+
                           lead_level.2.Exposure*z.day_prior + 
                           lead_level.3.High.Exposure*z.day_prior + 
                           lead_level.4.Clinically.affected*z.day_prior|raw_id), 
                      data = t.data,
                      REML = F,
                      weights = n_points,
                      control = lmerControl(optCtrl=list(maxfun= 50000), optimizer = "nloptwrap"))
summary(full.wac)$varcor
#' boundary is singular fit, but correlation structure looks okay
#' proceeding with model with all correlations

#### (9) model diagnostics ----
diagnostics.plot(full.wac)
#' not the best fit, some structure in the data
ranef.diagn.plot(full.wac)
#' BLUPs are normally distributed

#' VIF
#note the exclusion of the interaction for vif: 
xx=lme4::lmer(log1p.mcp.2~
                lead_level+z.day_prior + 
                z.age+
                (1 + z.day_prior +
                   z.age | id_test)+
                (1 + lead_level.2.Exposure + 
                   lead_level.3.High.Exposure + 
                   lead_level.4.Clinically.affected +
                   z.day_prior+
                   z.age+
                   lead_level.2.Exposure*z.day_prior + 
                   lead_level.3.High.Exposure*z.day_prior + 
                   lead_level.4.Clinically.affected*z.day_prior|raw_id), 
              data = t.data,
              REML = F,
              weights = n_points,
              control = lmerControl(optCtrl=list(maxfun= 10000), optimizer = "nloptwrap"))
vif(xx)
vif(xx)[, 3]^2
#' no obvious problem

#### (10) model stability ----
#' for (G)LMMs, stability is usually assessed by excluding levels of grouping 
#' factors rather than individual cases
full.stab=glmm.model.stab(model.res=full.wac, contr=NULL,
                          para=F, data=NULL)

table(full.stab$detailed$lme4.warnings) # 2 models failed to converge
table(full.stab$detailed$opt.warnings)

is.re=grepl(x=rownames(full.stab$summary), pattern="@")
m.stab.plot(full.stab$summary[!is.re, -1])
round(full.stab$summary[!is.re, -1], 3)
m.stab.plot(full.stab$summary[is.re, -1]) 
# wrt random effects, the model is not very stable. but just a control

#### (11) summary ----
summary(full.wac)$varcor
summary(full.wac)$coefficients
round(summary(full.wac)$coefficients, 3)

#' AIC
summary(full.wac)$AICtab["AIC"]
#' LogLik
logLik(full.wac)
#' no. levels per grouping factor 
summary(full.wac)$ngrps
#' number random effects
length(summary(full.wac)$varcor)
#' total no. of estimated effects
length(fixef(full.wac)) +
  length(summary(full.wac)$varcor) + 1

#' observations per estimated term
length(residuals(full.wac))/
  (length(fixef(full.wac)) + length(summary(full.wac)$varcor) + 1)
#' approx 363.5 observations per estimated term

#### (12) null model----
#' interested in the effects of lead exposure on mean step length and how that changes with time
null=lme4::lmer(log1p.mcp.2~
                  z.age+
                  (1 + z.day_prior +
                     z.age | id_test)+
                  (1 + lead_level.2.Exposure + 
                     lead_level.3.High.Exposure + 
                     lead_level.4.Clinically.affected +
                     z.day_prior+
                     z.age+
                     lead_level.2.Exposure*z.day_prior + 
                     lead_level.3.High.Exposure*z.day_prior + 
                     lead_level.4.Clinically.affected*z.day_prior|raw_id), 
                data = t.data,
                REML = F,
                weights = n_points,
                control = lmerControl(optCtrl=list(maxfun= 10000), optimizer = "nloptwrap"))
as.data.frame(anova(null, full.wac, test="Chisq"))

#' full null model comparison significant
#' lead level and/or day prior and/or the interaction do have an effect

#### (13) fixed effects inference ----
full.reml=lmerTest::lmer(log1p.mcp.2~
                           lead_level*z.day_prior + 
                           z.age+
                           (1 + z.day_prior +
                              z.age | id_test)+
                           (1 + lead_level.2.Exposure + 
                              lead_level.3.High.Exposure + 
                              lead_level.4.Clinically.affected +
                              z.day_prior+
                              z.age+
                              lead_level.2.Exposure*z.day_prior + 
                              lead_level.3.High.Exposure*z.day_prior + 
                              lead_level.4.Clinically.affected*z.day_prior|raw_id), 
                         data = t.data,
                         REML = T,
                         weights = n_points,
                         control = lmerControl(optCtrl=list(maxfun= 10000), optimizer = "nloptwrap"))

summary(full.reml)
round(summary(full.reml)$coefficients, 3)

#' getting more realiable p-values from drop1
tests=as.data.frame(drop1(full.wac, test="Chisq"))
round(tests, 3)
#' the interaction is not significant

#' removing the interaction term 
summary(xx)
round(summary(xx)$coefficients, 3)

tests.wni=as.data.frame(drop1(xx, test="Chisq"))
round(tests.wni, 3)

emmeans::emmeans(xx, pairwise ~ lead_level)

xx.reml=lmerTest::lmer(log1p.mcp.2~
                         lead_level+z.day_prior + 
                         z.age+
                         (1 + z.day_prior +
                            z.age | id_test)+
                         (1 + lead_level.2.Exposure + 
                            lead_level.3.High.Exposure + 
                            lead_level.4.Clinically.affected +
                            z.day_prior+
                            z.age+
                            lead_level.2.Exposure*z.day_prior + 
                            lead_level.3.High.Exposure*z.day_prior + 
                            lead_level.4.Clinically.affected*z.day_prior|raw_id), 
                       data = t.data,
                       REML = T,
                       weights = n_points,
                       control = lmerControl(optCtrl=list(maxfun= 10000), optimizer = "nloptwrap"))

#### (14) bootstrapping ----
#boot.full=boot.lmer(m=full.wac, discard.warnings=F,
#                  nboots=1000, para=T, n.cores=5, resol=1000, level=0.95)

m.stab.plot(boot.full$ci.estimates)
round(boot.full$ci.estimates, 3)

#### (15) effect sizes----
#' random effects
MuMIn::r.squaredGLMM(object=full.wac)

####---- (16) plotting ----
# Calculate typical values for other predictors
typical_values <- data.frame(
  lead_level = factor(c("1-Background", "2-Exposure", "3-High Exposure", "4-Clinically affected")),
  z.day_prior = mean(t.data$z.day_prior, na.rm = TRUE),
  z.age = mean(t.data$z.age, na.rm = TRUE),
  weights = mean(t.data$n_points, na.rm = TRUE)
)

# Expand typical values to match the length of z.doy_seq
new_data <- typical_values
preds.1 <- predict(xx, newdata = new_data, se.fit = TRUE,
                   re.form = NA, type = "link")

# Calculate confidence intervals on the response scale
conf_intervals <- data.frame(
  fit = (preds.1$fit),
  lwr = (preds.1$fit - 1.96 * preds.1$se.fit),  # 95% CI lower bound
  upr = (preds.1$fit + 1.96 * preds.1$se.fit)   # 95% CI upper bound
)

# Plot 
preds.final <- cbind(new_data, ci) #use bootstrapped confidence estimates for final plot 

break_seq = seq(11, 19, 0.5)

my_comparisons <- list( c("1-Background", "4-Clinically affected"), 
                        c("2-Exposure", "4-Clinically affected"))

quartz(height = 5, width = 5.5)
ggplot(preds.final, aes(x = lead_level, y = (orig), ymin = (X2.5.), ymax = (X97.5.))) +
  geom_errorbar(width = 0.2) +
  geom_point(stat = "identity", color = "black", size = 3) +
  theme_minimal() +
  labs(x = "Lead level",
       y = expression("MCP in km"^2*"")) +
  scale_y_continuous(breaks = break_seq, 
                     labels = sprintf("%.2f", round(expm1(break_seq)/1000000, digits = 2)),
                     limits = c(11,19))+
  geom_segment(aes(x = "1-Background", xend = "4-Clinically affected", y = 18.245, yend = 18.245), 
               linetype = "solid", lwd = 0.2)+
  geom_segment(aes(x = "2-Exposure", xend = "4-Clinically affected", y = 18.6, yend = 18.6), 
               linetype = "solid", lwd = 0.2) +
  annotate("text", x = 2.5, y = 18.445, label = "*", size = 3) +
  annotate("text", x = 3, y = 18.82, label = "*", size = 3) 

dev.copy2pdf(file="./plot_new_mcp_revision.pdf")

#### (17) final model reportings----
#' estimate and standard error
round(summary(full.wac)$coefficients, 3)
#' confidence intervals
round(boot.full$ci.estimates, 3)
#' chisq effect sizes 
round(tests, 3)
#' min max vals from stability
round(full.stab$summary[!is.re, -1], 3)
#' random effects
summary(full.wac)$varcor

#' removing the interaction term 
#' estimate and standard error
round(summary(xx)$coefficients, 3)
#' chisq effect sizes 
round(tests.wni, 3)

#' pairwise comparisons
emmeans::emmeans(xx, pairwise ~ lead_level)$contrasts

#save.image('mcp.RData')
#load('mcp.RData')
