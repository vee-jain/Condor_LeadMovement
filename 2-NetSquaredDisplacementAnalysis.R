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
library(lme4)
library(car)

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
hist(xdata$m.nsd)

#' obs per ID_test combo
plot(table(table(xdata$id_test)))
plot(table(table(xdata$raw_id)))

#### (5) look at the response ----
hist(xdata$m.nsd)
range(xdata$m.nsd)
#' data has a positive skew and zeroes

#' log1p tranformation is possible but seems like some distortions in data
hist(log1p(xdata$m.nsd))

#' transforming by half the minimun non-zero value
xx = signif(0.5*sort(unique(xdata$m.nsd))[2], 2)
hist(log(xdata$m.nsd+xx))

#' this is equivalent to adding a small constant and log transforming
xdata$m.nsd.2 <- xdata$m.nsd + min(xdata$m.nsd[xdata$m.nsd>0])/2
hist(log(xdata$m.nsd.2))

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

#### (6) theoretically identifiable random nsdopes ----
xx.fe.re = fe.re.tab(fe.model = "m.nsd~lead_level*day_prior + age",
                     re = "(1|id_test) +(1|raw_id) ", 
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


xx=xx.fe.re$detailed[4]
xx
xx.fe.re$summary[5]
#' for all but 2 raw_id we have at least 3 unique values of day_prior
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
#' there are 156 individuals in which the levels of factor lead_level
#' is not associated with at least two unique values of day_prior with at least 
#' two observations

xx=xx.fe.re$detailed[9]
xx
xx.fe.re$summary[13]
#' there are more than two unique values of day_prior per level 
#' of lead_level. we need >2 unique values of day_prior in at least two levels of the lead_level
#' include random nsdope

xx=xx.fe.re$detailed[10]
xx
xx.fe.re$summary[14]
#' there are 28 raw_id in which the levels of the factor lead_level does not occur in combination 
#' with at least two unique values of day_prior, each associated with at least two observations

#### (7) setting up model ----
str(xx.fe.re$data)
t.data=xx.fe.re$data

#' centre the dummy variables
t.data$lead_level.2.Exposure=
  t.data$lead_level.2.Exposure-mean(t.data$lead_level.2.Exposure)
t.data$lead_level.3.High.Exposure =
  t.data$lead_level.3.High.Exposure -mean(t.data$lead_level.3.High.Exposure)
t.data$lead_level.4.Clinically.affected=
  t.data$lead_level.4.Clinically.affected-mean(t.data$lead_level.4.Clinically.affected)

#' z-transform the covariates
range(t.data$age)
t.data$z.age=as.vector(scale(t.data$age))
range(t.data$day_prior)
t.data$z.day_prior=as.vector(scale(t.data$day_prior))

#' weight term allows to give more weight to observations that are more 
#' accurate (i.e., based on a larger sample)
plot(t.data$m.nsd, t.data$n_points)

#' response
t.data$m.nsd.2 <- t.data$m.nsd + min(t.data$m.nsd[t.data$m.nsd>0])/2
range(t.data$m.nsd.2)
t.data$log.m.nsd.2 = log(t.data$m.nsd.2)
hist((t.data$log.m.nsd.2))
hist(log(t.data$m.nsd.2))

#### (8) running the model----
full.wac = lme4::lmer(log.m.nsd.2~
                        lead_level*z.day_prior + #testing
                        z.age+ #control
                        (1 + z.day_prior+
                           z.age|id_test)+
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
#' boundary is singular fit, check correlation structure of random effects
#' NaNs produced, model likely too complex

#' remove correlation parameters
full = lme4::lmer(log.m.nsd.2~
                    lead_level*z.day_prior + #testing
                    z.age+ #control
                    (1 + z.day_prior+
                       z.age||id_test)+
                    (1 + lead_level.2.Exposure + 
                       lead_level.3.High.Exposure + 
                       lead_level.4.Clinically.affected +
                       z.day_prior+
                       z.age+
                       lead_level.2.Exposure*z.day_prior + 
                       lead_level.3.High.Exposure*z.day_prior + 
                       lead_level.4.Clinically.affected*z.day_prior||raw_id), 
                  data = t.data,
                  REML = F,
                  weights = n_points,
                  control = lmerControl(optCtrl=list(maxfun= 50000), optimizer = "nloptwrap"))
summary(full)$varcor
#' boundary is singular arises from some random effects being estimated at 0

#' compare log likelihoods of full model with all correlations with that of none 
logLik(full.wac)
logLik(full) 
#'removal of the estimates of all the correlations lead to only a small decrease in max likelihood
#' the boundary (singular) fit message indicates that some random effects are 
#' estimated to be close to their boundary which then implies that the 
#' estimation of their contribution might be un- reliable
#' random intercepts and nsdopes there is a lower bound of zero
#' and correlation parameters are bound between -1 and 1
#' since interest is mainly in fixed effects, leave as is

#### (9) model diagnostics ----
diagnostics.plot(full)
#' not the best fit, some structure in the data
ranef.diagn.plot(full)
#' BLUPs are normally distributed

#' VIF
#note the exclusion of the interaction and polynomial for vif: 
xx=lme4::lmer(log.m.nsd.2~
                lead_level+z.day_prior + #testing
                z.age+ #control
                (1 + z.day_prior+
                   z.age||id_test)+
                (1 + lead_level.2.Exposure + 
                   lead_level.3.High.Exposure + 
                   lead_level.4.Clinically.affected +
                   z.day_prior+
                   z.age+
                   lead_level.2.Exposure*z.day_prior + 
                   lead_level.3.High.Exposure*z.day_prior + 
                   lead_level.4.Clinically.affected*z.day_prior||raw_id), 
              data = t.data,
              REML = F,
              weights = n_points,
              control = lmerControl(optCtrl=list(maxfun= 50000), optimizer = "nloptwrap"))
vif(xx)
vif(xx)[, 3]^2
#' no obvious problem

#### (10) model stability ----
#' for (G)LMMs, stability is usually assessed by excluding levels of grouping 
#' factors rather than individual cases
full.stab=glmm.model.stab(model.res=full, contr=NULL,
                          para=F, data=NULL)

table(full.stab$detailed$lme4.warnings) # no models failed to converge
table(full.stab$detailed$opt.warnings)

is.re=grepl(x=rownames(full.stab$summary), pattern="@")
m.stab.plot(full.stab$summary[!is.re, -1])
m.stab.plot(full.stab$summary[is.re, -1]) 
# wrt random effects, the model is not very stable. but just a control

#### (11) summary ----
summary(full)$varcor
summary(full)$coefficients

#' AIC
summary(full)$AICtab["AIC"]
#' LogLik
logLik(full)
#' no. levels per grouping factor 
summary(full)$ngrps
#' number random effects
length(summary(full)$varcor)
#' total no. of estimated effects
length(fixef(full)) +
  length(summary(full)$varcor) + 1

#' observations per estimated term
length(residuals(full))/
  (length(fixef(full)) + length(summary(full)$varcor) + 1)
#' approx 198 observations per estimated term

#### (12) null model----
#' interested in the effects of lead exposure on mean step length and how that changes with time
null=lme4::lmer(log.m.nsd.2~
                  z.age+ #control
                  (1 + z.day_prior+
                     z.age||id_test)+
                  (1 + lead_level.2.Exposure + 
                     lead_level.3.High.Exposure + 
                     lead_level.4.Clinically.affected +
                     z.day_prior+
                     z.age+
                     lead_level.2.Exposure*z.day_prior + 
                     lead_level.3.High.Exposure*z.day_prior + 
                     lead_level.4.Clinically.affected*z.day_prior||raw_id), 
                data = t.data,
                REML = F,
                weights = n_points,
                control = lmerControl(optCtrl=list(maxfun= 50000), optimizer = "nloptwrap"))
as.data.frame(anova(null, full, test="Chisq"))
#' full null model comparison not significant
#' lead level and day prior do not have an effect
#' no further interpretation

#save.image('nsd.RData')
#load('nsd.RData')
