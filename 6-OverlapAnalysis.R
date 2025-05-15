#' Tracking solutions to a persistent threat: 
#' Spatial movement patterns reflect lead exposure in critically endangered 
#' California Condors

#' Script by: Varalika Jain

#' The final processed dataset 'overlap_final' is available on the Movebank data
#' repository associated with this study
#' This code relies on custom functions currently not provided with the code - 
#' please contact me for more information

####----LIBRARIES----####
library(dplyr)

#### (1) Prepare working space----
rm(list = ls())
dev.off()

#### (2) Load data ----
xdata <- read.csv("overlap_final.csv", stringsAsFactors = T)

#### (3) check NA and str ----
any(is.na(xdata))
str(xdata)

#### (4) look at predictors ----
#' dyads
xdata$to_raw_id <- as.numeric(as.factor(xdata$to_raw_id))
xdata$from_raw_id <- as.numeric(as.factor(xdata$from_raw_id))

xdata$to_raw_id <- sprintf("%03d", xdata$to_raw_id)
xdata$from_raw_id <- sprintf("%03d", xdata$from_raw_id)

xdata$dyad <- apply(xdata[, c("from_raw_id", "to_raw_id")], 1, function(x) paste0(sort(x), collapse = ""))

plot(table(table(xdata$dyad)))
plot(table(table(xdata$lead_dyad)))

plot(xdata$lead_dyad)

#### (5) look at the response ----
hist(xdata$overlap_95)

#' beta distribution, transform data 
xdata$tr.overlap_95=
  (xdata$overlap_95*(nrow(xdata) - 1) + 0.5)/nrow(xdata)

#' in decently sized data sets, the transformation has hardly any effect
plot(xdata$overlap_95, xdata$tr.overlap_95)
abline(a=0, b=1)

#### (6) grouping factors----
#' observations per subject and testing date
plot(table(xdata$dyad))
#' 2 obs per dyad

xx = table(xdata$dyad)
table(apply(X=xx>2, MARGIN=1, FUN=sum))

#### (6) theoretically identifiable random slopes ----
xx.fe.re = fe.re.tab(fe.model = "overlap_95~lead_dyad",
                     re = "(1|dyad)",
                     data = xdata)

xx=xx.fe.re$detailed[1]
xx
xx.fe.re$summary[1] 
#' random slope pf lead_dyad in dyad

#### (7) setting up model ----
str(xx.fe.re$data)
t.data=xx.fe.re$data
xx=xx.fe.re$data

center_dummy <- function(x) {
  mean_value <- mean(x)  # Calculate the mean of the column
  difference <- x - mean_value  # Calculate the difference of each element from the mean
  return(difference)  # Return the differences
}

t.data[,4:18] <- apply(t.data[,4:18], MARGIN=2, FUN=center_dummy)

#' transform response in this dataset
t.data$tr.overlap_95=
  (t.data$overlap_95*(nrow(t.data) - 1) + 0.5)/nrow(t.data)
#' in decently sized data sets, the transformation has hardly any effect
plot(t.data$overlap_95, t.data$tr.overlap_95)
abline(a=0, b=1)

#### (8) running the model----
library(glmmTMB)
#full.wac=glmmTMB(tr.overlap_95 ~ lead_dyad+
#                  (1 + lead_dyad.1.Background.2.Exposure+
#                    lead_dyad.1.Background.3.High.Exposure+
#                   lead_dyad.1.Background.4.Clinically.affected+
#                  lead_dyad.2.Exposure.1.Background +
#                 lead_dyad.2.Exposure.2.Exposure+
#                lead_dyad.2.Exposure.3.High.Exposure+
#               lead_dyad.2.Exposure.4.Clinically.affected+
#              lead_dyad.3.High.Exposure.1.Background+
#             lead_dyad.3.High.Exposure.2.Exposure+
#            lead_dyad.3.High.Exposure.3.High.Exposure +
#           lead_dyad.3.High.Exposure.4.Clinically.affected+
#          lead_dyad.4.Clinically.affected.1.Background+
#         lead_dyad.4.Clinically.affected.2.Exposure +
#        lead_dyad.4.Clinically.affected.3.High.Exposure+
#       lead_dyad.4.Clinically.affected.4.Clinically.affected|dyad),
# family=beta_family, data=t.data)
#' the full model fails to converge
summary(full.wac)$varcor

full = glmmTMB(tr.overlap_95 ~ lead_dyad+
                 (1 + lead_dyad.1.Background.2.Exposure+
                    lead_dyad.1.Background.3.High.Exposure+
                    lead_dyad.1.Background.4.Clinically.affected+
                    lead_dyad.2.Exposure.1.Background +
                    lead_dyad.2.Exposure.2.Exposure+
                    lead_dyad.2.Exposure.3.High.Exposure+
                    lead_dyad.2.Exposure.4.Clinically.affected+
                    lead_dyad.3.High.Exposure.1.Background+
                    lead_dyad.3.High.Exposure.2.Exposure+
                    lead_dyad.3.High.Exposure.3.High.Exposure +
                    lead_dyad.3.High.Exposure.4.Clinically.affected+
                    lead_dyad.4.Clinically.affected.1.Background+
                    lead_dyad.4.Clinically.affected.2.Exposure +
                    lead_dyad.4.Clinically.affected.3.High.Exposure+
                    lead_dyad.4.Clinically.affected.4.Clinically.affected||dyad),
               family=beta_family, data=t.data)

logLik(full) 

#### (9) model diagnostics ----
overdisp.test(full)

#no overdispersion
ranef.diagn.plot(full)
#' BLUPs are normally distributed

#' VIF
#' only one term so not needed

#### (10) model stability ----
#' for (G)LMMs, stability is usually assessed by excluding levels of grouping 
#' factors rather than individual casesin
stab.full=glmmTMB.stab(model.res=full, para=T, n.cores=5,
                       data=t.data)

table(stab.full$detailed$converged)
m.stab.plot(stab.full$summary[, -1])
round(stab.full$summary[1:16, -1], 3)
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
#' total no. of estimated e
#' \ffects
length(fixef(full)) +
  length(summary(full)$varcor) + 1

#' observations per estimated term
length(residuals(full))/
  (length(fixef(full)) + length(summary(full)$varcor) + 1)
#'  5430 observations per estimated term

#### (12) null model----
null=glmmTMB(tr.overlap_95 ~ 1 +
               (1 + lead_dyad.1.Background.2.Exposure+
                  lead_dyad.1.Background.3.High.Exposure+
                  lead_dyad.1.Background.4.Clinically.affected+
                  lead_dyad.2.Exposure.1.Background +
                  lead_dyad.2.Exposure.2.Exposure+
                  lead_dyad.2.Exposure.3.High.Exposure+
                  lead_dyad.2.Exposure.4.Clinically.affected+
                  lead_dyad.3.High.Exposure.1.Background+
                  lead_dyad.3.High.Exposure.2.Exposure+
                  lead_dyad.3.High.Exposure.3.High.Exposure +
                  lead_dyad.3.High.Exposure.4.Clinically.affected+
                  lead_dyad.4.Clinically.affected.1.Background+
                  lead_dyad.4.Clinically.affected.2.Exposure +
                  lead_dyad.4.Clinically.affected.3.High.Exposure+
                  lead_dyad.4.Clinically.affected.4.Clinically.affected||dyad),
             family=beta_family, data=t.data)
as.data.frame(anova(null, full, test="Chisq"))
#' lead level dyad has an effect

#### (13) fixed effects inference ----
round(summary(full)$coefficients$cond,3)
coefs=fixef(full)$cond
coefs
xx=c(coefs["(Intercept)"],
     coefs["(Intercept)"]+coefs["lead_dyad1-Background-2-Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad1-Background-3-High Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad1-Background-4-Clinically affected"],
     coefs["(Intercept)"]+coefs["lead_dyad2-Exposure-1-Background"],
     coefs["(Intercept)"]+coefs["lead_dyad2-Exposure-2-Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad2-Exposure-3-High Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad2-Exposure-4-Clinically affected"],
     coefs["(Intercept)"]+coefs["lead_dyad3-High Exposure-1-Background"],
     coefs["(Intercept)"]+coefs["lead_dyad3-High Exposure-2-Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad3-High Exposure-3-High Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad3-High Exposure-4-Clinically affected"],
     coefs["(Intercept)"]+coefs["lead_dyad4-Clinically affected-1-Background"],
     coefs["(Intercept)"]+coefs["lead_dyad4-Clinically affected-2-Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad4-Clinically affected-3-High Exposure"],
     coefs["(Intercept)"]+coefs["lead_dyad4-Clinically affected-4-Clinically affected"])
round(exp(xx)/(1+exp(xx)),3)

#' background, exposure, high exposed & cr. affected birds overlap most with cr. affected birds 
#' exposure, high exposure ainnd critically affected birds overlap least with background

#' pairwise comparisons
emmeans::emmeans(full, pairwise ~ lead_dyad)$contrasts

#### (14) bootstrapping ----
#' bootstrapping (truncation formula)
beta.tr<-function(x){
  if(any(is.na(x))){warning("x comprises NAs")}
  return((x*(length(x) - 1) + 0.5)/length(x))
}

#' simulate 100 new responses
sim.resp = data.frame(matrix(NA, nrow = 32580, ncol = 100))

# Run the simulation function 100 times
for (i in 1:100) {
  # Simulate data and store it in a new column of results_df
  sim.resp[, i] <- simulate(object = full, seed = 1)[, 1]
}

#' what happens when truncation is applied
xy = as.data.frame(apply(sim.resp, 2, beta.tr))
plot(xy[,1], sim.resp[,1])

#' Applying the truncation as responses are not 0<y<1
range(sim.resp)
sim.resp <- as.data.frame(apply(sim.resp, 2, beta.tr))
range(sim.resp)

refit_with_full <- function(column_data) {
  glmmTMB::refit(full, column_data) #apply the refit for model 'full' on each column
}

sim.models <- lapply(sim.resp, refit_with_full)

summary(sim.models[[1]])$coefficients$cond[,1]

##extract results:
all.res=sim.models[!sapply(sim.models, is.null)]
##estimates, fixef effects:
all.fe=lapply(all.res, function(x){summary(x)$coefficients$cond[,1]})

# Convert the list into a dataframe
all.fe.df <- do.call(rbind, all.fe)

# Convert row names to a column
all.fe.df <- data.frame(row.names = rownames(all.fe.df), all.fe.df, check.names = FALSE)

# Reset row names
rownames(all.fe.df) <- NULL

#confidence limits
ci=apply(all.fe.df, 2, quantile, prob=c((1-0.95)/2, 1-(1-0.95)/2), na.rm=T)
ci = as.data.frame(ci)

##fixed effects:
orig=as.data.frame(summary(full)$coefficients$cond[,1])
orig$names <- rownames(orig)

xy=ci[, match(orig$names, colnames(ci))]

ci.fe=data.frame(orig=orig$`summary(full)$coefficients$cond[, 1]`, t(ci))

ci.fe

#### (15) effect sizes----
#' random effects
MuMIn::r.squaredGLMM(object=full)

####---- (16) plotting ----
x1 = levels(t.data$lead_dyad)
typical_values <- data.frame(
  lead_dyad = x1)
preds.1 <- predict(full, newdata = typical_values, se.fit = TRUE,
                   re.form = NA, type = "link")
ci.1 <- ci.fe
ci.1$lead_dyad <- rownames(ci.1)

ci.1$dyad <-  strtrim(ci.1$lead_dyad,10)
ci.1$dyad[1] <- ci.1$dyad[2]

labels <- as.data.frame(strsplit(as.character(ci.1$lead_dyad),'-'))

ci.3 <- ci.1 %>%
  tidyr::separate(lead_dyad, c("lead", "from", "num", "to"), "-")

ci.3$from[1] <- ci.3$from[2]
ci.3$to[1] <- ci.3$from[1]

names <- rownames(ci.3)
ci.3$lead_dyad <- names

ci.3$from <- factor(ci.3$from, 
                    levels = c("Background", "Exposure", "High Exposure", "Clinically affected"))

ci.3$to <- factor(ci.3$to, 
                  levels = c("Background", "Exposure", "High Exposure", "Clinically affected"))


break_seq <- seq(-0.5,0.6,0.15)

quartz(height = 6, width = 7)
ggplot(ci.3, aes(x = to, y = (orig), ymin = (`X2.5.`), ymax = (`X97.5.`))) +
  geom_errorbar(width = 0.2) +
  geom_point(stat = "identity", color = "black", size = 3) +
  theme_minimal() +
  labs(x = "Lead level",
       y = expression("Proportion overlap")) +
  scale_y_continuous(breaks = break_seq, 
                     labels = sprintf("%.2f",round(exp(break_seq)/(1+exp(break_seq)), 2)),
                     limits = c(-0.5,0.63))+
  facet_grid(~from, scales = "free_x")+ 
  theme(axis.text.x = element_text(angle=90))+
  geom_segment(aes(x = "Background", xend = "Exposure", y = 0.50, yend = 0.50), linetype = "solid", lwd = 0.2)+
  geom_segment(aes(x = "Background", xend = "High Exposure", y = 0.55, yend = 0.55), linetype = "solid", lwd = 0.2) +
  geom_segment(aes(x = "Background", xend = "Clinically affected", y = 0.60, yend = 0.60), linetype = "solid", lwd = 0.2) +
  annotate("text", x = 1.5, y = 0.515, label = "*") +
  annotate("text", x = 2, y = 0.565, label = "*") +
  annotate("text", x = 2.5, y = 0.615, label = "*")

dev.copy2pdf(file="./overlap_new_plot.pdf")

save.image('overlap_analysis.RData')
load('overlap_analysis.RData')




