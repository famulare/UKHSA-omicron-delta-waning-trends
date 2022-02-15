# plot VE on logit axis to look at omicron vs delta waning

library(tidyverse)
library(mgcv)
library(matrixStats)
library(mvtnorm)

# load data grabbed from UKHSA 2022 Week 4 vaccine surveillance report 
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf

d <- read.csv('UKHSA-VE-grab.csv',stringsAsFactors = TRUE)

d$weeks <- ordered(d$weeks, levels=c('1','2-4','5-9','10-14','15-19','20-24','25+'))
d$weeks_numeric <- as.numeric(d$weeks) # for modeling. The bins are close enough to uniformly-spaced time intervals for descriptive modeling

d$vax_variant <- interaction(d$vax,d$variant) # for plotting

d$VEnorm <- d$VE/100 # for betareg
d$logitVE <-  qlogis(d$VE/100) # for plotting and log-linear modeling


# plot grabbed data
p1 <- ggplot(d) + 
        geom_point(aes(x=weeks, y=VE, group=vax, color=vax)) +
        geom_line(aes(x=weeks, y=VE, group=vax, color=vax)) +
        facet_wrap('variant') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
p1
ggsave('VE_vax_var.png',units='in',width=6,height=3,device='png')

p2 <- ggplot(d) + 
        geom_point(aes(x=weeks, y=logitVE, group=vax, color=vax)) +
        geom_line(aes(x=weeks, y=logitVE, group=vax, color=vax)) +
        facet_wrap('variant') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
p2
ggsave('VE_logit_vax_var.png',units='in',width=6,height=3,device='png')

# modeling strategy
# 1) compare no-interaction model to model with linear weeks-variant interaction, see interaction has lower AIC and is preffered model
# 2) compare slope-variant interaction to nonlinear weeks-variant, see nonlinear model is preferred
# conclude that it would be great to see this replicated in other settings to mitigate against UKHSA biases, but
# omicron immunity is less stable than delta immunity from wt vaccination. 


# prepping for linear modeling

# plot grouped by variant on one panel
p3 <- ggplot(d) + 
        geom_point(aes(x=weeks, y=VE, group=vax, color=variant)) +
        geom_line(aes(x=weeks, y=VE, group=vax_variant, color=variant)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
p3
p4 <- ggplot(d) + 
        geom_point(aes(x=weeks, y=logitVE, group=vax, color=variant)) +
        geom_line(aes(x=weeks, y=logitVE, group=vax_variant, color=variant)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
p4

# linear in weeks, no interactions
summary(
  m_lin_noint <- gam(VEnorm ~ weeks_numeric + variant + s(vax,bs='re'),
                     family=betar(link="logit"),
                     method='REML',
                     data=d %>% filter(VE>0 & weeks_numeric>1)) # exclude week 1
  )


# linear in weeks, week:variant interaction
summary(
  m_lin_weekvariant <- gam(VEnorm ~ weeks_numeric + variant + weeks_numeric:variant + s(vax,bs='re'),
                     family=betar(link="logit"),
                     method='REML',
                     data=d %>% filter(VE>0 & weeks_numeric>1))  # exclude week 1
)
# weeks_numeric:omicron interaction is negative and significant

# interaction model is preferred (deltaAIC=8)
AIC(m_lin_noint,m_lin_weekvariant)


# plot model fit for modal vaccine (vax effect = 0)
pd <- expand.grid(weeks=levels(d$weeks)[-c(1)], variant=levels(d$variant))
pd$weeks <- ordered(pd$weeks)
pd$weeks_numeric <- as.numeric(pd$weeks)+1
pd$vax <- 'model: variant:weeks interaction'
pd$vax_variant <- interaction(pd$vax,pd$variant)

fit<-predict(m_lin_weekvariant,newdata=pd, se=TRUE, exclude='s(vax)')
pd$logitVE <- fit$fit
pd$logitVE_lower <- fit$fit-2*fit$se.fit
pd$logitVE_upper <- fit$fit+2*fit$se.fit

pd$VE <-plogis(pd$logitVE)*100
pd$VE_lower <-plogis(pd$logitVE_lower)*100
pd$VE_upper <-plogis(pd$logitVE_upper)*100

p4 +
  geom_line(data=pd,aes(x=weeks,y=logitVE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd,aes(x=weeks_numeric,ymin=logitVE_lower,ymax=logitVE_upper,group=variant),alpha=0.2) 
ggsave('VE_logit_var_model.png',units='in',width=4,height=3,device='png')

p3 +
  geom_line(data=pd,aes(x=weeks,y=VE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd,aes(x=weeks_numeric,ymin=VE_lower,ymax=VE_upper,group=variant),alpha=0.2) 
ggsave('VE_var_model.png',units='in',width=4,height=3,device='png')


p2 +
  geom_line(data=pd,aes(x=weeks,y=logitVE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd,aes(x=weeks_numeric,ymin=logitVE_lower,ymax=logitVE_upper,group=variant),alpha=0.2) 
ggsave('VE_logit_vax_var_model.png',units='in',width=6,height=3,device='png')

p1 +
  geom_line(data=pd,aes(x=weeks,y=VE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd,aes(x=weeks_numeric,ymin=VE_lower,ymax=VE_upper,group=variant),alpha=0.2) 
ggsave('VE_vax_var_model.png',units='in',width=6,height=3,device='png')


## isolate time trend difference from level difference
## shift to compare trends starting from same intercept
pd2<-pd

# shift between omicron and delta
shift_intercept <- max(pd2$logitVE[pd2$variant=='omicron']) - max(pd2$logitVE[pd2$variant=='delta'])
  # this is slightly different than coef(m_lin_noint)['variantomicron'] because we're plotting a zeroed out random effect

# zero reference 
maxLogitVE <- max(pd2$logitVE)

# shift by omicron intercept 
pd2$logitVE[pd2$variant=='omicron'] <- pd2$logitVE[pd2$variant=='omicron'] - shift_intercept
pd2$logitVE_lower[pd2$variant=='omicron'] <- pd2$logitVE_lower[pd2$variant=='omicron'] - shift_intercept
pd2$logitVE_upper[pd2$variant=='omicron'] <- pd2$logitVE_upper[pd2$variant=='omicron'] - shift_intercept

#shift to zero reference
pd2$logitVE <- pd2$logitVE - maxLogitVE
pd2$logitVE_lower <- pd2$logitVE_lower - maxLogitVE
pd2$logitVE_upper <- pd2$logitVE_upper - maxLogitVE

# VE scale
# (this doesn't really mean anything, but I was curious to see it anyway)
pd2$VE <-plogis(pd2$logitVE)*100
pd2$VE_lower <-plogis(pd2$logitVE_lower)*100
pd2$VE_upper <-plogis(pd2$logitVE_upper)*100


ggplot(pd2) +
  geom_line(aes(x=weeks,y=logitVE,group=variant,color=variant)) +
  geom_ribbon(aes(x=weeks_numeric-1,ymin=logitVE_lower,ymax=logitVE_upper,group=variant,fill=variant),alpha=0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
ggsave('VE_logit_relative_var_model.png',units='in',width=4,height=3,device='png')


ggplot(pd2) +
  geom_line(aes(x=weeks,y=VE,group=variant,color=variant)) +
  geom_ribbon(aes(x=weeks_numeric-1,ymin=VE_lower,ymax=VE_upper,group=variant,fill=variant),alpha=0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
ggsave('VE_relative_var_model.png',units='in',width=4,height=3,device='png')


### in playing with the start week above, there is some evidence of nonlinearity. Let's gam!

# linear in weeks, week:variant interaction
summary(
  m_nonlin_weekvariant <- gam(VEnorm ~ s(weeks_numeric, by=variant, k=5) + variant + s(vax,bs='re'),
                           family=betar(link="logit"),
                           method='REML',
                           data=d %>% filter(VE>0 & weeks_numeric>1))  # exclude week 1
)
# nonlinear interaction model is preferred (deltaAIC=6.6)
AIC(m_lin_noint,m_lin_weekvariant,m_nonlin_weekvariant)


# plot model fit for modal vaccine (vax effect = 0)
pd3 <- expand.grid(weeks=levels(d$weeks)[-c(1)], variant=levels(d$variant))
pd3$weeks <- ordered(pd3$weeks)
pd3$weeks_numeric <- as.numeric(pd3$weeks)+1
pd3$vax <- 'model: variant:weeks nonlinear interaction'
pd3$vax_variant <- interaction(pd3$vax,pd3$variant)

fit<-predict(m_nonlin_weekvariant,newdata=pd3, se=TRUE, exclude='s(vax)')
pd3$logitVE <- fit$fit
pd3$logitVE_lower <- fit$fit-2*fit$se.fit
pd3$logitVE_upper <- fit$fit+2*fit$se.fit

pd3$VE <-plogis(pd3$logitVE)*100
pd3$VE_lower <-plogis(pd3$logitVE_lower)*100
pd3$VE_upper <-plogis(pd3$logitVE_upper)*100

p4 +
  geom_line(data=pd3,aes(x=weeks,y=logitVE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd3,aes(x=weeks_numeric,ymin=logitVE_lower,ymax=logitVE_upper,group=variant),alpha=0.2) 
ggsave('VE_logit_var_model_nonlinear.png',units='in',width=4,height=3,device='png')

p3 +
  geom_line(data=pd3,aes(x=weeks,y=VE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd3,aes(x=weeks_numeric,ymin=VE_lower,ymax=VE_upper,group=variant),alpha=0.2) 
ggsave('VE_var_model_nonlinear.png',units='in',width=4,height=3,device='png')


p2 +
  geom_line(data=pd3,aes(x=weeks,y=logitVE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd3,aes(x=weeks_numeric,ymin=logitVE_lower,ymax=logitVE_upper,group=variant),alpha=0.2) 
ggsave('VE_logit_vax_var_model_nonlinear.png',units='in',width=6,height=3,device='png')

p1 +
  geom_line(data=pd3,aes(x=weeks,y=VE, group=variant),size=1,linetype='dashed') +
  geom_ribbon(data=pd3,aes(x=weeks_numeric,ymin=VE_lower,ymax=VE_upper,group=variant),alpha=0.2) 
ggsave('VE_vax_var_model_nonlinear.png',units='in',width=6,height=3,device='png')



## plot to compare trends starting from same intercept
pd4<-pd3

# shift between omicron and delta
shift_intercept <- max(pd4$logitVE[pd4$variant=='omicron']) - max(pd4$logitVE[pd4$variant=='delta'])
# this is slightly different than coef(m_nonlin_weekvariant)['variantomicron'] because we're plotting a zeroed out random effect

# zero reference 
maxLogitVE <- max(pd4$logitVE)

# shift by omicron intercept 
pd4$logitVE[pd4$variant=='omicron'] <- pd4$logitVE[pd4$variant=='omicron'] - shift_intercept
pd4$logitVE_lower[pd4$variant=='omicron'] <- pd4$logitVE_lower[pd4$variant=='omicron'] - shift_intercept
pd4$logitVE_upper[pd4$variant=='omicron'] <- pd4$logitVE_upper[pd4$variant=='omicron'] - shift_intercept

#shift to zero reference
pd4$logitVE <- pd4$logitVE - maxLogitVE
pd4$logitVE_lower <- pd4$logitVE_lower - maxLogitVE
pd4$logitVE_upper <- pd4$logitVE_upper - maxLogitVE

# VE scale
pd4$VE <-plogis(pd4$logitVE)*100
pd4$VE_lower <-plogis(pd4$logitVE_lower)*100
pd4$VE_upper <-plogis(pd4$logitVE_upper)*100



# observe that omicron efficacy keeps falling exponentially with time while delta stabilizes
ggplot(pd4) +
  geom_line(aes(x=weeks,y=logitVE,group=variant,color=variant)) +
  geom_ribbon(aes(x=weeks_numeric-1,ymin=logitVE_lower,ymax=logitVE_upper,group=variant,fill=variant),alpha=0.2) +
  theme_bw() 
ggsave('VE_logit_relative_var_model_nonlinear.png',units='in',width=4,height=3,device='png')


ggplot(pd4) +
  geom_line(aes(x=weeks,y=VE,group=variant,color=variant)) +
  geom_ribbon(aes(x=weeks_numeric-1,ymin=VE_lower,ymax=VE_upper,group=variant,fill=variant),alpha=0.2) +
  theme_bw()
ggsave('VE_relative_var_model_nonlinear.png',units='in',width=4,height=3,device='png')


# let's look at (logitVE_omi,1 - logitVE_omi,M)/(logitVE_delta,1 - logitVE_delta,M) 
# as a proxy for the relationship between antibody titers and VE demonstrated in
# Khoury et al https://www.nature.com/articles/s41591-021-01377-8

# to get ratio of differences of smooths, we have to sample from the model

# sampling smooths
# mix of https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
# and https://gist.github.com/noamross/8bf1fc5b2f629b3a7e1eb6b4572e8388
pdat <- expand.grid(weeks = levels(pd4$weeks),#c('15-19','20-24','25+'),
                    variant = c('delta', 'omicron'),
                    vax='modal')
pdat$weeks_numeric <- pd4$weeks_numeric[levels(pd4$weeks)==pdat$weeks]
xp <- predict(m_nonlin_weekvariant, newdata = pdat, type = 'lpmatrix',unconditional = TRUE)

# which cols of xp relate to splines of interest?
c1 <- grepl('delta', colnames(xp))
c2 <- grepl('omicron', colnames(xp))
## which rows of xp relate to variants of interest?
r1 <- with(pdat, variant == 'delta')
r2 <- with(pdat, variant == 'omicron')

X <- xp 

## filter out cols of X related to splines for other vaccines
X <- X[, (c1 | c2)]
## filter out the parametric cols
X[, !grepl('^s\\(', colnames(X))] <- 0

# get mean coefs of the relevant smooths
mean_coef <- coef(m_nonlin_weekvariant)[c1|c2]

# Get the variance-covariance matrix of coefficients, accounting for smoothing
# uncertainty, for the relevant smooths
vcov_mat <- vcov(m_nonlin_weekvariant, unconditional = TRUE)
vcov_mat <- vcov_mat[c1|c2,c1|c2]


# Draw samples from the posterior and make predictions from them
N<-1e3
coefs <- rmvnorm(N, mean = mean_coef, sigma = vcov_mat)
preds <- X %*% t(coefs)

# tidy predictions
preds <- as.data.frame(preds)
colnames(preds)<-paste('rep',1:N,sep='')
preds$variant = c(rep('delta',6),rep('omicron',6))

# diff predictions (logitVE_1 - logitVE_M)
preds[1:6,1:N] <- preds[rep(1, 6),1:N] - preds[1:6,1:N] 
preds[6+(1:6),1:N] <- preds[rep(6+1, 6),1:N] - preds[6+(1:6),1:N]

# ratio preds: omicron over delta
preds <- preds[6+(1:6),1:N]/preds[1:6,1:N]
preds[1,] <- 1 # limit at zero time difference is 1

# index time 
preds$weeks_numeric <- 1:6
preds$weeks <- factor(levels(pd4$weeks),levels=levels(pd4$weeks))

preds$mean = rowMeans(preds[,1:N])
preds$lower = rowQuantiles(as.matrix(preds[,1:N]),probs=0.025)
preds$upper = rowQuantiles(as.matrix(preds[,1:N]),probs=0.975)

titer_drop_ratio <- preds %>% select(weeks, weeks_numeric,mean,lower, upper)

# this shows that at the 25+ week time point, logitVE for omicron falls 1.3 (1.1,1.6) 
# more than delta. Under the Khoury model, this should correspond to a 1.3x lower drop in NAb titers for omicron relative to delta

# here's the evidence from the only paper I've seen with that comparison at 4-6 months post-booster
# compare to zhao fig 1 panels I & J https://www.nejm.org/doi/full/10.1056/NEJMc2119426

zhao_drop = data.frame(median_NAb_1mo = c(2133,516),
                       median_NAb_4to6mo = c(331,51),
                       variant = c('delta','omicron'))

zhao_drop$log_diff <- log(zhao_drop$median_NAb_1mo)-log(zhao_drop$median_NAb_4to6mo)
zhao_drop$waning_drop_ratio <- c(NaN,zhao_drop$log_diff[2]/zhao_drop$log_diff[1])
zhao_drop$weeks <- '20-24' # middle match for 4-6 months (18-26 weeks)
zhao_drop$label <- 'Zhao: omicron-delta'
zhao_drop

# this ratio is solidly inside the predicted confidence interval at 6 months
# and if the best comparison is really the 20-24 week bin, then it's still just inside the upper C1

# Pajon et al https://www.nejm.org/doi/full/10.1056/NEJMc2119912 provides a similar
# comparison for omicron vs D614G.  The observed ratio there should overestimate
# the difference with delta, since delta is antigenically further from WT than D614G.
# WT

# boost taken directly from the paper
pajon_boost_drop = data.frame(median_NAb_1mo = c(2423,850),
                       median_NAb_6mo = c(1067,136),
                       variant = c('D614G','omicron'))

pajon_boost_drop$log_diff <- log(pajon_boost_drop$median_NAb_1mo)-log(pajon_boost_drop$median_NAb_6mo)
pajon_boost_drop$waning_drop_ratio <- c(NaN,pajon_boost_drop$log_diff[2]/pajon_boost_drop$log_diff[1])
pajon_boost_drop$weeks <- '25+' # middle match for 4-6 months (18-26 weeks)
pajon_boost_drop$label <- 'Pajon: omicron-D614G: WA1 boost'
pajon_boost_drop

# this is outside the upper end of the titer_drop_ratio interval, 
# again consistent with the model, although not a tight bound

# 2-dose digitized to account for censoring, since not reported in the paper

library(EnvStats)

pd <-  read.csv('pajon-2dose-grab.csv',stringsAsFactors = TRUE)
pd$weeks[pd$weeks==4] <- '2-4'
pd$weeks[pd$weeks==30] <- '25+'

# checking figure digitization matches paper for delta 2 dose
idx <- pd$weeks=='2-4' & pd$variant=='D614G'
wt1<-10^mean(pd$log10_titer[idx]) # paper reports 1496, so good agreement!
wt1 <- 1496

idx <- pd$weeks=='25+' & pd$variant=='D614G'
wt7<-10^mean(pd$log10_titer[idx]) # paper reports 193, so good agreement!
wt7 <- 193

# censoring-aware estimates for omicron
idx <- pd$weeks=='2-4' & pd$variant=='omicron'
omicron1 <- 10^enormCensored(pd$log10_titer[idx],pd$censored[idx])$parameters[1] # paper reports 43 with only 3 censored variables

idx <- pd$weeks=='25+' & pd$variant=='omicron'
omicron7 <- 10^enormCensored(pd$log10_titer[idx],pd$censored[idx])$parameters[1] # paper reports 23 but 45% of data are censored

pajon_primary_drop = data.frame(median_NAb_1mo = c(wt1,omicron1),
                              median_NAb_6mo = c(wt7,omicron7),
                              variant = c('D614G','omicron'))

pajon_primary_drop$log_diff <- log(pajon_primary_drop$median_NAb_1mo)-log(pajon_primary_drop$median_NAb_6mo)
pajon_primary_drop$waning_drop_ratio <- c(NaN,pajon_primary_drop$log_diff[2]/pajon_primary_drop$log_diff[1])
pajon_primary_drop$weeks <- '25+' # middle match for 4-6 months (18-26 weeks)
pajon_primary_drop$label <- 'Pajon: omicron-D614G: WA1 primary'
pajon_primary_drop

# this shows maturation of omicron NAbs after two doses, relative to WT waning, 
# which is what I expect.  Why then does UKHSA data show the waning pattern it does?


ggplot(titer_drop_ratio) +
  geom_line(aes(x=weeks,y=mean,group='all')) +
  geom_ribbon(aes(x=as.numeric(weeks),ymin=lower,ymax=upper),alpha=0.2) +
  geom_hline(aes(yintercept=1),linetype='dashed')+
  geom_segment(data=zhao_drop,aes(x=5-0.4,xend=6,y=waning_drop_ratio, yend=waning_drop_ratio, color='Zhao: omicron-delta: boost'),color='red') +
  geom_label(data=zhao_drop,aes(1,y=2.1, label=label,color=label),color='red',hjust=0,label.size=NA) +
  geom_point(data=pajon_boost_drop,aes(x=weeks,y=waning_drop_ratio, color='Pajon: omicron-D614G'),color='blue') +
  geom_label(data=pajon_boost_drop,aes(x=1,y=1.96, label='Pajon: omicron-D614G: boost',color='Pajon: omicron-D614G'),color='blue',hjust=0,label.size=NA) +
  geom_point(data=pajon_primary_drop,aes(x=weeks,y=waning_drop_ratio, color='Pajon: omicron-D614G'),color='green') +
  geom_label(data=pajon_primary_drop,aes(x=1,y=1.82, label='Pajon: omicron-D614G: primary',color='Pajon: omicron-D614G'),color='green',hjust=0,label.size=NA) +
  theme_bw() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('ratio of log titer difference') 

ggsave('ratio_of_waning_ratios_model_nonlinear.png',units='in',width=4,height=3,device='png')

############

# taken together, these data are consistent with omicron waning following WT vaccination
# being less stable than delta waning after 3-4 months, both after 2-dose and booster vaccination.