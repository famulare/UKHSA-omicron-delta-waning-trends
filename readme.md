# Comparison of vaccine efficacy waning trends over time for omicron vs delta using UHKSA estimates
Mike Famulare

[Institute for Disease Modeling](www.idmod.org) | Global Health | Bill \& Melinda Gates Foundation

4 February 2022

## Summary

Here, I take a look at how the trend of waning vaccine efficacy over time in data from UKHSA differs between omicron and delta, and show that it is consistent with limited data on neutralizing antibody titer data. Taken together, I think we are seeing evidence that serologic immunity to omicron from wild-type (WT, late-2019-era SARS-CoV-2) vaccination is less stable than against delta (and likely other variants less distantly related to WT than omicron).

With estimates of vaccine efficacy against symptomatic disease grabbed from [UKHSA COVID-19 vaccine surveillance report -- Week 4 -- 27 Jan 2022](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf), I used mixed-effects beta-regression modeling with logit link to isolate the shared trend over time from the differences in baseline protection level across various vaccine combinations. The isolated temporal trends indicate that waning of protection against symptomatic disease occurs at similar rates for omicron and delta during the first 19 weeks post-vaccination. In weeks 20-24 and 25+ after vaccination, the trends for delta and omicron diverge, with delta efficacy remaining more stable as omicron efficacy continues to fall exponentially toward zero. While this difference is only supported by VE estimates for people with two-dose vaccination because post-booster data at that time interval is not yet available, this data suggests that protection from symptomatic disease following WT vaccination is less stable for omicron than for delta and likely other variants more antigenically similar to WT.  

The observation that efficacy against symptomatic disease following WT vaccination is less stable for omicron is consistent with evidence from neutralization studies in boosted populations. For example, [Pajon et al 2022](https://t.co/lBxMqBjWpL) and [Zhao et al 2022](https://t.co/mswdI6Lqxk)) show that neutralizing antibody titers following WT vaccination wane further for omicron four to six months after booster dose than do titers against variants that are more antigenically similar to WT. I use the logit-linear relationship between vaccine efficacy and neutralizing antibody titers first demonstrated in [Khoury et al 2021](https://www.nature.com/articles/s41591-021-01377-8) to show that the observed difference in UKHSA vaccine efficacy waning is consistent with the observed differences in both neutralization studies. 

Taken together, it is plausible that repeated WT vaccination will not be capable of producing stable protection against symptomatic disease for omicron(-like) variants more than 3-4 months post-vaccination.  Thus, omicron-specific (or omicron-descendant) vaccines may be necessary to stimulate more durable immunity, even if WT vaccine offers similar protection shortly after vaccination. 

## Analysis

Everything necessary to reproduce this analysis is available at [Github/famulare/UKHSA-omicron-delta-waning-trends](https://github.com/famulare/UKHSA-omicron-delta-waning-trends).

Estimates of vaccine efficacy (VE) against symptomatic disease were grabbed from [UKHSA COVID-19 vaccine surveillance report -- Week 4 -- 27 Jan 2022](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf) using [Automeris Web Plot Digitizer](https://automeris.io/WebPlotDigitizer/) and is available at [UKHSA-VE-grab.csv](UKHSA-VE-grab.csv).

To model the waning trends for omicron and delta efficacy as a function of weeks since last vaccination, starting from weeks 2-4, I used mixed-effects beta regression with logit link using the [mcgv](https://rdrr.io/cran/mgcv/) package in R. Beta regression is appropriate to model outcomes on a (0,1) scale. The outcome was vaccine efficacy and the covariates (weeks since last vaccination), variant (omicron or delta), and vaccine schedule (mix of AZ, Pfizer, and Moderna in 2 and 3 dose schedules). The differences in protection level among the various vaccines were modeled as random intercepts. After performing model selection among linear/no-interaction, linear/weeks-variant interaction, and non-linear/weeks-variant interaction models, I found that the best fitting model was the nonlinear interaction model. (See [UKHSA-VE-grab.R](UKHSA-VE-grab.R) for full model specification and workflow.) The model results shown in the figures below represent the prediction and 95\% empirical Bayesian credible interval for efficacy vs weeks since vaccition for the "modal vaccine" (an archetypal vaccine with the random intercept equal to zero). 

![VE_vax_var.png](VE_vax_var_model_nonlinear.png)
**Figure 1.** Vaccine efficacy estimates for various vaccine combos and omicron and delta variants. Colored dots and lines are directly captured from [UKHSA COVID-19 vaccine surveillance report -- Week 4 -- 27 Jan 2022](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf) using [Automeris Web Plot Digitizer](https://automeris.io/WebPlotDigitizer/). Dashed black line and gray ribbon shows best fit and 95\% empirical bayesian credible interval for the modal vaccine under the model.

![VE_vax_var.png](VE_logit_vax_var_model_nonlinear.png)
**Figure 2.** Logit of vaccine efficacy estimates shown in Figure 1. The logit transform linearizes vaccine efficacy estimates and is also proportional to the mechanistic relationship between vaccine efficacy and log neutralizing antibody titers, as shown by [Khoury et al 2021](https://www.nature.com/articles/s41591-021-01377-8) and others.

Using the model, we can isolate the waning trend from the difference in intercept between the variants. Below, we see that the waning trends are very similar through 19 weeks, but diverge thereafter. VE for omicron continues to fall linearly in logit space while the waning rate for delta slows down. If this pattern reflects biology and isn't just an artifact of some uncorrected confounding, then we are seeing evidence that immunity against omicron (and omicron-like variants) from the WT-based vaccines currently in use will be less stable than for delta. 

![VE_logit_relative_var_model_nonlinear.png](VE_logit_relative_var_model_nonlinear.png)
**Figure 3.** This shows the isolated time trends on the logit scale for delta and omicron. Through 19 weeks, the waning trends are essentially identical. From weeks 20-24 on, delta waning appears to slow down while omicron waning continues to fall linearly in logit space (and thus exponentially in vaccine efficacy). 

To bring in immune correlate data from boosted populations and independently verify that the difference in waning rates estimate above is consistent with immunological evidence and in boosted populations, I modeled the relationship between efficacy and antibody titer using the logit-linear relationship first demonstrated in [Khoury et al 2021](https://www.nature.com/articles/s41591-021-01377-8). With some algebra, one can derive a parameter-free prediction for the relative difference in NAb waning between omicron and delta from the observed VE difference. 

Starting with the fundamental relationship,
![logit(VE) \sim log(NAb)](https://latex.codecogs.com/svg.latex?\textrm{logit}(V\!E)&space;\sim&space;\log(N\!Ab)) 
one can derive a difference-in-difference prediction of the expected ratio of waning ratios from 1 to M months after vaccination using only the time trends of VE and zero free parameters:

![relative titer drop](https://latex.codecogs.com/svg.latex?\frac{N\!Ab_{omicron}^{(1)}/N\!Ab_{omicron}^{(M)}}{N\!Ab_{delta}^{(1)}/N\!Ab_{delta}^{(M)}}&space;=&space;\exp\left((\textrm{logit}(V\!E_{omicron}^{(1)})-\textrm{logit}(V\!E_{omicron}^{(M)}))&space;-&space;(\textrm{logit}(V\!E_{delta}^{(1)})-\textrm{logit}(V\!E_{delta}^{(M)})\right).)

When this ratio of titers, which depends only on differences within variant and not across variants, is greater than 1, it indicates that omicron is waning at a faster rate than delta over the time interval. 

This prediction for the ratio of antibody titer waning ratios is compared to data from two relevant papers. [Zhao et al 2022](https://t.co/mswdI6Lqxk)) look how titers from Zifivax (a subunit vaccine in trials in China) after a booster schedule (second dose 1 month after first, third dose 4 months after second) taken 1 month after booster and 4-6 months after booster compare for delta and omicron.  Median delta titers go from 2133 at 1 month to 516 at 4-6 months, and omicron from 331 to 51.  The ratio of waning ratio is thus (2133/516)/(331/51) = 1.58. This ratio tells us that, in the point estimate, omicron titers are waning at a faster rate than delta titers over that time interval. Similarly, from [Pajon et al 2022](https://t.co/lBxMqBjWpL), we can calculate an upper bound for this ratio for the Moderna vaccine by comparing titers for WT-D614G and omicron 1 and 6 months post-boost. The ratio of reported waning ratios for omicron over delta is 6.3/2.3=2.74.  In Figure 4, we show that the Zhao estimate overlaps the predictions from VE in the time period studied, and the Pajon estimate is a reasonable upper bound.

![ratio_of_waning_ratios_model_nonlinear.png](ratio_of_waning_ratios_model_nonlinear.png)
**Figure 4.** Ratio of waning ratios for omicron relative to delta. Line and ribbon are predicted from the difference in waning trends from UKHSA efficacy estimates. Red line shows the ratio of medians for Zifivax evaluated from 1 to 4-6 months post-booster. Blue shows the upper bound for the Moderna vaccine comparing omicron to D614G.  

**Concordance between these two very different types of evidence indicates that the VE waning is primarily biological and not due to confounding, and that differences in NAb waning that may not seem significant in any one study are epidemiologically meaningful.**

## Conclusion

Early data on vaccine efficacy and neutralizing antibody titers indicate that serologic protection against omicron (and likely omicron-like variants) from our current WT-based vaccines is less stable than protection against delta, and that this difference only begins to manifest 4 months after the most recent vaccination. The very limited data are consistent with similar waning rate differences after two primary doses and boosters.  **Thus, it will be important to consider longer-term immune response when considering updates to vaccine strains. It is not sufficient to only look at immune response one month after the most recent vaccination, as this differene in durability will not be seen. **
