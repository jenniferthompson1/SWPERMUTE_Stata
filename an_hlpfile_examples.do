

use  TBconfirmed, clear

* use permute on mixed effect model
swpermute _b[arm], intervention(arm) period(study_month) cluster(lab) reps(100): melogit confirmed arm i.study_month || lab :

*Calculate the risk of the outcome in each cluster-period{p_end}
collapse (sum) cases=confirmed (count) N=confirmed, by(lab study_month arm)
gen risk = cases / N

*Run a within period analysis on cluster-summaries to calculate a risk difference{p_end}
swpermute r(mu_2) - r(mu_1), intervention(arm) period(study_month) cluster(lab) withinperiod weightperiod(variance r(se)^2)  reps(100): ttest risk, by(arm)

*Test different null values to generate a confidence interval{p_end}
swpermute r(mu_2) - r(mu_1), outcome(risk) null(0.05 0.15) intervention(arm) period(study_month) cluster(lab) withinperiod weightperiod(variance r(se)^2)  reps(100): ttest risk, by(arm)
