cd "H:\My Documents\PhD\Stata_command\progs\paper_code"

*******************************************
*  describe data                          *
*******************************************

*bits in the paper
use TBdiagnostic, clear
sjlog using logdescribe, replace
describe

list in 1/5
sjlog close, replace


*******************************************
*  Overall analysis                       *
*******************************************

cap sjlog close
sjlog using logoverall, replace
swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ reps(1000) seed(20255) nodots: /*
	*/ melogit confirmed i.study_month arm || lab : 
sjlog close, replace

*******************************************
*  within-period example                    *
*******************************************

cap sjlog close
sjlog using logwithinnull1, replace
collapse (mean) risk_confirmed = confirmed , by(lab study_month arm)

swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000): regress risk_confirmed arm
sjlog close, replace


cap sjlog close
sjlog using logwithinnull2, replace
* Estimate standard error
display .1052611 / invnorm( 1 - 0.0500 / 2 )

* Initial lower bound of 95% CI:
display .1052611 - 1.96 * .0537
* Initial upper bound of 95% CI:
display .1052611 + 1.96 * .0537
sjlog close, replace

	
cap sjlog close
sjlog using logwithinnull3, replace
swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(0.211) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm
sjlog close, replace
	
swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(0.200) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm

swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(0.197) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm

swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(0.198) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm

swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(0.199) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm
		
swpermute _b[arm], cluster(lab) period(study_month) intervention(arm) /*
	*/ seed(9845) withinperiod weightperiod(variance _se[arm]^2) nodots /*
	*/ reps(1000) null(-0.001) outcome(risk_confirmed): /*
	*/ regress risk_confirmed arm
	

/* confidence interval is  (0, 19.8%) */


