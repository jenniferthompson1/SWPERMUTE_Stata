
********************************************************************************
* Compare permute with swpermute and check setting seed
********************************************************************************

set more off

webuse cancer, clear
gen id = _n
gen time = 1

keep if drug==1 | drug == 2

replace drug = drug - 1

noi di _n "run a poisson model:"
poisson _d drug, exposure(_t) 
scalar regress_b = _b[drug]
matrix temp = r(table)
scalar regress_p = temp[4,2]

noi di _n "Run with cluster permute"
set seed 300
swpermute _b[drug], cluster(id) intervention(drug) period(time) reps(500) : poisson _d drug, exposure(_t) 
matrix temp = r(p)
scalar sw_p = temp[1,4]
scalar sw_b = r(obs_value)

assert reldif(regress_b,sw_b) < 1e-14


noi di _n "Check that permute gives a similar result"
set seed 300
permute drug z=_b[drug], reps(500) : poisson _d drug, exposure(_t) 
matrix temp = r(p)
scalar p_p = temp[1,1]
matrix temp = r(b)
scalar p_b = temp[1,1]

assert reldif(sw_p,p_p)<0.1



*repeat with the same seed
noi di _n "Repeat cluster permute with the same seed set to check this repeats as before"
set seed 300
swpermute _b[drug], cluster(id) intervention(drug) period(time) reps(500) : poisson _d drug, exposure(_t) 
matrix temp = r(p)
scalar sw1_p = temp[1,4]
scalar sw1_b = r(obs_value)

assert reldif(sw_p,sw1_p) < 1e-14
assert reldif(regress_b,sw1_b) < 1e-14


*check specifiying seed in the command itself
noi di _n "Specify seed within the command"
swpermute _b[drug], cluster(id) intervention(drug) period(time) reps(500) seed(300): poisson _d drug, exposure(_t) 
matrix temp = r(p)
scalar sw2_p = temp[1,4]
scalar sw2_b = r(obs_value)

assert reldif(sw_p,sw2_p) < 1e-14
assert reldif(regress_b,sw2_b) < 1e-14


*check it changes with a different seed
noi di _n "Specify a different seed that should give a different answer"
swpermute _b[drug], cluster(id) intervention(drug) period(time) reps(1000) seed(315670): poisson _d drug, exposure(_t) 
matrix temp = r(p)
scalar sw3_p = temp[1,4]
scalar sw3_b = r(obs_value)


assert reldif(regress_b,sw3_b) < 1e-14
assert reldif(sw2_p,sw3_b) > 1e-14




********************************************************************************
* PKdata
********************************************************************************

webuse pkdata3, clear
replace treat = treat-1
*note, ignoring carryover

gen strata = 0
replace strata = 1 if _n > _N/2
bysort id: replace strata = strata[1]


tostring id, gen(idstr)
tostring treat, gen(treatstr)
tostring period, gen(periodstr)


mixed outcome i.treat period || id:
scalar regress_all_b = _b[1.treat]
test 1.treat == -20
scalar regress_p20 = r(p)

regress outcome i.treat if period == 1
scalar regress_1_b = _b[1.treat]
scalar regress_1_se = _se[1.treat]
testparm i.treat
scalar regress_1_p = r(p)
test 1.treat == -2
scalar regress_1_p2 = r(p)

regress outcome i.treat if period == 2
scalar regress_2_b = _b[1.treat]
scalar regress_2_se = _se[1.treat]

noi di _n "Check overall result is right"
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
matrix p_pk = r(p)
assert reldif(`=`r(obs_value)'', regress_all_b) <1e-14

noi di _n "Check it doesn't change with the same seed"
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
matrix p_pk2 = r(p)
assert reldif(p_pk[1,4],p_pk2[1,4]) < 1e-14

noi di _n "Check it does change with a different seed"
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(9385482): /*
	*/	mixed outcome i.treat period || id:
matrix p_pk2 = r(p)

assert reldif(p_pk[1,4],p_pk2[1,4]) > 1e-14

noi di _n "Check no error with strata()"
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) strata(strata) reps(500) seed(35542): /*
	*/	mixed outcome i.treat period || id:


noi di _n "Check withinstep results are right- no weights"
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(1) withinperiod weightperiod(none) seed (01734): /*
	*/ mixed outcome i.treat period
	
matrix withinstep = r(obs_period)

assert reldif(withinstep[1,2], 0.5) < 1e-14
assert reldif(withinstep[2,2], 0.5) < 1e-14
assert reldif(withinstep[1,1], regress_1_b) < 1e-14
assert reldif(withinstep[2,1], regress_2_b) < 1e-14
assert reldif(`=`r(obs_value)'', (regress_2_b+regress_1_b)/2) < 1e-14

noi di _n "Check withinstep results are right- number of cluster weights"

swpermute _b[1.treat] , /*
	*/ cluster(id) intervention(treat) period(period) reps(1) withinperiod weightperiod(N): /*
	*/ regress outcome i.treat 
matrix withinstep = r(obs_period)

assert reldif(withinstep[1,2], 0.5) < 1e-14
assert reldif(withinstep[2,2], 0.5) < 1e-14
assert reldif(withinstep[1,1], regress_1_b) < 1e-14
assert reldif(withinstep[2,1], regress_2_b) < 1e-14
assert reldif(`=`r(obs_value)'', /*
	*/ (withinstep[1,2] * withinstep[1,1] + withinstep[2,2] * withinstep[2,1]) / (withinstep[1,2] + withinstep[2,2])) < 1e-14

noi di _n "Check withinstep results are right- variance weights"

swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(1) withinperiod weightperiod(variance   _se[1.treat]^2): /*
	*/ regress outcome i.treat 
matrix withinstep = r(obs_period)

assert reldif(withinstep[1,1], regress_1_b) < 1e-14
assert reldif(withinstep[2,1], regress_2_b) < 1e-14
assert reldif(withinstep[1,2], (regress_1_se^-2)/((regress_1_se^-2)+(regress_2_se^-2))) < 1e-14
assert reldif(withinstep[2,2], (regress_2_se^-2)/((regress_1_se^-2)+(regress_2_se^-2))) < 1e-14
assert reldif(`=`r(obs_value)'', /*
	*/ (withinstep[1,2] * withinstep[1,1] + withinstep[2,2] * withinstep[2,1]) / (withinstep[1,2] + withinstep[2,2])) < 1e-14

* Check p values
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(737027) null(0 -2) outcome(outcome): /*
	*/ regress outcome i.treat if period == 1

matrix p = r(p)

assert p[1,5] < p[1,4] & p[1,4] < p[1,6]
assert p[2,5] < p[2,4] & p[2,4] < p[2,6]

** P-vale should be close to the regression value
assert reldif(p[1,4], regress_1_p ) < 0.1
assert reldif(p[2,4], regress_1_p2 ) < 0.1


*check right and left
 noi di _n "Check right and left p values are different"
 
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382) left: /*
	*/	mixed outcome i.treat period || id:
matrix p_pkleft = r(p)
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382) right: /*
	*/	mixed outcome i.treat period || id:
matrix p_pkright = r(p)

assert p_pkleft[1,2] != p_pkright[1,2] 
assert p_pkleft[1,2] != p_pk[1,2]
assert p_pk[1,2] != p_pkright[1,2] 


*check changing the level changes the CI
 noi di _n "Check right and left p values are different and cumm to the total"
 
swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382) level(90): /*
	*/	mixed outcome i.treat period || id:
matrix p_pk2 = r(p)

assert p_pk2[1,5] < p_pk[1,6]

********************************************************************************
* Check error messages
********************************************************************************

noi di _n _dup(79) "="

noi di _n "Check error messages"

noi di "Invalid command given"
cap noi swpermute _b[1.treat], /*
	*/	cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/ notacommand outcome i.treat period
assert _rc!=0

noi di "No command given"
cap noi swpermute _b[1.treat], /*
	*/	cluster(id) intervention(treat) period(period) reps(500) seed(382): 
assert _rc!=0
	
noi di _n "Invalid statistic given"
cap noi swpermute _b[nota1.treat], /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0
	
noi di _n "No statistic given"
cap noi swpermute , /*
	*/ cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0
	
noi di _n "Strata is cluster variable"
cap noi swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) strata(id) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0
	
noi di _n "Strata is by variable"
cap noi swpermute  _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) strata(treat) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0
	
noi di _n "Strata is time variable"
cap noi swpermute  _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) strata(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

	
noi di _n "Strata is something else"
cap noi swpermute  _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(period) strata(5) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "No variance given for withinstep"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382) withinperiod weightperiod(variance): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "2 statistics given"
cap noi swpermute _b[1.treat] _se[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Cluster variable is a string"
cap noi swpermute _b[1.treat], /*
	*/ cluster(idstr) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Cluster missing"
cap noi swpermute _b[1.treat], /*
	*/  intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Cluster not a variable"
cap noi swpermute _b[1.treat], /*
	*/  cluster(otherword) intervention(treat) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Time variable is a string"
cap noi swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treat) period(periodstr) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Time variable is a non integer numeric variable- should run"
preserve
	replace period = period + 0.5
	swpermute _b[1.treat], /*
		*/ cluster(id) intervention(treat) period(period) reps(500) seed(382): /*
		*/	mixed outcome i.treat period || id:
restore

noi di _n "Time missing"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Time not a variable"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(otherword) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "By variable is a string"
cap noi swpermute _b[1.treat], /*
	*/ cluster(id) intervention(treatstr) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "By missing"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) period(period) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "By not a variable"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(periodtest) reps(500) seed(382): /*
	*/	mixed outcome i.treat period || id:
assert _rc!=0

noi di _n "Check periods with all clusters in one condition have a warning that they aren't included in a within step analysis"
preserve
	replace treat = 0 if period == 1
	cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382) withinperiod weightperiod(N): /*
	*/	regress outcome i.treat 
restore

noi di _n "An invalid period weight given"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382) withinperiod weightperiod(sequence): /*
	*/	regress outcome i.treat 
assert _rc != 0

noi di _n "Check abbreviation of variance works"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(1) seed(382) withinperiod weightperiod(var _se[1.treat]^2): /*
	*/	regress outcome i.treat 
assert _rc == 0
local withinest1 = r(obs_value)
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(1) seed(382) withinperiod weightperiod(variance _se[1.treat]^2): /*
	*/	regress outcome i.treat 
assert reldif(`withinest1', r(obs_value))  < 1e-14
	
noi di _n "No outcome specified with null values"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382) null(0 2 3): /*
	*/	regress outcome i.treat 
assert _rc!=0

noi di _n "Invalid null values given"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(500) seed(382) null(potato) outcome(outcome): /*
	*/	regress outcome i.treat 
assert _rc!=0

noi di _n "Invalid number of reps given"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(potato) seed(382) null(potato) outcome(outcome): /*
	*/	regress outcome i.treat 
assert _rc!=0

cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) reps(38 43) seed(382) : /*
	*/	regress outcome i.treat 
assert _rc!=0

noi di _n "Invalid seed"
cap noi swpermute _b[1.treat], /*
	*/  cluster(id) intervention(treat) period(period) seed(anumber) : /*
	*/	regress outcome i.treat 
assert _rc!=0

noi di _n "intervention not constant in each cluster-period"
preserve
	replace treat = 1 if _n==1
	cap noi swpermute _b[1.treat], /*
		*/  cluster(id) intervention(treat) period(period) seed(anumber) : /*
		*/	regress outcome i.treat 
	assert _rc!=0
restore

********************************************************************************
* Handling of missing data
********************************************************************************


use  "$ProgDir\testing_data", clear

noi di _n "Check data is the same"
swpermute _b[rx], intervention(rx) period(time) cluster(cluster) /*
	*/ reps(50) withinperiod seed(877587): /*
	*/ logit outcome1 rx  
assert reldif(`=`r(obs_value)'', 0.2813965) < 1e-7

noi di _n "Check what happens if a period is missing"
replace time = . if time == 3
swpermute _b[rx], intervention(rx) period(time) cluster(cluster) /*
	*/ reps(50) withinperiod  seed(4844): /*
	*/ logit outcome1 rx 
assert reldif(`=`r(obs_value)'', .304327) < 1e-7
	
noi di _n "Check what happens if intervention is missing of a few individuals	"
use  "$ProgDir\testing_data", clear
replace rx = . if _n==352
replace rx = . if _n==473
swpermute _b[rx], intervention(rx) period(time) cluster(cluster) /*
	*/ reps(50) withinperiod seed(908948): /*
	*/ logit outcome1 rx  
assert reldif(`=`r(obs_value)'', .2865429) < 1e-7

noi di _n "Check what happens if intervention is missing for a whole cluster-period	"
use  "$ProgDir\testing_data", clear
replace rx = . if group == 1 & time == 5
replace rx = . if group == 2 & time == 4
replace rx = . if group == 3 & time == 3
replace rx = . if group == 4 & time == 2
swpermute _b[rx], intervention(rx) period(time) cluster(cluster) /*
	*/ reps(50) withinperiod seed(908948): /*
	*/ logit outcome1 rx  
assert reldif(`=`r(obs_value)'', 0.3378973) < 1e-7

noi di _n "check what happens if a cluster name is missing	"
use  "$ProgDir\testing_data", clear

replace cluster = . if cluster == 3

swpermute _b[rx], intervention(rx) period(time) cluster(cluster) /*
	*/ reps(50) withinperiod seed(8295421): /*
	*/ logit outcome1 rx 
assert reldif(`=`r(obs_value)'', 0.2991831) < 1e-7
