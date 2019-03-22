set seed 753877
*generate an SWT with 4 groups, 10 clusters per group, 6 time points (standard design)
clear
set obs  4

gen group = _n

expand 10

gen cluster = _n

expand 5

bysort cluster: gen time = _n

gen rx=0

forvalues i = 1/4{
	replace rx=1 if group == `i' & time > 5 - `i'
}

* baseline p will be 0.6 so odds of 1.5 (log(odds)=1.9)
* Intervention will increase this to 0.7 so odds of 2.33
* We want an ICC of 0.03
* so btw would be 0.03/(0.6*0.4)=0.125 (see Turner 2001)

gen btw1 = rnormal()*0.125

*also incuding a continuous outcome:
gen btw2 = rnormal()*0.1

* 40 patients per period
gen N = ceil(exp(rnormal(2.7,0.3)))

*work out odds
*adding a period effect that decreased the probabilities by about 0.02 per period
*equivalent to a log(or)=-0.08

gen logodds = 0.41 - 0.08 * (time - 1) + 0.44 * rx + btw1

gen p = exp(logodds)/(1+exp(logodds))

gen cases = rbinomial(N,p)

expand N

gen outcome1 = 0
bysort group cluster time: replace outcome1 = 1 if _n <= cases

gen outcome2 = 0.9*rnormal() + btw2 + 0.1 * (time - 1) + 0.2 * rx

drop p logodds btw1 btw2 cases N

save "C:\Users\JENNIFER\Documents\PhD\Stata_command\progs\swexample", replace
