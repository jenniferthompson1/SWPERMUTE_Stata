set seed 47804

clear
set obs 14

gen byte lab = _n

gen byte arm1 = 0

gen byte arm2 = 0 
replace arm2 = 1 if inlist(lab, 1 ,13)==1

gen byte arm3 = 0 
replace arm3 = 1 if inlist(lab, 1 ,13,7,14)==1

gen byte arm4 = 0 
replace arm4 = 1 if inlist(lab, 1 ,13,7,14,2,12)==1

gen byte arm5 = 0 
replace arm5 = 1 if inlist(lab, 1 ,13,7,14,2,12,3,4)==1

gen byte arm6 = 0 
replace arm6 = 1 if inlist(lab, 1 ,13,7,14,2,12,3,4,5)==1 | inlist(lab,9)==1

gen byte arm7 = 0 
replace arm7 = 1 if inlist(lab, 1 ,13,7,14,2,12,3,4,5)==1 | inlist(lab,9,8,10)==1

gen byte arm8 = 1

gen varb = rnormal(0,sqrt(0.16))

reshape long arm, i(lab) j(study_month)

gen periodeffect = 0 if study_month==1
replace periodeffect = -0.19 if study_month==2
replace periodeffect = -0.04 if study_month==3
replace periodeffect = -0.31 if study_month==4
replace periodeffect = -0.24 if study_month==5
replace periodeffect = -0.07 if study_month==6
replace periodeffect = -0.20 if study_month==7
replace periodeffect = -0.12 if study_month==8

gen logodds = 0.88 + periodeffect + 0.42 * arm + varb

gen p = exp(logodds)/(1 + exp(logodds))

gen N = round(rnormal(35,10),1)

gen cases = rbinomial(N, p)

expand N

gen byte confirmed=0
bysort lab study_month: replace confirm = 1 if _n <= cases

drop varb periodeffect logodds p N cases

label variable lab "Laboratory"
label variable study_month "Study period"
label variable arm "Smear microscopy (0) or Xpert (1)"
label variable confirmed "Clinical (0) or bacterially confirmed (1)"

label define lbl_arm 0 "Smear" 1 "Xpert", replace
label values arm lbl_arm

label define lbl_conf 0 "Clinical" 1 "Confirmed", replace
label values confirmed lbl_conf

save "$ProgDir\tbconfirmed", replace
