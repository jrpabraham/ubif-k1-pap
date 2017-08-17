clear all

import delim "K1_Data.csv", varn(1) colr(2:)

encode condition, gen(treatment)
ta treatment, gen(treat)

loc depvars "vidnum savamt msgdec"
loc covariates

fvset base 3 treatment

foreach depvar of varlist `depvars' {

    reg `depvar' i.treatment, vce(cl surveyid)

}
