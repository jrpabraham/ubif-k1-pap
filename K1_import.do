clear all

adopath + "."

import delim "K1_Data.csv", varn(1) colr(2:)

encode condition, gen(treatment)
ta treatment, gen(treat)

gen temp = sesunemp == "TRUE" if ~mi(sesunemp)
drop sesunemp
ren temp sesunemp

fvset base 3 treatment

/* Balance tests */

loc balancevars "socfem socpri socage sesunemp socincwinsln socconwinsln socsav"

foreach depvar of varlist `balancevars' {

    reg `depvar' i.treatment, vce(cl surveyid)
    test 1.treatment = 2.treatment

}

/* Primary hypothesis tests with covariates */

loc outcomes "vidnum savamt msgdec"

foreach depvar of varlist `outcomes' {

    reg `depvar' i.treatment, vce(cl surveyid)
    test 1.treatment = 2.treatment

    reg `depvar' i.treatment##c.socfemc i.treatment##c.socpric i.treatment##c.socagec i.treatment##c.sesunempc, vce(cl surveyid)
    test 1.treatment = 2.treatment

}

/* Primary hypothesis tests with MI correction */

mat H1 = J(3, 1, .)
mat H2 = J(3, 1, .)
mat H3 = J(3, 1, .)
loc row = 1

foreach depvar of varlist `outcomes' {

    reg `depvar' i.treatment, vce(cl surveyid)

    test 1.treatment = 0
    mat H1[`row', 1] = r(p)

    test 2.treatment = 0
    mat H2[`row', 1] = r(p)

    test 1.treatment = 2.treatment
    mat H3[`row', 1] = r(p)

    loc ++row

}

minq H1
minq H2
minq H3

/* Primary hypothesis tests with randomization inference */

/* Secondary hypothesis tests with covariates */

loc outcomes "selscorez stiscorez affscorez sesladnowz seslady2z"

foreach depvar of varlist `outcomes' {

    reg `depvar' i.treatment, vce(cl surveyid)
    test 1.treatment = 2.treatment

    reg `depvar' i.treatment##c.socfemc i.treatment##c.socpric i.treatment##c.socagec i.treatment##c.sesunempc, vce(cl surveyid)
    test 1.treatment = 2.treatment

}
