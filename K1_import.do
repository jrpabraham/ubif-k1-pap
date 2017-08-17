clear all

import delim "K1_Data.csv", varn(1) colr(2:)

encode condition, gen(treatment)
ta treatment, gen(treat)

gen temp = sesunemp == "TRUE" if ~mi(sesunemp)
drop sesunemp
ren temp sesunemp

fvset base 3 treatment

loc balancevars "socfem socpri socage sesunemp socincwinsln socconwinsln socsav"

foreach depvar of varlist `balancevars' {

    reg `depvar' i.treatment, vce(cl surveyid)
    test 1.treatment = 2.treatment

}


loc outcomes "vidnum savamt msgdec"

foreach depvar of varlist `outcomes' {

    reg `depvar' i.treatment, vce(cl surveyid)
    test 1.treatment = 2.treatment

    reg `depvar' i.treatment##c.socfemc i.treatment##c.socpric i.treatment##c.socagec i.treatment##c.sesunempc, vce(cl surveyid)
    test 1.treatment = 2.treatment

}
