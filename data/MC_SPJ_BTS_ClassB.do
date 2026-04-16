* ============================================================
* MC_SPJ_BTS_ClassB.do
* SPJ + Bootstrap for Class B (two-way gravity FE)
* WZ2021: beta_SPJ = 2*beta_full - mean(beta_aa, beta_ab, beta_ba, beta_bb)
* ============================================================

use "ivppmlhdfe_ClassB.dta", clear

local B = 200

* --- Step 1: Full-sample estimate ---
ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair)
local b_full = _b[x1]

* --- Step 2: Country-split (4 subpanels) ---
preserve
keep exp
duplicates drop
gen byte half = (runiform() < 0.5)
rename exp cid
tempfile country_halves
save `country_halves'
restore

rename exp cid
merge m:1 cid using `country_halves', nogenerate keep(match)
rename cid exp
rename half half_exp

rename imp cid
merge m:1 cid using `country_halves', nogenerate keep(match)
rename cid imp
rename half half_imp

* Estimate on 4 subpanels
local b_sum = 0
forvalues ee = 0/1 {
    forvalues ii = 0/1 {
        qui ivppmlhdfe y x2 (x1 = z) if half_exp == `ee' & half_imp == `ii', ///
            absorb(exp#year imp#year) vce(cluster pair)
        local b_sum = `b_sum' + _b[x1]
    }
}
local b_spj = 2 * `b_full' - `b_sum' / 4

di "SPJ estimate: `b_spj'"

* --- Step 3: Bootstrap SE ---
tempfile bootresults
postfile pf b_spj_boot using `bootresults', replace

forvalues b = 1/`B' {
    preserve

    * Resample directed pairs with replacement
    bsample, cluster(pair) idcluster(new_pair)
    drop pair
    rename new_pair pair

    * Re-split countries on bootstrap sample
    tempfile bch
    tempvar tag
    gen byte `tag' = 0
    bysort exp: replace `tag' = 1 if _n == 1
    gen byte half_exp_b = (runiform() < 0.5) if `tag' == 1
    bysort exp: replace half_exp_b = half_exp_b[1]

    tempvar tag2
    gen byte `tag2' = 0
    bysort imp: replace `tag2' = 1 if _n == 1
    gen byte half_imp_b = (runiform() < 0.5) if `tag2' == 1
    bysort imp: replace half_imp_b = half_imp_b[1]

    * Estimate full + 4 subpanels
    local bb_ok = 1

    capture noisily qui ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair)
    if _rc != 0 local bb_ok = 0
    if `bb_ok' local bb_full = _b[x1]

    local bb_sum = 0
    if `bb_ok' {
        forvalues ee = 0/1 {
            forvalues ii = 0/1 {
                if `bb_ok' {
                    capture noisily qui ivppmlhdfe y x2 (x1 = z) ///
                        if half_exp_b == `ee' & half_imp_b == `ii', ///
                        absorb(exp#year imp#year) vce(cluster pair)
                    if _rc != 0 {
                        local bb_ok = 0
                    }
                    else {
                        local bb_sum = `bb_sum' + _b[x1]
                    }
                }
            }
        }
    }

    if `bb_ok' {
        local bb_spj = 2 * `bb_full' - `bb_sum' / 4
        post pf (`bb_spj')
    }

    restore
}
postclose pf

use `bootresults', clear
summ b_spj_boot
local se_boot = r(sd)
_pctile b_spj_boot, p(2.5 97.5)
di "SPJ estimate:  `b_spj'"
di "Bootstrap SE:  `se_boot'"
di "95% CI:        [`=r(r1)', `=r(r2)']"
