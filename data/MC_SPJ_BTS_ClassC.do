* ============================================================
* MC_SPJ_BTS_ClassC.do
* 8-Panel SPJ + Bootstrap for Class C (three-way gravity FE)
* beta_8p = 4*beta_full - 2*mean(country) - 2*mean(time) + mean(8cell)
* ============================================================

use "ivppmlhdfe_ClassC.dta", clear

local B = 1000
local fe "exp#imp exp#year imp#year"
local vce "vce(cluster pair)"

* --- Step 1: Full-sample estimate ---
ivppmlhdfe y x2 (x1 = z), absorb(`fe') `vce'
local b_full = _b[x1]

* --- Step 2: Country-split halves ---
preserve
keep exp
duplicates drop
gen byte half_c = (runiform() < 0.5)
rename exp cid
tempfile ch
save `ch'
restore
rename exp cid
merge m:1 cid using `ch', nogenerate keep(match)
rename cid exp
rename half_c half_exp
rename imp cid
merge m:1 cid using `ch', nogenerate keep(match)
rename cid imp
rename half_c half_imp

* --- Step 3: Time-split half ---
summ year, meanonly
local T_mid = floor((r(min) + r(max)) / 2)
gen byte half_t = (year <= `T_mid')

* --- Step 4: 4 country subpanels (full T) ---
local b_country = 0
local spj_ok = 1
forvalues ee = 0/1 {
    forvalues ii = 0/1 {
        capture qui ivppmlhdfe y x2 (x1 = z) if half_exp == `ee' & half_imp == `ii', ///
            absorb(`fe') `vce'
        if _rc != 0 local spj_ok = 0
        if `spj_ok' local b_country = `b_country' + _b[x1]
    }
}
if `spj_ok' local b_country = `b_country' / 4

* --- Step 5: 2 time subpanels (full N) ---
local b_time = 0
if `spj_ok' {
    forvalues tt = 0/1 {
        capture qui ivppmlhdfe y x2 (x1 = z) if half_t == `tt', absorb(`fe') `vce'
        if _rc != 0 local spj_ok = 0
        if `spj_ok' local b_time = `b_time' + _b[x1]
    }
}
if `spj_ok' local b_time = `b_time' / 2

* --- Step 6: 8 cross subpanels ---
local b_8cell = 0
local n_8 = 0
if `spj_ok' {
    forvalues ee = 0/1 {
        forvalues ii = 0/1 {
            forvalues tt = 0/1 {
                if `spj_ok' {
                    capture qui ivppmlhdfe y x2 (x1 = z) ///
                        if half_exp == `ee' & half_imp == `ii' & half_t == `tt', ///
                        absorb(`fe') `vce'
                    if _rc != 0 {
                        local spj_ok = 0
                    }
                    else {
                        local b_8cell = `b_8cell' + _b[x1]
                        local n_8 = `n_8' + 1
                    }
                }
            }
        }
    }
}
if `spj_ok' local b_8cell = `b_8cell' / `n_8'

* --- Step 7: 8-panel SPJ ---
if `spj_ok' {
    local b_spj = 4 * `b_full' - 2 * `b_country' - 2 * `b_time' + `b_8cell'
    di "8p-SPJ estimate: `b_spj'"
}
else {
    local b_spj = .
    di "8p-SPJ failed: subpanel too small. Try larger Nc."
}

* --- Step 8: Bootstrap SE ---
tempfile bootresults
postfile pf b_spj_boot using `bootresults', replace

forvalues b = 1/`B' {
    preserve

    * Resample directed pairs with replacement
    bsample, cluster(pair) idcluster(new_pair)
    drop pair
    rename new_pair pair

    * Re-split countries on bootstrap sample (no preserve/restore)
    tempvar tag1
    gen byte `tag1' = 0
    bysort exp: replace `tag1' = 1 if _n == 1
    gen byte half_exp_b = (runiform() < 0.5) if `tag1' == 1
    bysort exp: replace half_exp_b = half_exp_b[1]

    tempvar tag2
    gen byte `tag2' = 0
    bysort imp: replace `tag2' = 1 if _n == 1
    gen byte half_imp_b = (runiform() < 0.5) if `tag2' == 1
    bysort imp: replace half_imp_b = half_imp_b[1]

    gen byte half_t_b = (year <= `T_mid')

    * Full sample on bootstrap data
    local bb_ok = 1

    capture noisily qui ivppmlhdfe y x2 (x1 = z), absorb(`fe') `vce'
    if _rc != 0 local bb_ok = 0
    if `bb_ok' local bb_full = _b[x1]

    * 4 country subpanels
    local bb_country = 0
    if `bb_ok' {
        forvalues ee = 0/1 {
            forvalues ii = 0/1 {
                if `bb_ok' {
                    capture noisily qui ivppmlhdfe y x2 (x1 = z) ///
                        if half_exp_b == `ee' & half_imp_b == `ii', ///
                        absorb(`fe') `vce'
                    if _rc != 0 {
                        local bb_ok = 0
                    }
                    else {
                        local bb_country = `bb_country' + _b[x1]
                    }
                }
            }
        }
    }
    if `bb_ok' local bb_country = `bb_country' / 4

    * 2 time subpanels
    local bb_time = 0
    if `bb_ok' {
        forvalues tt = 0/1 {
            if `bb_ok' {
                capture noisily qui ivppmlhdfe y x2 (x1 = z) ///
                    if half_t_b == `tt', absorb(`fe') `vce'
                if _rc != 0 {
                    local bb_ok = 0
                }
                else {
                    local bb_time = `bb_time' + _b[x1]
                }
            }
        }
    }
    if `bb_ok' local bb_time = `bb_time' / 2

    * 8 cross subpanels
    local bb_8cell = 0
    local bb_n8 = 0
    if `bb_ok' {
        forvalues ee = 0/1 {
            forvalues ii = 0/1 {
                forvalues tt = 0/1 {
                    if `bb_ok' {
                        capture noisily qui ivppmlhdfe y x2 (x1 = z) ///
                            if half_exp_b == `ee' & half_imp_b == `ii' & half_t_b == `tt', ///
                            absorb(`fe') `vce'
                        if _rc != 0 {
                            local bb_ok = 0
                        }
                        else {
                            local bb_8cell = `bb_8cell' + _b[x1]
                            local bb_n8 = `bb_n8' + 1
                        }
                    }
                }
            }
        }
    }

    if `bb_ok' {
        local bb_8cell = `bb_8cell' / `bb_n8'
        local bb_spj = 4 * `bb_full' - 2 * `bb_country' - 2 * `bb_time' + `bb_8cell'
        post pf (`bb_spj')
    }

    restore
}
postclose pf

use `bootresults', clear
summ b_spj_boot
local se_boot = r(sd)
_pctile b_spj_boot, p(2.5 97.5)
di "8p-SPJ estimate:  `b_spj'"
di "Bootstrap SE:     `se_boot'"
di "95% CI:           [`=r(r1)', `=r(r2)']"
