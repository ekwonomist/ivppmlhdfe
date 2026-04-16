* ============================================================
* MC_SPJ_BTS_ClassA.do
* SPJ + Bootstrap for Class A (individual + time FE)
* FVW 2016: beta_SPJ = 3*beta_full - beta_T/2 - beta_N/2
* ============================================================

use "ivppmlhdfe_ClassA.dta", clear

local B = 200

* --- Step 1: Full-sample estimate ---
ivppmlhdfe y x2 (x1 = z), absorb(id year) vce(robust)
local b_full = _b[x1]

* --- Step 2: Time-split (split T into two halves) ---
summ year, meanonly
local T_mid = floor((r(min) + r(max)) / 2)

ivppmlhdfe y x2 (x1 = z) if year <= `T_mid', absorb(id year) vce(robust)
local b_T1 = _b[x1]

ivppmlhdfe y x2 (x1 = z) if year > `T_mid', absorb(id year) vce(robust)
local b_T2 = _b[x1]

local b_Thalf = (`b_T1' + `b_T2') / 2

* --- Step 3: Individual-split (split N into two halves) ---
tempvar half_id
bysort id: gen byte `half_id' = (runiform() < 0.5) if _n == 1
bysort id: replace `half_id' = `half_id'[1]

ivppmlhdfe y x2 (x1 = z) if `half_id' == 1, absorb(id year) vce(robust)
local b_N1 = _b[x1]

ivppmlhdfe y x2 (x1 = z) if `half_id' == 0, absorb(id year) vce(robust)
local b_N2 = _b[x1]

local b_Nhalf = (`b_N1' + `b_N2') / 2

* --- Step 4: SPJ point estimate ---
local b_spj = 3 * `b_full' - `b_Thalf' - `b_Nhalf'
di "SPJ estimate: `b_spj'"

* --- Step 5: Bootstrap SE ---
tempfile bootresults
postfile pf b_spj_boot using `bootresults', replace

forvalues b = 1/`B' {
    preserve

    * Resample N individuals with replacement
    bsample, cluster(id) idcluster(new_id)
    drop id
    rename new_id id

    * Individual-split on bootstrap sample
    tempvar bhalf
    bysort id: gen byte `bhalf' = (runiform() < 0.5) if _n == 1
    bysort id: replace `bhalf' = `bhalf'[1]

    * Recompute SPJ (capture errors from ivppmlhdfe)
    local bb_ok = 1
    capture noisily qui ivppmlhdfe y x2 (x1 = z), absorb(id year) vce(robust)
    if _rc != 0 local bb_ok = 0
    if `bb_ok' local bb_full = _b[x1]

    if `bb_ok' {
        capture noisily qui ivppmlhdfe y x2 (x1 = z) if year <= `T_mid', absorb(id year) vce(robust)
        if _rc != 0 local bb_ok = 0
    }
    if `bb_ok' local bb_T1 = _b[x1]

    if `bb_ok' {
        capture noisily qui ivppmlhdfe y x2 (x1 = z) if year > `T_mid', absorb(id year) vce(robust)
        if _rc != 0 local bb_ok = 0
    }
    if `bb_ok' local bb_T2 = _b[x1]

    if `bb_ok' {
        capture noisily qui ivppmlhdfe y x2 (x1 = z) if `bhalf' == 1, absorb(id year) vce(robust)
        if _rc != 0 local bb_ok = 0
    }
    if `bb_ok' local bb_N1 = _b[x1]

    if `bb_ok' {
        capture noisily qui ivppmlhdfe y x2 (x1 = z) if `bhalf' == 0, absorb(id year) vce(robust)
        if _rc != 0 local bb_ok = 0
    }
    if `bb_ok' local bb_N2 = _b[x1]

    if `bb_ok' {
        local bb_Thalf = (`bb_T1' + `bb_T2') / 2
        local bb_Nhalf = (`bb_N1' + `bb_N2') / 2
        local bb_spj = 3 * `bb_full' - `bb_Thalf' - `bb_Nhalf'
        post pf (`bb_spj')
    }

    restore
}
postclose pf

* --- Step 6: Compute SE and CI ---
use `bootresults', clear
summ b_spj_boot
local se_boot = r(sd)
_pctile b_spj_boot, p(2.5 97.5)
local ci_lo = r(r1)
local ci_hi = r(r2)

di "SPJ estimate:  `b_spj'"
di "Bootstrap SE:  `se_boot'"
di "95% CI:        [`ci_lo', `ci_hi']"
