* =============================================================
* DGP.do — Generate example datasets for ivppmlhdfe
* =============================================================
*
* Model:  y_g = exp(x'beta + psi_g) * v_g
*
* Parameters:
*   beta1   = 0.5    (true coefficient on endogenous x1)
*   beta2   = 0.3    (true coefficient on exogenous x2)
*   sigma_v = 1.5    (multiplicative error SD)
*   rho_ev  = 0.5    (endogeneity: correlation between x1 and v)
*   pi_z    = 0.8    (first-stage instrument strength)
*   sigma_a = 0.3    (FE SD)
*   rho_ax  = 0.4    (Mundlak correlation: FE correlated with x1)
*
* Endogeneity:
*   x1 = pi_z * z + e               (first stage)
*   ln(v) = -sigma_v^2/2 + sigma_v * (rho_ev * e + sqrt(1-rho_ev^2) * u)
*   => Cov(x1, v) != 0 when rho_ev != 0
*
* Instrument z ~ iid N(0,1), independent of all errors and FE.
*
* Output:
*   data/ivppmlhdfe_ClassA.dta  (N=100, T=10,   1,000 obs)
*   data/ivppmlhdfe_ClassB.dta  (Nc=30, T=20,  17,400 obs)
*   data/ivppmlhdfe_ClassC.dta  (Nc=30, T=20,  17,400 obs)
* =============================================================

clear all
set more off
set seed 1235

* --- DGP parameters ---
local beta1   = 0.5
local beta2   = 0.3
local sigma_v = 1.5
local rho_ev  = 0.5
local pi_z    = 0.8
local sigma_a = 0.3
local rho_ax  = 0.4


* =============================================================
* Class A: Individual + Time FE (alpha_i + gamma_t)
*   alpha_i correlated with x1 via Mundlak device
* =============================================================

local N = 100
local T = 10
local NT = `N' * `T'

drop _all
set obs `NT'

gen int id   = ceil(_n / `T')
gen int year = mod(_n - 1, `T') + 1

* Shocks
gen double z  = rnormal()
gen double x2 = rnormal() * 0.5
gen double e  = rnormal()
gen double w  = rnormal()

* First stage
gen double x1 = `pi_z' * z + e

* Individual FE correlated with x1 (Mundlak)
bysort id: egen double mean_x1 = mean(x1)
gen double xi = rnormal() * `sigma_a'
bysort id: replace xi = xi[1]
bysort id: replace mean_x1 = mean_x1[1]
gen double alpha_i = `rho_ax' * mean_x1 + xi

* Time FE
gen double gamma_t = rnormal() * 0.1
bysort year: replace gamma_t = gamma_t[1]

* Linear predictor
gen double eta = `beta1' * x1 + `beta2' * x2 + alpha_i + gamma_t

* Multiplicative error (endogenous)
gen double lnv = -`sigma_v'^2 / 2 + `sigma_v' * (`rho_ev' * e + sqrt(1 - `rho_ev'^2) * w)
gen double y = exp(eta) * exp(lnv)

keep id year y x1 x2 z
save "ivppmlhdfe_ClassA.dta", replace


* =============================================================
* Class B: Two-Way Gravity FE (alpha_it + gamma_jt)
*   Exporter-year and importer-year FE
* =============================================================

local Nc = 30
local T  = 20
local N_pairs = `Nc' * (`Nc' - 1)
local NT_g = `N_pairs' * `T'

drop _all
set obs `NT_g'

* Panel structure: directed pairs x time
gen int pair = ceil(_n / `T')
gen int year = mod(_n - 1, `T') + 1

* Recover exp and imp from pair index
gen int exp = 1 + floor((pair - 1) / (`Nc' - 1))
gen int imp_idx = mod(pair - 1, `Nc' - 1)
gen int imp = cond(imp_idx >= exp - 1, imp_idx + 2, imp_idx + 1)
drop imp_idx

* Shocks
gen double z  = rnormal()
gen double x2 = rnormal() * 0.5
gen double e  = rnormal()
gen double w  = rnormal()
gen double x1 = `pi_z' * z + e

* Exporter-year FE
egen int exp_year = group(exp year)
gen double delta_it = rnormal() * `sigma_a'
bysort exp_year: replace delta_it = delta_it[1]

* Importer-year FE
egen int imp_year = group(imp year)
gen double gamma_jt = rnormal() * `sigma_a'
bysort imp_year: replace gamma_jt = gamma_jt[1]

* Linear predictor
gen double eta = `beta1' * x1 + `beta2' * x2 + delta_it + gamma_jt

* Multiplicative error
gen double lnv = -`sigma_v'^2 / 2 + `sigma_v' * (`rho_ev' * e + sqrt(1 - `rho_ev'^2) * w)
gen double y = exp(eta) * exp(lnv)

keep exp imp year pair y x1 x2 z
save "ivppmlhdfe_ClassB.dta", replace


* =============================================================
* Class C: Three-Way Gravity FE (alpha_ij + alpha_it + gamma_jt)
*   Pair FE correlated with x1 via Mundlak + directional FE
* =============================================================

local Nc = 30
local T  = 20
local N_pairs = `Nc' * (`Nc' - 1)
local NT_g = `N_pairs' * `T'

drop _all
set obs `NT_g'

* Panel structure
gen int pair = ceil(_n / `T')
gen int year = mod(_n - 1, `T') + 1

gen int exp = 1 + floor((pair - 1) / (`Nc' - 1))
gen int imp_idx = mod(pair - 1, `Nc' - 1)
gen int imp = cond(imp_idx >= exp - 1, imp_idx + 2, imp_idx + 1)
drop imp_idx

* Shocks
gen double z  = rnormal()
gen double x2 = rnormal() * 0.5
gen double e  = rnormal()
gen double w  = rnormal()
gen double x1 = `pi_z' * z + e

* Pair FE correlated with x1 (Mundlak)
bysort pair: egen double mean_x1_pair = mean(x1)
gen double xi_pair = rnormal() * `sigma_a'
bysort pair: replace xi_pair = xi_pair[1]
bysort pair: replace mean_x1_pair = mean_x1_pair[1]
gen double alpha_ij = `rho_ax' * mean_x1_pair + xi_pair

* Exporter-year FE
egen int exp_year = group(exp year)
gen double delta_it = rnormal() * 0.1
bysort exp_year: replace delta_it = delta_it[1]

* Importer-year FE
egen int imp_year = group(imp year)
gen double gamma_jt = rnormal() * 0.1
bysort imp_year: replace gamma_jt = gamma_jt[1]

* Linear predictor
gen double eta = `beta1' * x1 + `beta2' * x2 + alpha_ij + delta_it + gamma_jt

* Multiplicative error
gen double lnv = -`sigma_v'^2 / 2 + `sigma_v' * (`rho_ev' * e + sqrt(1 - `rho_ev'^2) * w)
gen double y = exp(eta) * exp(lnv)

keep exp imp year pair y x1 x2 z
save "ivppmlhdfe_ClassC.dta", replace

di _n
di "DGP.do complete. Datasets saved."
