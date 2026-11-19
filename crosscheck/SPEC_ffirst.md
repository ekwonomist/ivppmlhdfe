# ffirst — first-stage statistics: exact conventions for all four ports

Stata track, 2026-11-19. Pinned by matching ivreghdfe 1.1.4 / ivreg2 4.1.12 / ranktest 2.0.04 on Stata 18.5 MP
(desktop EKWONOMIST, job `ivppml010_ref`, unit `ffirst`). Reference values: `results\stata_ffirst_reference.json`
(key `ffirst` of each case, same keys as the ffirst block of SPEC 7). A numpy implementation of the formulas below,
written from this file and not from the GPL code, is `stata\scripts\ffirst_reference.py`; it reproduces ivreghdfe's
numbers to **≤ 1.1e-14 relative on all 15 IV cases** (robust, one-way and two-way cluster, overidentified, two
endogenous regressors, pweight, fweight, offset, separation, singletons).

**In one line:** the ffirst statistics of ivppmlhdfe are exactly what
`ivreghdfe z X1 (X2 = Z) [aw = w·μ], absorb(<same FEs>) <same vce> ffirst` reports at the converged IRLS
solution, where z is the working dependent variable (z does not enter any first-stage statistic).

## 1. Inputs

At convergence (SPEC 3):
- `μ_i` the fitted mean; `ω_i = w_i·μ_i` the IRLS weight (w = user pweight, or 1; for fweights see below).
- the columns partialled out of the fixed effects under ω (accel: the refreshed columns of SPEC 2; noaccel: the
  last IRLS iteration's columns together with that iteration's weights, as the v0.9.4 VCE uses them):
  `X̃` endogenous (N_obs × Ke), `W̃` included exogenous (N_obs × Kw), `Z̃` excluded instruments (N_obs × Lz).
  No constant: it is absorbed by the fixed effects.
- `R = [Z̃, W̃]` all instruments, L = Lz + Kw columns.
- Weights in the two roles:
  - cross products use ω: `A'ΩB = Σ_i ω_i a_i b_i'`;
  - the robust meat uses ψ: `Σ_i ψ_i s_i s_i'`, with ψ_i = ω_i² (no weights, pweights);
  - **fweights** f_i (rows replicated f_i times): ω_i = f_i·μ_i, ψ_i = f_i·μ_i².
  - Any rescaling of ω cancels in every statistic (ivreg2 normalises aweights to sum to N; that does not matter).
- `N` = number of observations in the estimation sample (= e(N)); with fweights N = Σ f_i.
- Clusters: one-way G = number of clusters. Two-way (c1, c2): G = min(G1, G2), and every meat is
  M = M(c1) + M(c2) − M(c1∩c2) (Cameron-Gelbach-Miller), **no** per-term G/(G−1) factor and **no** PSD fix
  (unlike the main ivppmlhdfe VCE).
- `sdof` = the absorbed degrees of freedom e(df_a) (section 4), **or 1 when e(df_a) = 0** (ivreghdfe:
  "adjustment to match ivreg2 and old reghdfe, happens if absvar is nested in cluster"). dofminus = 0.

## 2. Meats and Wald statistics (no small-sample factor inside)

For a residual vector e (or matrix) and regressor matrix R, with per-observation score s_i = e_i·r_i:
- robust: `M = Σ_i ψ_i s_i s_i'`
- cluster: `M = Σ_g u_g u_g'`, `u_g = Σ_{i∈g} ω_i s_i` (two-way: combine as above).
- `V0 = (R'ΩR)⁻¹ M (R'ΩR)⁻¹` (HC0 / CR0).
- `W0 = π̂_Z' [V0_ZZ]⁻¹ π̂_Z`, the Wald statistic that the Lz excluded-instrument coefficients are zero.

Turning a Wald statistic W0 with q restrictions into an F statistic (ivreg2's rule; L counts the included
exogenous regressors, not the constant):

| VCE | F | df2 |
|---|---|---|
| robust | `W0 / q · (N − L − sdof) / N` | `N − L − sdof` |
| cluster | `W0 / q · (N − L − sdof)/(N − 1) · (G − 1)/G` | `G − 1` |

p-values: `Ftail(q, df2, F)` (upper tail of F); chi-squared statistics: `chi2tail(df, stat)`.

## 3. The statistics

For each endogenous regressor x_k (row k of e(first)):
1. **F, df1, df2, p.** Weighted LS of x̃_k on R: π̂ = (R'ΩR)⁻¹R'Ωx̃_k, e = x̃_k − Rπ̂; W0 as in section 2;
   F by the table with q = Lz; df1 = Lz.
2. **partialR2** = 1 − e'Ωe / e_r'Ωe_r, where e_r is the residual of x̃_k on W̃ alone (x̃_k itself if Kw = 0).
   Independent of the VCE and of sdof.
3. **Sanderson-Windmeijer (2016).** If Ke = 1: SWF = F, SWdf1 = df1, SWdf2 = df2, SWp = p. If Ke ≥ 2:
   - X̂ = R·Π̂ (first-stage fitted values of all endogenous regressors);
   - S = [X̂_{−k}, W̃] (fitted values of the OTHER endogenous regressors, plus W̃), S_act = [X̃_{−k}, W̃];
   - b1 = (S'ΩS)⁻¹ S'Ω x̃_k;  e1 = x̃_k − S_act·b1 (actual, not fitted, regressors);
   - weighted LS of e1 on R: coefficients b2, residual e2; W0 from b2's Z block with the meat of e2;
   - q = SWdf1 = Lz − Ke + 1; SWF and SWdf2 by the table; SWp = Ftail(SWdf1, SWdf2, SWF).
   (The meat ignores the estimation error in b1, as ivreg2 does.)

For the system (W̃ partialled out of X̃ and Z̃ by Ω-weighted LS: Ȳ = M_W X̃, Z̄ = M_W Z̃):
- A = Z̄'ΩZ̄, B = Ȳ'ΩȲ, C = Z̄'ΩȲ, Π̂ = A⁻¹C (Lz × Ke). Lower Cholesky factors A = L_A L_A', B = L_B L_B'.
- Θ̂ = L_A' Π̂ (L_B')⁻¹; full SVD Θ̂ = U diag(s) V' (s descending, U Lz×Lz, V Ke×Ke). s are the canonical
  correlations.
- Kleibergen-Paap (2006) rk statistic for H0: rank(Π) = kk with kk = Ke − 1:
  - partition U = [·, ·] at kk: u12 = U[1:kk, kk+1:Lz], u22 = U[kk+1:Lz, kk+1:Lz]; same for V: v12 = V[1:kk, kk+1:Ke],
    v22 = V[kk+1:Ke, kk+1:Ke] (empty blocks when kk = 0);
  - A⊥ = [u12; u22]·u22⁻¹·(u22 u22')^{1/2},  B⊥ = (v22 v22')^{1/2}·(v22')⁻¹·[v12; v22]'
    (^{1/2} = symmetric PSD square root; with kk = 0, A⊥ = I_Lz and B⊥ = I_1);
  - T = B⊥ ⊗ A⊥',  λ = T·vec(Θ̂)  (vec stacks columns), J = L_B⁻¹ ⊗ L_A⁻¹;
  - for a residual matrix v (N_obs × Ke): scores s_i = v_i ⊗ z̄_i (row-wise Kronecker: all z̄ for v_1, then for v_2,
    ...), meat M(v) as in section 2;
  - rk(v) = λ' [T J M(v) J' T']⁻¹ λ.
- **idstat** (KP rk LM, underidentification) = rk(Ȳ) (LM: the "residual" is Ȳ itself). iddf = Lz − Ke + 1;
  idp = chi2tail(iddf, idstat). Independent of sdof.
- **widstat** (KP rk Wald F, weak identification) from rk_W = rk(Ȳ − Z̄Π̂):
  robust `rk_W / N · (N − L − sdof) / Lz`; cluster `rk_W/(N − 1) · (N − L − sdof) · (G − 1)/G / Lz`.
- **cdf** (Cragg-Donald Wald F, always i.i.d.) = min_j s_j²/(1 − s_j²) · (N − sdof − L) / Lz.

With Ke = 1 the rk formulas reduce to the robust (or cluster) score test (LM) and Wald test that Π = 0, and
widstat = F.

## 4. sdof: the absorbed degrees of freedom (reghdfe's rule, = ivppmlhdfe's e(df_a))

Intercept FEs g = 1..F in the order of absorb(), G_g levels in the estimation sample; M_g = redundant levels.
1. If vce is cluster: every FE that is **nested within a cluster variable** (each FE level has a single value of
   that cluster variable; an FE identical to a cluster variable is nested) gets M_g = G_g. Call them nested.
2. Let idx = the remaining FEs, in order. If len(idx) > 1, set M_g = 1 for every g in idx except the first; if
   some FE was nested and len(idx) ≥ 1, set M_g = 1 for every g in idx including the first.
3. Pairwise mobility groups: for every pair g < h in idx, m = number of connected components of the bipartite
   graph whose nodes are the levels of FE g and of FE h and whose edges are the observations;
   M_h = max(M_h, m).
4. df_a = Σ_g (G_g − M_g); **sdof = df_a, or 1 if df_a = 0**.

This is reghdfe 5.9's `estimate_dof()` (dof(all) = pairwise clusters continuous), which ivppmlhdfe v0.9.4 already
posts as e(df_a); ivreghdfe 1.1.4 (reghdfe 6) gives the same counts on every test case and adds the "1 if 0" rule.
It is a conservative count: in class C (pair + exporter-year + importer-year, 20 countries × 10 years) the exact
rank deficiency is 49 levels, reghdfe's rule counts 40 (df_a = 740 robust). Values per case: `ffirst_meta.df_a`
and `ffirst_meta.sdof` in the reference JSON (A: 109; A cluster id: 9; B: 590 robust and cluster pair; B_sep: 589;
C robust: 740; C cluster pair: 389; C two-way exp imp: 0 → sdof 1).

## 5. Names (Stata)

- `e(first)`: one row per endogenous regressor (rownames = the variable names), columns
  `F df1 df2 p partialR2 SWF SWdf1 SWdf2 SWp` (the SPEC 7 keys). ivreg2's own e(first) is the transpose with
  more rows; the mapping is F↔F, df1↔df, df2↔df_r, p↔pvalue, partialR2↔pr2, SWF↔SWF, SWdf1↔SWFdf1,
  SWdf2↔SWFdf2, SWp↔SWFp.
- `e(idstat) e(iddf) e(idp) e(cdf) e(widstat)` as in ivreg2; also `e(sdofminus)` = sdof.
- SPEC 7 JSON: `ffirst = {"<endog>": {"F","df1","df2","p","partialR2","SWF","SWdf1","SWdf2","SWp"}, "idstat",
  "iddf", "idp", "cdf", "widstat"}`; df1/df2/iddf as numbers.

## 6. Where each reference number came from

- `ffirst` in `results\stata_ffirst_reference.json`: ivreghdfe 1.1.4's e(first) and e(idstat), e(iddf), e(idp),
  e(cdf), e(widstat), run on each IV case of `testdata\cases.json` at the v0.9.4 solution:
  `predict mu`, ω = w·μ (B_pw: `[pw=ω]`; others `[aw=ω]`; A_fw: data `expand`ed by fw, then `[aw=μ]`),
  working depvar z = ln μ − offset − 1 + y/μ, `absorb(<FEs>) <vce> ffirst tolerance(1e-12)`.
- Cross-checks in the same job: ivreg2 on the data partialled at tolerance 1e-14 (`reghdfe v [aw=ω],
  absorb(...) residuals()`, exported in %21x), `nocons`, with sdofminus(0) and sdofminus(df_a): numpy agrees
  with both to ≤ 1.1e-14; ivreghdfe and ivreg2-with-sdofminus(df_a) agree to ≤ 1.2e-14 where sdof = df_a.
- `ffirst_numpy`: the numpy values (same keys); `ffirst_meta`: N, G, L, Ke, Lz, df_a, sdof, the canonical
  correlations, and the numpy-vs-ivreghdfe maximum relative difference.
- The reference is taken at the converged μ (as the accel path's refreshed columns). The noaccel path uses the
  last iteration's lagged weights, so its ffirst numbers can differ from the reference by the same order as
  the v0.9.4 SE lag (up to about 1e-6 relative; measured per case in `stata\REPORT.md`).

## 7. Reference values (x1 row; full precision in the JSON)

| case | Ke | Lz | N | G | sdof | F (x1) | partial R2 (x1) | SW F (x1) | idstat (df) | cdf | widstat |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A_rob | 1 | 1 | 1000 | - | 109 | 286.567027 | 0.377333 | 286.567027 | 102.477582 (1) | 538.729455 | 286.567027 |
| A_cl | 1 | 1 | 1000 | 100 | 9 | 180.177782 | 0.377333 | 180.177782 | 29.970058 (1) | 599.328943 | 180.177782 |
| B_rob | 1 | 1 | 8700 | - | 590 | 2634.910052 | 0.366616 | 2634.910052 | 1049.289309 (1) | 4693.088222 | 2634.910052 |
| B_cl | 1 | 1 | 8700 | 870 | 590 | 2636.786587 | 0.366616 | 2636.786587 | 499.816387 (1) | 4693.088222 | 2636.786587 |
| C_rob | 1 | 1 | 3800 | - | 740 | 1122.188731 | 0.336482 | 1122.188731 | 562.413834 (1) | 1550.768318 | 1122.188731 |
| C_cl | 1 | 1 | 3800 | 380 | 389 | 813.857893 | 0.336482 | 813.857893 | 178.077737 (1) | 1728.766905 | 813.857893 |
| C_cl2 | 1 | 1 | 3800 | 20 | 1 | 503.868789 | 0.336482 | 503.868789 | 9.495901 (1) | 1925.528876 | 503.868789 |
| A_overid | 1 | 2 | 1000 | - | 109 | 143.728574 | 0.377377 | 143.728574 | 104.454708 (2) | 269.111936 | 143.728574 |
| C_overid | 1 | 2 | 3800 | 380 | 389 | 408.168352 | 0.338873 | 408.168352 | 178.448113 (2) | 873.416197 | 408.168352 |
| C_k2 | 2 | 3 | 3800 | 380 | 389 | 297.596474 | 0.340954 | 477.940499 | 180.963302 (2) | 475.386638 | 264.918821 |
| B_pw | 1 | 1 | 8700 | 870 | 590 | 2539.818223 | 0.365745 | 2539.818223 | 502.806274 (1) | 4675.504279 | 2539.818223 |
| A_fw | 1 | 1 | 1999 | - | 109 | 612.587968 | 0.360651 | 612.587968 | 200.790964 (1) | 1065.004982 | 612.587968 |
| C_off | 1 | 1 | 3800 | - | 740 | 1119.204551 | 0.337391 | 1119.204551 | 559.157775 (1) | 1557.086835 | 1119.204551 |
| B_sep | 1 | 1 | 8671 | - | 589 | 2493.980398 | 0.368821 | 2493.980398 | 984.534035 (1) | 4721.432150 | 2493.980398 |
| A_single | 1 | 1 | 1000 | - | 109 | 286.567027 | 0.377333 | 286.567027 | 102.477582 (1) | 538.729455 | 286.567027 |

C_k2's x3 row and every p-value are in the JSON.

## 8. Two facts the ports should know (found while pinning this)

- **B_sep counts:** Stata v0.9.4 reports N = 8671, num_singletons = 29, num_separated = 0. The 29 observations of the
  all-zero exporter-year group are dropped by reghdfe's iweight trick and counted as singletons, not as
  separated. SPEC 4 requires these counts to match exactly, so a port that reports 29 separated will FAIL
  compare.py against Stata on B_sep. (Not changed in v0.10.0: noaccel must equal v0.9.4 bit for bit.)
- **ppmlhdfe on the desktop is 2.3.3 (02nov2025)**; the copy archived in `program\archive\existing_stata_code`
  is 2.3.0. Both default to a warm-started partial-out (use_exact_partial = 0) with standardised data
  (standardize_data = 1), which v0.9.4 does not do, so v0.9.4 equals ppmlhdfe 2.3.3 to 2.3e-15 in coefficients
  but only 1.4e-10 in SEs on A_ppml / C_ppml (`results\stata_ppmlhdfe.json`), not bit for bit.
