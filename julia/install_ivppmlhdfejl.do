* install_ivppmlhdfejl.do
* One-time setup for the ivppmlhdfejl Julia backend
*
* Usage:
*   1. Unzip the ivppmlhdfe folder
*   2. In Stata:
*        cd "C:/your_path/ivppmlhdfe"
*        do julia/install_ivppmlhdfejl.do
*
* Prerequisites: Julia 1.9+, ssc install julia

* Source: julia/ subfolder (relative to CWD)
local src "julia"

* Destination: Stata's PLUS ado directory
local dest "`c(sysdir_plus)'"

* Download ado files
copy "`src'/ivppmlhdfejl.ado" "`dest'i/ivppmlhdfejl.ado", replace
copy "`src'/ivppmlhdfejl_load.ado" "`dest'i/ivppmlhdfejl_load.ado", replace
copy "`src'/ivppmlhdfejl_project.toml" "`dest'i/ivppmlhdfejl_project.toml", replace

* Create subdirectory for Julia package
cap mkdir "`dest'IVPPMLFixedEffectModels"
cap mkdir "`dest'IVPPMLFixedEffectModels/src"

* Copy Julia package files
copy "`src'/IVPPMLFixedEffectModels/Project.toml" "`dest'IVPPMLFixedEffectModels/Project.toml", replace
copy "`src'/IVPPMLFixedEffectModels/src/IVPPMLFixedEffectModels.jl" "`dest'IVPPMLFixedEffectModels/src/IVPPMLFixedEffectModels.jl", replace

di as txt ""
di as txt "ivppmlhdfejl v0.9.4 installed successfully."
di as txt "Files placed in: `dest'"
di as txt ""
di as txt "Requirements: Julia 1.9+ and {stata ssc install julia}"
di as txt "First run will precompile Julia packages (~30 seconds one-time cost)."
