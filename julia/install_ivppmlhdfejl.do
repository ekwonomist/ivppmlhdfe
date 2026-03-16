* install_ivppmlhdfejl.do
* One-time setup for the ivppmlhdfejl Julia backend
* Run this after: net install ivppmlhdfe, from("https://raw.githubusercontent.com/ekwonomist/ivppmlhdfe/main/") replace

local url "https://raw.githubusercontent.com/ekwonomist/ivppmlhdfe/main/julia"

* Destination: Stata's PLUS ado directory
local dest "`c(sysdir_plus)'"

* Download ado files
copy "`url'/ivppmlhdfejl.ado" "`dest'i/ivppmlhdfejl.ado", replace
copy "`url'/ivppmlhdfejl_load.ado" "`dest'i/ivppmlhdfejl_load.ado", replace
copy "`url'/ivppmlhdfejl_project.toml" "`dest'i/ivppmlhdfejl_project.toml", replace

* Create subdirectory for Julia package
cap mkdir "`dest'IVPPMLFixedEffectModels"
cap mkdir "`dest'IVPPMLFixedEffectModels/src"

* Download Julia package files
copy "`url'/IVPPMLFixedEffectModels/Project.toml" "`dest'IVPPMLFixedEffectModels/Project.toml", replace
copy "`url'/IVPPMLFixedEffectModels/src/IVPPMLFixedEffectModels.jl" "`dest'IVPPMLFixedEffectModels/src/IVPPMLFixedEffectModels.jl", replace

di as txt ""
di as txt "ivppmlhdfejl installed successfully."
di as txt "Files placed in: `dest'"
di as txt ""
di as txt "Requirements: Julia 1.9+ and {stata ssc install julia}"
