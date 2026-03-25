*! ivppmlhdfejl_load  0.8.0 24mar2026
*! Loads Julia environment and IVPPMLFixedEffectModels package
*! Pattern follows reghdfejl_load.ado (jl SetEnv + using)

cap program drop ivppmlhdfejl_load
program define ivppmlhdfejl_load
  version 15

  if `"$ivppmlhdfejl_loaded"' == "" {
    cap jl version
    if _rc {
      di as err `"Cannot access Julia. {cmd:ivppmlhdfejl} requires the {cmd:jl} command."'
      di as err `"Install via {stata ssc install julia}."'
      di as err `"Julia must also be installed; see {help jl##installation:help jl}."'
      exit 198
    }

    * Set up Julia package environment (installs/updates deps from Project.toml)
    qui findfile ivppmlhdfejl_project.toml
    local projectfile `r(fn)'
    qui jl SetEnv @ivppmlhdfejl, update project(`projectfile')

    * Dev-install the local IVPPMLFixedEffectModels package
    qui findfile "IVPPMLFixedEffectModels/Project.toml"
    local pkgpath `r(fn)'
    local pkgpath = subinstr(`"`pkgpath'"', "\", "/", .)
    local pkgdir  = subinstr(`"`pkgpath'"', "/Project.toml", "", .)
    _jl: import Pkg; if !haskey(Pkg.project().dependencies, "IVPPMLFixedEffectModels"); Pkg.develop(path=raw"`pkgdir'"); end;

    * Load packages (precompiled — near-zero JIT after first install)
    _jl: using IVPPMLFixedEffectModels;

    global ivppmlhdfejl_loaded 1
    if c(noisily) di as txt "(ivppmlhdfejl: IVPPMLFixedEffectModels loaded)"
  }
end
