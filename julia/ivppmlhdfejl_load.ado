*! ivppmlhdfejl_load  0.10.0 19nov2026
*! Loads Julia environment and IVPPMLFixedEffectModels package
*! Pattern follows reghdfejl_load.ado (jl SetEnv + using)
*! v0.10.0: re-develops the package when the environment points at another copy
*! (e.g. an older version found earlier on the adopath), and checks the version.

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

    * Dev-install the local IVPPMLFixedEffectModels package (the copy found on the
    * adopath). v0.10.0: if the environment already points at a DIFFERENT copy,
    * re-develop this one, so the version that ships with this ado is the one loaded.
    qui findfile "IVPPMLFixedEffectModels/Project.toml"
    local pkgpath `r(fn)'
    local pkgpath = subinstr(`"`pkgpath'"', "\", "/", .)
    local pkgdir  = subinstr(`"`pkgpath'"', "/Project.toml", "", .)
    _jl: import Pkg; let _u = get(Pkg.project().dependencies, "IVPPMLFixedEffectModels", nothing); _src = _u === nothing ? nothing : Pkg.dependencies()[_u].source; if _src === nothing || lowercase(normpath(_src)) != lowercase(normpath(raw"`pkgdir'")); Pkg.develop(path=raw"`pkgdir'"); end; end;

    * Load packages (precompiled — near-zero JIT after first install)
    _jl: using IVPPMLFixedEffectModels;

    * Version check: the bridge needs the v0.10.0 keywords (accel, ffirst, weight_type)
    tempname _v
    _jl: st_numscalar("`_v'", isdefined(IVPPMLFixedEffectModels, :IVPPML_VERSION) && IVPPMLFixedEffectModels.IVPPML_VERSION == "0.10.0" ? 1.0 : 0.0);
    if scalar(`_v') != 1 {
      di as err "ivppmlhdfejl 0.10.0 needs IVPPMLFixedEffectModels 0.10.0, but another version is loaded"
      di as err "  (package found at `pkgdir'). Restart Stata so Julia reloads the package."
      exit 198
    }

    global ivppmlhdfejl_loaded 1
    if c(noisily) di as txt "(ivppmlhdfejl: IVPPMLFixedEffectModels 0.10.0 loaded from `pkgdir')"
  }
end
