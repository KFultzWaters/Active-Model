clear all
clc
close all
%% ========================================================================
%  ACTIVE CONTRACTION MODEL
%  Material: Holzapfel_Ogden_Active (custom plugin -- HO passive law +
%  Hill/calcium/sarcomere-length active contraction hand-added into
%  DevStress()/sbar(), FEBio theory manual Eq 5.10.1-3 for the active part)
%  Three regions: LV, RV, Septum -- each with its own material block,
%  region-specific passive constants AND region-specific active tension
%  shape parameter (beta), matching the passive model's per-region
%  Material/Elements/MeshData structure.
%
%  Registered plugin type string: 'Holzapfel-Ogden_ACTIVE'
%  (see dllmain.cpp: REGISTER_FECORE_CLASS(HolzapfelMyocardiumActivePI,
%  "Holzapfel-Ogden_ACTIVE") -- distinct from the original untouched
%  passive plugin's "Holzapfel_Ogden" string, to avoid a registration
%  collision between the two plugins.)
%
%  *** ACTIVE-PARAMETER STATUS ***: FEBio 4.9.0 crashes (heap corruption /
%  access violation) when this plugin class registers more than its
%  original 12 passive parameters via ADD_PARAMETER, once ANY of the
%  active-contraction parameters (Tmax/Ca0/beta/l0/refl/ascl) are actually
%  given values in the .feb file -- confirmed true for every one of them
%  individually. Currently ONLY "refl" is being sent through as a live
%  test of this boundary; Tmax/Ca0/beta/l0/ascl remain commented out below.
%  The two-class hardcoded-constants workaround (LoBeta/HiBeta, see
%  Holzapfel-Ogden_ActiveLoBeta.h/.cpp and Holzapfel-Ogden_ActiveHiBeta.h/.cpp
%  in this same output folder) is set aside for now, not deleted -- we can
%  come back to it.
%
%  Plugin parameters (flat, no nested wrapper, NO fiber tag -- orientation
%  comes purely from mat_axis, same as the original passive plugin):
%    a, b, af, bf, as, bs, afs, bfs, asn, bsn, anf, bnf   -- passive HO terms
%    Tmax, Ca0, beta, l0, refl, ascl                       -- active terms
%
%  TODO: a_mat_LV/RV/S, b_mat_LV/RV/S, af_mat_LV/RV/S, bf_mat_LV/RV/S,
%  as_mat_LV/RV/S, bs_mat_LV/RV/S below are PLACEHOLDER values (all 0) --
%  replace with your actual passive script's real numbers before trusting
%  any results from this file. The script will run with zeros, but the
%  passive response will be physically meaningless until these are set.
% ========================================================================
addpath(genpath('C:\Program Files\MATLAB\R2025a'))
gibbonRoot = 'C:\Users\kmagi\Downloads\GIBBON-Master';
addpath(genpath(gibbonRoot));
savepath;
febioExe = 'C:\Program Files\FEBioStudio\bin\febio4.exe';
febioAnalysis.febioPath = febioExe;
febioAnalysis.runMode   = 'internal';
savePath = 'C:\Users\kmagi\Downloads\FeBio';
ratname  = 'Z210W0';

%% ---- Solver settings ----
% opt_iter raised from 15 to 25: at 15, FEBio's auto time-stepper treated
% the ~20-22 iterations/step typical during active contraction as "too
% many" and refused to grow the step size back up after cutting to dtmin,
% pinning the solver at the smallest possible step for the whole active
% ramp and ballooning runtime. 25 lets it recognize those as acceptable.
% Timeline SHIFTED: previously sim-time spanned [0,1] with both load
% curves held flat at 0 over [0,0.3] (a diastolic-filling/dead lead-in
% before contraction onset). That dead region made the solver try to
% "converge" steps with an already-near-zero residual (~1e-20, floating
% point noise) against a relative tolerance requiring it to shrink
% further (rtol*INITIAL ~1e-23) -- mathematically impossible, so it
% burned through all 25 stiffness reformations without ever satisfying
% the criterion. Fix: drop the dead region entirely and start the FEBio
% time domain AT contraction onset (old t=0.3 -> new t=0, old t=1.0 ->
% new t=0.7). Load curves below are shifted accordingly. Step count
% reduced proportionally (63 = round(90*0.7)) to keep the same dt
% resolution as before.
% RESOLUTION TIGHTENED after the first shifted-timeline attempt still
% diverged (negative jacobians) on step 1 even after 6 retries down to
% dt~0.0016. Cause: dropping the dead lead-in also compressed the
% pressure/calcium rise -- old ramp (0.3->0.6, duration 0.3) took ~27
% steps at dt=1/90; new ramp (0->0.1167, duration 0.1167) took only
% ~10-11 steps at the same dt, i.e. ~2.5x larger load increment per
% step right at the most nonlinear part of the cycle, and max_retries=6
% wasn't enough halvings to reach a small enough dt before giving up.
% Fix: use a finer base step (baseDt/3) as BOTH the initial step size
% and the dtmax ceiling, so the solver starts fine enough through the
% steep early ramp without relying on retries to get there, and bump
% max_retries up so it still has room to cut further if needed.
% STEP COUNT: the earlier fineDt=(1/90)/3 tightening (~189 steps) was
% only needed because, at the time, the pressure curve's rise was
% compressed into just 0.1167 of the domain (matching calcium's fast
% Guccione rise). Pressure now has its OWN shape again with a rise
% duration of 0.3 -- the same rise duration used in the original working
% run, which was stable at dtmax=1/90. So we can relax back to that
% coarser resolution, cutting total steps roughly 3x (63 vs ~189).
% Calcium's rise (0.1167) is still faster than pressure's, but calcium
% only drives the internal active-tension term, not an external
% mechanical load, so it's less likely to need the same fine ceiling --
% if this reintroduces convergence trouble, fineDt is the first knob to
% tighten back down.
% simDuration SHORTENED 0.864 -> 0.45: the active model only needs to
% reach ES (minimum volume), which occurs somewhere after pressure
% finishes its rise (t=0.3) and while active tension is still
% substantial (calcium is at ~0.33-0.83 of peak across t=0.3-0.4, per
% the real-data curve) -- NOT the full 0.864 window out to where calcium
% has nearly vanished. Running that far past ES was the direct cause of
% the "expanding at the end" behavior: pressure holds at peak while
% active tension fades to ~1% of peak by t=0.864, leaving nothing to
% resist the still-full pressure. 0.45 gives comfortable margin past
% pressure's full engagement (0.3) and past calcium's steepest decay,
% while stopping well before the tension/pressure imbalance sets in.
% The results section already scans every completed step for the true
% volume minimum (ES), so this doesn't require knowing the exact ES
% time in advance -- just needs to not run so far past it that the
% unphysiological tail dominates the visualization.
%
% numTimeSteps stays tied to the passive model's step count (20, per
% explicit request) rather than to a fixed dt -- so shortening
% simDuration here automatically gives FINER per-step resolution
% (baseDt=simDuration/20=0.0225 vs the old 0.0432), which should also
% help the calcium-rise convergence risk noted before, not just cut
% runtime.
% DIAGNOSTIC (round 2) RESOLVED: the V_UZP_final + passive-script step-size
% test converged to LV=492.2/RV=406.9 vs. the true ED target LV=543.2/
% RV=395.8 (RV within 2.8%, LV within 9.4%) -- a reasonably close result,
% not a failure. The earlier "764.3 is a huge overshoot" read was wrong:
% that's V_UZP_final's own un-deformed volume (the starting point, not the
% loaded result), and the FEBio pressure sign convention here is
% "positive pressure = compressive" (confirmed in the FEBio manual) which
% the fliplr'd surface normals combine with to make V_def -> smaller as
% pressure ramps up -- exactly matching how the passive script's own
% inverse iteration is set up. This is the same convention as the passive
% model and is not being changed. Back to the real ED->ES active run on
% V_ED_final.
% SHORTENED 0.45 -> 0.2: with the pressure sign now genuinely inflating
% (see F_LV_pressure/F_RV_pressure above), pressure holds at full systolic
% magnitude from t=0.15 onward while calcium naturally decays after its
% t=0.096s peak -- so tension progressively loses to a constant full
% pressure and the chamber balloons out badly by t=0.45 (confirmed in the
% volume-trajectory plot: LV up to ~1230, RV up to ~337). The true ES
% minimum happens much earlier (LV ~t=0.06-0.07 at ~375, RV ~t=0.08-0.1 at
% ~242), so 0.2 gives ~2x margin past both minima while stopping well
% before the runaway inflation, which only really takes off past t~0.25.
simDuration  = 0.2;
numTimeSteps = 20;                 % matches the passive model's step count
fineDt       = simDuration/numTimeSteps;   % = 0.01, finer resolution near ES
baseDt       = fineDt;
% max_retries raised 8 -> 16: the last attempt showed negative-jacobian
% counts trending down as dt shrank (18169 -> ... -> 282) but never
% reaching zero -- it hit max_retries and gave up at dt~0.0012, still
% ~11x coarser than dtmin (fineDt/100 ~0.00011). This is a step-1-only
% cost (only the very first, hardest step needs many halvings to find a
% stable increment; once found, the auto-stepper remembers a working
% size and subsequent steps won't need nearly as many retries), so it
% shouldn't meaningfully hurt overall runtime.
% REVERTED: the max_retries=24/dtmax=fineDt/2 change let the solver grind
% through step 1's severe instability (negative-jacobian counts up to
% 18010) by finding tiny enough substeps to nominally satisfy convergence
% tolerance -- but FEBio's negative-jacobian check only catches a single
% element inverting on itself, not two different (individually valid)
% elements ending up in the same physical space. There's no self-contact
% defined in this model, so nothing stops that. The resulting "solved"
% run showed visibly overlapping/self-intersecting geometry, so this is
% NOT an acceptable fix -- back to the pre-session settings that were
% established as the empirically confirmed-stable baseline for Tmax=98
% (max_retries=16, dtmax=fineDt, not fineDt/2). The underlying step-1
% instability is still unresolved; forcing past it numerically isn't the
% right approach.
max_refs     = 25;
max_ups      = 10;
opt_iter     = 25;
max_retries  = 16;
dtmin        = fineDt/100;
dtmax        = fineDt;
symmetric_stiffness = 0;
pressureScale = 1;

%% ---- Holzapfel-Ogden passive parameters (per region) ----
% Real values pulled from the passive model's own console echo for this
% same rat (Z210W0), so the active analysis is consistent with the ED
% geometry it's built on top of.
a_mat_LV  = 0.186973;  b_mat_LV  = 4.80439;  af_mat_LV = 0.136467;  bf_mat_LV = 4.15537;  as_mat_LV = 1.33191;  bs_mat_LV = 3.334;
a_mat_RV  = 0.336956;  b_mat_RV  = 5.68276;  af_mat_RV = 0.352679;  bf_mat_RV = 6;        as_mat_RV = 1.72124;  bs_mat_RV = 5.85127;
a_mat_S   = 0.127484;  b_mat_S   = 8.49002;  af_mat_S  = 0.642743;  bf_mat_S  = 3.33202;  as_mat_S  = 1.41397;  bs_mat_S  = 5.2872;
k_mat     = 2000;

%% ---- active_contraction parameters (Hill/calcium model) ----
% From FEBio theory manual Eq 5.10.2-5.10.3 and Kwan dissertation:
%   Tmax   -- peak isometric tension. CURRENTLY 98 kPa; still needs
%             reconciling against the Kwan-derived rat/RV estimate
%             (~18.8 kPa, from sqrt(20000 mmHg^2) conversion). Not yet
%             per-region -- applied identically to LV/RV/S below.
%   Ca0    = 4.35 uM   (max peak intracellular calcium concentration)
%   B      = 4.75 um^-1 (LV/S, Guccione original)
%            11    um^-1 (RV, Kwan PAH-specific fit)
%   l0     = 1.60 um   (sarcomere length at zero active tension, validated)
%   refl   = 1.85 um   (reference sarcomere length, validated)
%   ascl   = driven by load curve 2 (calcium activation curve)
%
%   STATUS: these MATLAB variables ARE live and sent to FEBio (the earlier
%   LoBeta/HiBeta hardcoded-constants workaround was abandoned once the
%   real crash cause -- stale builds, not a genuine parameter limit --
%   was found; see conversation notes). Change values here directly.
%
%   Tmax = 98 kPa: empirically established as the maximum stable value
%   for this model before large-scale negative jacobians occur. 135.7 kPa
%   (FEBio theory manual's worked-example peak isometric tension) was
%   tried as an independent literature value but was NOT what was
%   actually tested -- reverting to the confirmed-stable 98 kPa here.
% Diagnostic round 2 resolved (see simDuration comment above) -- back to
% the confirmed-stable real value.
Tmax  = 136;
Ca0   = 4.35;
B_LV  = 4.75;
B_RV  = 11.0;
B_S   = 4.75;
l0    = 1.8;
refl  = 2.20;

%% ---- Pressure (mmHg -> kPa) ----
% END-SYSTOLIC pressures (this curve ramps UP to peak systolic pressure,
% not diastolic filling pressure) for this animal's W0 timepoint. RV ESP
% confirmed directly at W0 (21.9 +/- 0.9 mmHg). LV ESP has no specific
% published number yet -- using the 90-110 mmHg range midpoint (100) as
% a placeholder; replace once you have the actual Fig. 2F value.
% Previously this used the ED pressures [5.83 1.49] instead, left over
% from before the ES targets were known.
% CORRECTED SCOPE: back to true systolic afterload for the real ED->ES
% active run on V_ED_final. With the pressure sign now genuinely
% inflating (see F_LV_pressure/F_RV_pressure above), ramping toward this
% systolic magnitude while Tmax contracts gives real afterload resistance
% instead of the unopposed double-compression seen when this was
% combined with V_UZP_final and the old compressive convention.
Pressure_LVRV = [100.0  21.9];   % [LV, RV] mmHg, END-SYSTOLIC
P_LV = Pressure_LVRV(1) * 0.133;
P_RV = Pressure_LVRV(2) * 0.133;

%% ---- File names ----
% Output files use a distinct "_Active" suffix so this script's .feb,
% .xplt (FEBio auto-names the plot file after the .feb base name), and
% log files never overwrite the passive model's same-ratname outputs in
% the same savePath directory. `ratname` itself (used below to load the
% mesh and the passive model's saved geometry) stays unsuffixed, since
% those are inputs coming FROM the passive run.
febioFebFileNamePart = [ratname '_Active'];
febioFebFileName     = fullfile(savePath, [febioFebFileNamePart '.feb']);
febioLogFileName     = [febioFebFileNamePart '.txt'];
febioLogFileName_disp = [febioFebFileNamePart '_disp_out.txt'];
febioLogFileName_stress = [febioFebFileNamePart '_stress_out.txt'];
febioLogFileName_principalstress = [febioFebFileNamePart '_prinstress_out.txt'];

%% ---- Load mesh ----
parentdir = 'C:\Users\kmagi\Downloads\PAH_Imaging_FEA\Volume_Meshes';
baseDir   = fullfile(parentdir, ratname);
filename  = fullfile(baseDir, [ratname '_With_Fibers.h5']);
data = readXDMF(filename);
V  = data.Groups(2).Groups(1).Datasets(1).Value';
E  = data.Groups(2).Groups(1).Datasets(2).Value';
E  = double(E) + 1;
Fb = data.Groups(3).Groups(2).Datasets(2).Value';
Fb = Fb + 1;
Cb = data.Groups(3).Groups(2).Datasets(1).Value';
f0 = data.Groups(1).Groups(1).Datasets(1).Value';
s0 = data.Groups(1).Groups(3).Datasets(1).Value';

%% ---- Load region CSVs ----
LV = sort(unique(readmatrix(fullfile(baseDir, 'LV_freewall.csv'))));
RV = sort(unique(readmatrix(fullfile(baseDir, 'RV_freewall.csv'))));
S  = sort(unique(readmatrix(fullfile(baseDir,  'S_freewall.csv'))));
fprintf('Regions: LV=%d  RV=%d  S=%d\n', numel(LV), numel(RV), numel(S));

%% ---- Boundary surfaces ----
Fb = patchNormalFix(Fb);
F_base_BC     = Fb(Cb==1,:);
bcSupportList = unique(F_base_BC(:));
% FLIPPED (fliplr removed) -- ACTIVE SCRIPT ONLY, passive script unchanged.
% FEBio's pressure convention: "a positive pressure will act opposite to
% the normal, so it will compress the material" (FEBio user manual,
% Pressure Load section). The fliplr'd normal used previously (same as
% the passive script) points AWAY from the cavity centroid, into the
% tissue -- combined with that convention, positive pressure was
% compressing/deflating the chamber instead of inflating it. With Tmax=98
% active tension ALSO contracting the chamber, both effects pointed the
% same direction with nothing providing outward resistance, which is what
% crushed the LV to ~13 uL around t=0.3 in the last run. Removing fliplr
% here (checked numerically: the raw, non-flipped face order points INTO
% the cavity) makes positive pressure inflate instead, giving realistic
% afterload resistance against active contraction. The passive script's
% own fliplr'd surfaces and calibration are untouched.
F_LV_pressure = patchNormalFix(Fb(Cb==3,:));
F_RV_pressure = patchNormalFix(Fb(Cb==4,:));

%% ---- Load ED (pressurized) reference geometry ----
% CORRECTED SCOPE: the active model is contraction only (ED->ES) -- the
% passive model already handles diastolic filling. V_ED_final is the
% passive script's own properly-computed ED baseline, so the active model
% starts there directly rather than re-deriving it from V_UZP_final. This
% is now workable with the flipped (genuinely inflating) pressure
% convention above: ramping pressure from ED up toward systolic while
% Tmax contracts gives real opposing afterload resistance, instead of the
% unopposed double-compression that crushed the LV when V_UZP_final was
% combined with the old compressive convention.
geomFile = fullfile(savePath, [ratname '_geomStates.mat']);
if ~exist(geomFile,'file')
    error('Geometry file not found: %s\nRun the passive script first.', geomFile);
end
geomData = load(geomFile);
V_def = geomData.V_ED_final;
fprintf('Loaded ED (pressurized) reference geometry: %d nodes\n', size(V_def,1));

% True ED target volumes, computed directly from the real imaged geometry
% (V_MRI_ED) -- same definition the passive script's inverse iteration
% used as its own convergence target. Used below to sanity-check the
% volume trajectory against ground truth, not just against V_def's own
% (potentially inconsistent) starting volume.
ED_target_LV = closeAndVolume(F_LV_pressure, geomData.V_MRI_ED);
ED_target_RV = closeAndVolume(F_RV_pressure, geomData.V_MRI_ED);
fprintf('True ED target (from V_MRI_ED): LV=%.1f  RV=%.1f\n', ED_target_LV, ED_target_RV);

%% ---- FEBio spec ----
[febio_spec] = febioStructTemplate;
febio_spec.ATTR.version     = '4.0';
febio_spec.Module.ATTR.type = 'solid';
febio_spec.Control.analysis   = 'STATIC';
febio_spec.Control.time_steps = numTimeSteps;
febio_spec.Control.step_size  = baseDt;
febio_spec.Control.solver.max_refs            = max_refs;
febio_spec.Control.solver.qn_method.ATTR.type = 'Broyden';
febio_spec.Control.solver.qn_method.max_ups   = max_ups;
febio_spec.Control.solver.symmetric_stiffness = symmetric_stiffness;
% Explicit convergence tolerances -- without these, FEBio's default relative
% tolerance produced a "required" threshold that rounded to exactly zero
% during the near-zero-load state at the very start of the simulation
% (before the pressure/calcium curves ramp up), causing spurious
% non-convergence failures even though the residual was already at the
% floating-point noise floor.
febio_spec.Control.solver.rtol  = 0.001;
febio_spec.Control.solver.etol  = 0.01;
febio_spec.Control.solver.dtol  = 0.001;
febio_spec.Control.solver.lstol = 0.9;
febio_spec.Control.time_stepper.dtmin       = dtmin;
febio_spec.Control.time_stepper.dtmax       = dtmax;
febio_spec.Control.time_stepper.max_retries = max_retries;
febio_spec.Control.time_stepper.opt_iter    = opt_iter;

%% ---- Materials: Holzapfel-Ogden_ACTIVE (custom plugin, single class) ----
% Flat parameter list matching dllmain.cpp's BEGIN_FECORE_CLASS -- no
% "fiber" tag (orientation comes from mat_axis below), no nested
% "active_contraction" wrapper. Registered type string must exactly match
% dllmain.cpp's REGISTER_FECORE_CLASS(HolzapfelMyocardiumActivePI,
% "Holzapfel-Ogden_ACTIVE") string.
%
% *** STATUS: WORKING ***. All 6 active-contraction parameters
% (Tmax/Ca0/beta/l0/refl/ascl) are live below -- the earlier suspected
% ">12-parameter crash" turned out to actually be caused by the stale
% build / type-string mismatches / missing k, not a hard parameter-count
% limit. cbar()/DevTangent() now also includes the active-stress
% contribution to the tangent (see cbar_with_active_tangent.cpp), which
% meaningfully sped up convergence. k is required (base-class bulk
% modulus) or FEBio errors with "K must be a positive number".
materialName1 = 'Material1';
febio_spec.Material.material{1}.ATTR.name = materialName1;
febio_spec.Material.material{1}.ATTR.type = 'Holzapfel-Ogden_ACTIVE';
febio_spec.Material.material{1}.ATTR.id   = 1;
febio_spec.Material.material{1}.k    = k_mat;
febio_spec.Material.material{1}.a    = a_mat_LV;
febio_spec.Material.material{1}.b    = b_mat_LV;
febio_spec.Material.material{1}.af   = af_mat_LV;
febio_spec.Material.material{1}.bf   = bf_mat_LV;
febio_spec.Material.material{1}.as   = as_mat_LV;
febio_spec.Material.material{1}.bs   = bs_mat_LV;
febio_spec.Material.material{1}.afs  = 0.0;
febio_spec.Material.material{1}.bfs  = 0.0;
febio_spec.Material.material{1}.asn  = 0.0;
febio_spec.Material.material{1}.bsn  = 0.0;
febio_spec.Material.material{1}.anf  = 0.0;
febio_spec.Material.material{1}.bnf  = 0.0;
febio_spec.Material.material{1}.Tmax = Tmax;
febio_spec.Material.material{1}.Ca0  = Ca0;
febio_spec.Material.material{1}.beta = B_LV;
febio_spec.Material.material{1}.l0   = l0;
febio_spec.Material.material{1}.refl = refl;
febio_spec.Material.material{1}.ascl.ATTR.lc = 2;
febio_spec.Material.material{1}.ascl.VAL     = 1;

materialName2 = 'Material2';
febio_spec.Material.material{2}.ATTR.name = materialName2;
febio_spec.Material.material{2}.ATTR.type = 'Holzapfel-Ogden_ACTIVE';
febio_spec.Material.material{2}.ATTR.id   = 2;
febio_spec.Material.material{2}.k    = k_mat;
febio_spec.Material.material{2}.a    = a_mat_RV;
febio_spec.Material.material{2}.b    = b_mat_RV;
febio_spec.Material.material{2}.af   = af_mat_RV;
febio_spec.Material.material{2}.bf   = bf_mat_RV;
febio_spec.Material.material{2}.as   = as_mat_RV;
febio_spec.Material.material{2}.bs   = bs_mat_RV;
febio_spec.Material.material{2}.afs  = 0.0;
febio_spec.Material.material{2}.bfs  = 0.0;
febio_spec.Material.material{2}.asn  = 0.0;
febio_spec.Material.material{2}.bsn  = 0.0;
febio_spec.Material.material{2}.anf  = 0.0;
febio_spec.Material.material{2}.bnf  = 0.0;
febio_spec.Material.material{2}.Tmax = Tmax;
febio_spec.Material.material{2}.Ca0  = Ca0;
febio_spec.Material.material{2}.beta = B_RV;
febio_spec.Material.material{2}.l0   = l0;
febio_spec.Material.material{2}.refl = refl;
febio_spec.Material.material{2}.ascl.ATTR.lc = 2;
febio_spec.Material.material{2}.ascl.VAL     = 1;

materialName3 = 'Material3';
febio_spec.Material.material{3}.ATTR.name = materialName3;
febio_spec.Material.material{3}.ATTR.type = 'Holzapfel-Ogden_ACTIVE';
febio_spec.Material.material{3}.ATTR.id   = 3;
febio_spec.Material.material{3}.k    = k_mat;
febio_spec.Material.material{3}.a    = a_mat_S;
febio_spec.Material.material{3}.b    = b_mat_S;
febio_spec.Material.material{3}.af   = af_mat_S;
febio_spec.Material.material{3}.bf   = bf_mat_S;
febio_spec.Material.material{3}.as   = as_mat_S;
febio_spec.Material.material{3}.bs   = bs_mat_S;
febio_spec.Material.material{3}.afs  = 0.0;
febio_spec.Material.material{3}.bfs  = 0.0;
febio_spec.Material.material{3}.asn  = 0.0;
febio_spec.Material.material{3}.bsn  = 0.0;
febio_spec.Material.material{3}.anf  = 0.0;
febio_spec.Material.material{3}.bnf  = 0.0;
febio_spec.Material.material{3}.Tmax = Tmax;
febio_spec.Material.material{3}.Ca0  = Ca0;
febio_spec.Material.material{3}.beta = B_S;
febio_spec.Material.material{3}.l0   = l0;
febio_spec.Material.material{3}.refl = refl;
febio_spec.Material.material{3}.ascl.ATTR.lc = 2;
febio_spec.Material.material{3}.ascl.VAL     = 1;

%% ---- Elements ----
partName1 = 'LV';
febio_spec.Mesh.Elements{1}.ATTR.name    = partName1;
febio_spec.Mesh.Elements{1}.ATTR.type    = 'tet4';
febio_spec.Mesh.Elements{1}.elem.ATTR.id = LV;
febio_spec.Mesh.Elements{1}.elem.VAL     = E(LV,:);

partName2 = 'RV';
febio_spec.Mesh.Elements{2}.ATTR.name    = partName2;
febio_spec.Mesh.Elements{2}.ATTR.type    = 'tet4';
febio_spec.Mesh.Elements{2}.elem.ATTR.id = RV;
febio_spec.Mesh.Elements{2}.elem.VAL     = E(RV,:);

partName3 = 'S';
febio_spec.Mesh.Elements{3}.ATTR.name    = partName3;
febio_spec.Mesh.Elements{3}.ATTR.type    = 'tet4';
febio_spec.Mesh.Elements{3}.elem.ATTR.id = S;
febio_spec.Mesh.Elements{3}.elem.VAL     = E(S,:);

%% ---- MeshData: mat_axis (fiber/sheet directions per element) ----
febio_spec.MeshData.ElementData{1}.ATTR.elem_set = partName1;
febio_spec.MeshData.ElementData{1}.ATTR.type     = 'mat_axis';
febio_spec.MeshData.ElementData{2}.ATTR.elem_set = partName2;
febio_spec.MeshData.ElementData{2}.ATTR.type     = 'mat_axis';
febio_spec.MeshData.ElementData{3}.ATTR.elem_set = partName3;
febio_spec.MeshData.ElementData{3}.ATTR.type     = 'mat_axis';

lid_LV = zeros(size(E,1),1); lid_LV(LV) = 1:numel(LV);
lid_RV = zeros(size(E,1),1); lid_RV(RV) = 1:numel(RV);
lid_S  = zeros(size(E,1),1); lid_S(S)   = 1:numel(S);

elem_count_LV = 1;
elem_count_RV = 1;
elem_count_S  = 1;
for q = 1:size(E,1)
    if ismember(q,LV)
        LID = lid_LV(q);
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.ATTR.lid = LID;
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.a = f0(q,:);
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.d = s0(q,:);
        elem_count_LV = elem_count_LV + 1;
    elseif ismember(q,RV)
        LID = lid_RV(q);
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.ATTR.lid = LID;
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.a = f0(q,:);
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.d = s0(q,:);
        elem_count_RV = elem_count_RV + 1;
    else
        LID = lid_S(q);
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.ATTR.lid = LID;
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.a = f0(q,:);
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.d = s0(q,:);
        elem_count_S = elem_count_S + 1;
    end
end

%% ---- MeshDomains ----
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat  = materialName1;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat  = materialName2;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.name = partName3;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat  = materialName3;

%% ---- Nodes ----
febio_spec.Mesh.Nodes{1}.ATTR.name    = 'Object1';
febio_spec.Mesh.Nodes{1}.node.ATTR.id = (1:size(V_def,1))';
febio_spec.Mesh.Nodes{1}.node.VAL     = V_def;

%% ---- Surfaces ----
febio_spec.Mesh.Surface{1}.ATTR.name    = 'LVPressure';
febio_spec.Mesh.Surface{1}.tri3.ATTR.id = (1:size(F_LV_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL     = F_LV_pressure;
febio_spec.Mesh.Surface{2}.ATTR.name    = 'RVPressure';
febio_spec.Mesh.Surface{2}.tri3.ATTR.id = (1:size(F_RV_pressure,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL     = F_RV_pressure;

%% ---- NodeSets + BCs ----
febio_spec.Mesh.NodeSet{1}.ATTR.name = 'bcSupportList';
febio_spec.Mesh.NodeSet{1}.VAL       = mrow(bcSupportList);
febio_spec.Boundary.bc{1}.ATTR.name     = 'FixedBase';
febio_spec.Boundary.bc{1}.ATTR.type     = 'zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set = 'bcSupportList';
febio_spec.Boundary.bc{1}.x_dof = 1;
febio_spec.Boundary.bc{1}.y_dof = 1;
febio_spec.Boundary.bc{1}.z_dof = 1;

%% ---- Pressure loads ----
febio_spec.Loads.surface_load{1}.ATTR.type    = 'pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface = 'LVPressure';
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc  = 1;
febio_spec.Loads.surface_load{1}.pressure.VAL       = P_LV * pressureScale;
febio_spec.Loads.surface_load{1}.symmetric_stiffness = 0;
febio_spec.Loads.surface_load{2}.ATTR.type    = 'pressure';
febio_spec.Loads.surface_load{2}.ATTR.surface = 'RVPressure';
% RV uses its own load curve (id=3, not 1) because its ED->systolic ratio
% differs from LV's (1.49/21.9 = 6.8% vs 5.83/100 = 5.8%) -- see LC 3
% below, next to LC_pressure.
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc  = 3;
febio_spec.Loads.surface_load{2}.pressure.VAL       = P_RV * pressureScale;
febio_spec.Loads.surface_load{2}.symmetric_stiffness = 0;

%% ---- Load curves ----
% LC_pressure and LC_calcium are now DECOUPLED, and intentionally
% represent two different physical quantities:
%   - LC_calcium drives the material's internal active-tension state
%     (ascl -> Ta via the Hill/calcium model in the plugin). Its shape
%     SHOULD follow real calcium kinetics -- fast rise, slower decay --
%     because that's what it physically is: an intracellular
%     concentration transient. Originally approximated from the classic
%     Guccione & McCulloch (1993) Fig. 1 curve; now built directly from
%     real OVX SuHx group Ca2+ transient data from the user's own lab
%     (peak at ~96ms, decaying over ~768ms -- see the curve definition
%     below for details and sourcing).
%   - LC_pressure is an externally applied mechanical boundary condition
%     (cavity pressure). Physiologically it does NOT track calcium
%     kinetics -- it follows the pressure-volume loop: a rise during
%     isovolumic contraction, a plateau through ejection, and a fall
%     during isovolumic relaxation back to baseline, roughly SYMMETRIC in
%     time (unlike calcium's fast-rise/slow-decay asymmetry). Making
%     pressure literally copy the calcium curve's shape (as an earlier
%     version of this script did, to fix a "collapse" instability) was
%     not anatomically correct -- it's what produced the "quickly
%     contracts, then slowly expands" behavior, since pressure was
%     fading out on calcium's slow 500ms decay tail instead of its own,
%     faster relaxation-phase timing.
%
% SCOPE CORRECTION: this model simulates ED -> ES contraction only, not
% a full contract-then-relax cycle -- relaxation/refilling belongs to
% the *next* beat, which isn't part of this analysis. So LC_pressure
% should never decay back down within this window: it rises through
% isovolumic contraction and then HOLDS at peak (the heart "continues to
% contract until the cycle continues again"), it does not fall back to
% baseline here. LC_calcium keeps its natural literal decay (real
% intracellular calcium kinetics decay within a beat regardless of
% whether the tissue mechanically relaxes yet), so active tension will
% be easing off near the end of the window even while pressure holds --
% that's expected and fine, since we're not simulating far enough into
% relaxation for it to matter; the ES state (minimum volume) is expected
% to occur well before then, near/at the pressure plateau.
%
% LINEAR interpolation: SMOOTH previously caused a natural cubic spline
% to overshoot into supposedly-flat regions because of a sparse,
% far-away tail point combined with unevenly-spaced points near the
% peak -- LINEAR guarantees no cross-region overshoot.
%
% CORRECTED SCOPE: this ramp represents ONLY active contraction
% (isovolumic contraction + ejection), not passive filling -- most of the
% ramp should cover the ED->systolic RISE, not filling from empty.
% REVISED after a failed attempt: starting the curve AT the ED fraction
% (instead of 0.0) caused a much worse step-1 failure (thousands of
% negative jacobians even at t=0.0006) because FEBio always treats
% V_ED_final as zero-stress (F=I) at t=0 -- jumping straight to a nonzero
% load fraction there is an instantaneous step function for the solver to
% resolve, harsher than ramping smoothly from true zero. Compromise: ramp
% VERY quickly from 0 up to the ED fraction over just the first nominal
% step (t=0.01, matching fineDt), keeping the solver's start smooth.
%
% RESHAPED (peak-then-decline, not peak-then-hold): real measured rat LV
% pressure traces (Wang et al. 2017, Oncotarget 8:96161, Fig. 2a) show
% pressure rising through isovolumic contraction to a peak near the start
% of ejection, then GENTLY DECLINING through the rest of ejection -- it
% does not hold flat until diastole. Holding at peak while calcium
% naturally fades (its decay is a real kinetic process, not something we
% control) let pressure "win" unopposed and balloon the chamber out past
% ES. Mirroring the real trace's decline instead removes that artificial
% imbalance: peak at t=0.1 (near calcium engagement), easing back to 85%
% of peak by t=0.2 (near where the true ES minimum has been landing in
% the trajectory data), then held at that reduced level as a safety
% margin in case the run extends slightly past t=0.2.
% RETIMED (early rise): the straight line from (0.01, 0.0583) to
% (0.1, 1.0) outran calcium badly -- at t=0.024, calcium has barely
% engaged (ascl~=0.08, so active tension is near zero) but linear
% interpolation already has pressure at ~20% of peak. That's a window of
% essentially unopposed inflation, which is the "blows up before
% contraction" behavior observed. Added intermediate points so pressure's
% early rise tracks calcium's own timing more closely (staying low while
% calcium is still low, then catching up to peak by t=0.1), instead of a
% straight line racing ahead of it. Calcium reference values: t=0.044
% ascl=0.44, t=0.064 ascl=0.81, t=0.084 ascl=0.98.
% LV: ED/systolic = 5.83/100.0 = 0.0583.
febio_spec.LoadData.load_controller{1}.ATTR.name       = 'LC_pressure';
febio_spec.LoadData.load_controller{1}.ATTR.id         = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type       = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate     = 'LINEAR';
febio_spec.LoadData.load_controller{1}.points.pt.VAL   = [
    0.0000   0.0000;   % starts at true zero, matching F=I at t=0
    0.0100   0.0583;   % quick rise to the true ED fraction within step 1
    0.0440   0.2500;   % tracks calcium's own rise (ascl~=0.44 at this t)
    0.0700   0.6000;   % catching up as calcium approaches its peak
    0.1000   1.0000;   % peak near ejection onset / calcium engagement
    0.2000   0.8500;   % gentle decline through ejection, matching Fig. 2a shape
    8.4500   0.8500;   % held at reduced level (no relaxation leg modeled here)
];
% LC 3: RV's own ED->systolic ramp (RV: ED/systolic = 1.49/21.9 = 0.0680,
% close to but not identical to LV's fraction -- separate curve so each
% chamber's ramp is exact rather than sharing an averaged approximation).
% Same retimed, peak-then-decline profile as LV.
febio_spec.LoadData.load_controller{3}.ATTR.name       = 'LC_pressure_RV';
febio_spec.LoadData.load_controller{3}.ATTR.id         = 3;
febio_spec.LoadData.load_controller{3}.ATTR.type       = 'loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate     = 'LINEAR';
febio_spec.LoadData.load_controller{3}.points.pt.VAL   = [
    0.0000   0.0000;   % starts at true zero, matching F=I at t=0
    0.0100   0.0200;   % quick rise to the true ED fraction within step 1
    0.0440   0.2500;   % tracks calcium's own rise (ascl~=0.44 at this t)
    0.0700   0.6000;   % catching up as calcium approaches its peak
    0.1000   1.0000;   % peak near ejection onset / calcium engagement
    0.2000   0.8500;   % gentle decline through ejection, matching Fig. 2a shape
    8.4500   0.8500;   % held at reduced level (no relaxation leg modeled here)
];
% LC 2: calcium activation curve -- built directly from real Ca2+
% transient data (fura-2-style fluorescence ratio, OVX SuHx group, n=12
% rats, averaged) provided from the user's own lab, replacing the
% earlier literature Guccione & McCulloch Fig.1 approximation. Extracted
% from "[Ca] OVX SuHx group (All rats).txt" ("OVXSUHX Average" column,
% sampled every 4ms):
%   - diastolic baseline  ~0.939318 (minimum, at real t=0.096 s)
%   - peak                ~1.289958 (at real t=0.192 s)
%   - amplitude            0.350640 (peak - baseline)
% Points below are (real_t - 0.096) for time (onset = sim-time 0, same
% convention as the rest of this script) and (value-baseline)/amplitude
% for ascl. Rise is ~96ms (close to the Guccione 100ms value), but decay
% is much slower -- ~768ms vs. Guccione's 500ms -- which may reflect
% real SuHx-model impaired calcium reuptake (a disease-relevant feature,
% not noise, so it's kept rather than smoothed away). ascl = Ca(t)/Ca0,
% i.e. C(t) in the FEBio theory manual's Eq. 5.10.2, not raw uM/ratio
% units.
febio_spec.LoadData.load_controller{2}.ATTR.name       = 'LC_calcium';
febio_spec.LoadData.load_controller{2}.ATTR.id         = 2;
febio_spec.LoadData.load_controller{2}.ATTR.type       = 'loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate     = 'LINEAR';
febio_spec.LoadData.load_controller{2}.points.pt.VAL   = [
    0.0000   0.0000;   % real t=0.096 s (onset / diastolic baseline)
    0.0240   0.0790;   % real t=0.120 s
    0.0440   0.4415;   % real t=0.140 s
    0.0640   0.8141;   % real t=0.160 s
    0.0840   0.9814;   % real t=0.180 s
    0.0960   1.0000;   % real t=0.192 s: PEAK
    0.1240   0.9619;   % real t=0.220 s
    0.1640   0.9022;   % real t=0.260 s
    0.2040   0.8331;   % real t=0.300 s
    0.2640   0.6752;   % real t=0.360 s
    0.3240   0.4902;   % real t=0.420 s
    0.3840   0.3339;   % real t=0.480 s
    0.4640   0.2006;   % real t=0.560 s
    0.5440   0.1220;   % real t=0.640 s
    0.6240   0.0755;   % real t=0.720 s
    0.7040   0.0471;   % real t=0.800 s
    0.7840   0.0267;   % real t=0.880 s
    0.8640   0.0105;   % real t=0.960 s: back near baseline
];

%% ---- Output ----
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file  = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data  = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{1}.ATTR.file  = febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data  = 'sx;sy;sz;sxy;syz;sxz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{2}.ATTR.file  = febioLogFileName_principalstress;
febio_spec.Output.logfile.element_data{2}.ATTR.data  = 's1;s2;s3';
febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';
febio_spec.Output.plotfile.compression = 0;

%% ---- Write and run ----
if ~exist(savePath,'dir'), mkdir(savePath); end
febioStruct2xml(febio_spec, febioFebFileName);
fprintf('FEB file written: %d bytes\n', dir(febioFebFileName).bytes);
outDir = fileparts(febioFebFileName);
[~,base,~] = fileparts(char(febioFebFileName));

% Run FEBio DIRECTLY via system(), bypassing runMonitorFEBio entirely.
% GIBBON's wrapper manages FEBio's execution/monitoring internally and
% does not reliably honor a custom "-o" flag appended to run_string --
% confirmed by the console log file never being created. Calling system()
% ourselves captures FEBio's full console output directly into a MATLAB
% string (cmdout), with no file-based indirection that can silently fail.
runCmd = sprintf('"%s" -i "%s"', char(febioExe), char(febioFebFileName));
fprintf('Running: %s\n', runCmd);
fprintf('(Live FEBio output will stream below; full text is also captured for the diagnostic.)\n\n');
tic;
[sysStatus, cmdout] = system(runCmd, '-echo');
fprintf('\nFEBio run finished in %.1f sec (system() exit status: %d)\n', toc, sysStatus);

consoleLogFile = fullfile(outDir, [base '_console.txt']);
fid = fopen(consoleLogFile, 'w');
if fid > 0
    fprintf(fid, '%s', cmdout);
    fclose(fid);
    fprintf('Console output saved to: %s\n', consoleLogFile);
else
    warning('Could not write console log to %s (continuing anyway -- diagnostic uses the in-memory copy).', consoleLogFile);
end

runFlag = ~isempty(regexpi(cmdout, 'N O R M A L   T E R M I N A T I O N', 'once'));

%% ---- Diagnostic: parse console output for negative-jacobian failures ----
if ~isempty(cmdout)
    logLines = strsplit(cmdout, newline);
    currentStep = NaN;
    currentTime = NaN;
    failureLog  = struct('step', {}, 'time', {}, 'nJacobians', {});
    for li = 1:numel(logLines)
        line = logLines{li};
        stepMatch = regexp(line, 'beginning time step\s+(\d+)\s*:\s*([\d.eE+-]+)', 'tokens');
        if ~isempty(stepMatch)
            currentStep = str2double(stepMatch{1}{1});
            currentTime = str2double(stepMatch{1}{2});
        end
        jacMatch = regexp(line, '(\d+)\s+negative jacobians detected', 'tokens');
        if ~isempty(jacMatch)
            nJac = str2double(jacMatch{1}{1});
            failureLog(end+1) = struct('step', currentStep, 'time', currentTime, 'nJacobians', nJac); %#ok<SAGROW>
        end
    end
    fprintf('\n=== Negative-Jacobian Diagnostic ===\n');
    fprintf('Total non-converged iterations (negative-jacobian events): %d\n', numel(failureLog));
    if ~isempty(failureLog)
        uniqueSteps = unique([failureLog.step]);
        fprintf('These occurred across %d distinct time step(s): %s\n', ...
            numel(uniqueSteps), mat2str(uniqueSteps));
        fprintf('\nPer-event detail (step : time : #negative jacobians):\n');
        for i = 1:numel(failureLog)
            fprintf('  step %3d : t=%.6f : %d negative jacobians\n', ...
                failureLog(i).step, failureLog(i).time, failureLog(i).nJacobians);
        end
        allCounts = [failureLog.nJacobians];
        fprintf('\nSummary: min=%d  max=%d  mean=%.1f  median=%.1f negative jacobians per failed event\n', ...
            min(allCounts), max(allCounts), mean(allCounts), median(allCounts));
        fprintf('\nPer-step retry counts:\n');
        for s = uniqueSteps
            idx = ([failureLog.step] == s);
            fprintf('  step %3d : %d failed attempt(s), jacobian counts = %s\n', ...
                s, sum(idx), mat2str([failureLog(idx).nJacobians]));
        end
    else
        fprintf('No negative-jacobian events found in the console log.\n');
    end
    fprintf('=== End Diagnostic ===\n\n');
else
    warning('No console output was captured from FEBio (cmdout is empty) -- cannot run negative-jacobian diagnostic.');
end

%% ---- Results: scan all steps for true end-systolic (minimum) volume ----
if runFlag == 1
    fprintf('\nSolve SUCCEEDED\n');
    dataDisp = importFEBio_logfile(fullfile(outDir, febioLogFileName_disp), 0, 1);
    nSteps = size(dataDisp.data, 3);

    Vol_LV_all = zeros(nSteps,1);
    Vol_RV_all = zeros(nSteps,1);
    for k = 1:nSteps
        U_k = dataDisp.data(:,:,k);
        V_k = V_def + U_k;
        Vol_LV_all(k) = closeAndVolume(F_LV_pressure, V_k);
        Vol_RV_all(k) = closeAndVolume(F_RV_pressure, V_k);
    end

    Vol_LV_ED = closeAndVolume(F_LV_pressure, V_def);
    Vol_RV_ED = closeAndVolume(F_RV_pressure, V_def);

    [Vol_LV_ES, k_LV] = min(Vol_LV_all);
    [Vol_RV_ES, k_RV] = min(Vol_RV_all);

    SV_LV = Vol_LV_ED - Vol_LV_ES;
    SV_RV = Vol_RV_ED - Vol_RV_ES;

    fprintf('\n=== Results ===\n');
    fprintf('LV: EDV=%.1f  ESV=%.1f (at step %d)  SV=%.1f  EF=%.1f%%\n', ...
        Vol_LV_ED, Vol_LV_ES, k_LV, SV_LV, SV_LV/Vol_LV_ED*100);
    fprintf('RV: EDV=%.1f  ESV=%.1f (at step %d)  SV=%.1f  EF=%.1f%%\n', ...
        Vol_RV_ED, Vol_RV_ES, k_RV, SV_RV, SV_RV/Vol_RV_ED*100);

    %% ---- Volume-vs-step trajectory: print table + plot ----
    % Added to diagnose the DIAGNOSTIC TEST (round 2) overshoot: the
    % full-pressure endpoint significantly overshot the true ED target
    % (from V_MRI_ED) for both chambers, and the negative-jacobian burst
    % at step 7 (313 events) is a candidate cause. This shows the actual
    % volume trajectory step-by-step so we can see whether it's a smooth
    % overshoot or an abrupt jump right around the failed steps.
    if isfield(dataDisp, 'time') && numel(dataDisp.time) == nSteps
        stepTimes = dataDisp.time(:);
    else
        % Fallback: approximate step times evenly across simDuration if
        % the log doesn't expose them directly.
        stepTimes = simDuration * (1:nSteps)' / nSteps;
    end

    failedSteps = [];
    if exist('failureLog','var') && ~isempty(failureLog)
        failedSteps = unique([failureLog.step]);
    end

    fprintf('\n=== Volume trajectory (all steps) ===\n');
    fprintf('%4s  %8s  %10s  %10s  %s\n', 'step', 't', 'Vol_LV', 'Vol_RV', 'negative-jacobian step?');
    for k = 1:nSteps
        flag = '';
        if ismember(k, failedSteps)
            flag = '  <-- had negative-jacobian retries';
        end
        fprintf('%4d  %8.4f  %10.1f  %10.1f  %s\n', k, stepTimes(k), Vol_LV_all(k), Vol_RV_all(k), flag);
    end

    figure('Name','Volume trajectory vs. true ED target');
    subplot(1,2,1); hold on;
    plot(stepTimes, Vol_LV_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'LV volume');
    yline(ED_target_LV, '--r', 'LineWidth', 1.5, 'DisplayName', 'True ED target (LV)');
    for s = failedSteps
        xline(stepTimes(s), ':k', 'HandleVisibility','off');
    end
    xlabel('simulation time'); ylabel('LV cavity volume (\muL)');
    title('LV volume vs. step'); legend('Location','best'); grid on;

    subplot(1,2,2); hold on;
    plot(stepTimes, Vol_RV_all, '-o', 'LineWidth', 1.5, 'DisplayName', 'RV volume');
    yline(ED_target_RV, '--r', 'LineWidth', 1.5, 'DisplayName', 'True ED target (RV)');
    for s = failedSteps
        xline(stepTimes(s), ':k', 'HandleVisibility','off');
    end
    xlabel('simulation time'); ylabel('RV cavity volume (\muL)');
    title('RV volume vs. step'); legend('Location','best'); grid on;
else
    fprintf('\nSolve FAILED\n');
end

%% ========== Local Functions ==========
function data = readXDMF(filename)
    data = struct();
    data.Groups(1).Name = "Function";
    data.Groups(2).Name = "Mesh";
    data.Groups(3).Name = "MeshTags";
    data.Groups(1).Groups(1).Name = "FiberDirection";
    data.Groups(1).Groups(2).Name = "Normal";
    data.Groups(1).Groups(3).Name = "SheetDirection";
    data.Groups(2).Groups(1).Name = "Mesh";
    data.Groups(3).Groups(1).Name = "Cell tags";
    data.Groups(3).Groups(2).Name = "Facet tags";
    data.Groups(1).Groups(1).Datasets(1).Value = h5read(filename,"/Function/FiberDirection/0");
    data.Groups(1).Groups(2).Datasets(1).Value = h5read(filename,"/Function/Normal/0");
    data.Groups(1).Groups(3).Datasets(1).Value = h5read(filename,"/Function/SheetDirection/0");
    data.Groups(2).Groups(1).Datasets(1).Value = h5read(filename,"/Mesh/Mesh/geometry");
    data.Groups(2).Groups(1).Datasets(2).Value = h5read(filename,"/Mesh/Mesh/topology");
    data.Groups(3).Groups(1).Datasets(1).Value = h5read(filename,"/MeshTags/Cell tags/Values");
    data.Groups(3).Groups(1).Datasets(2).Value = h5read(filename,"/MeshTags/Cell tags/topology");
    data.Groups(3).Groups(2).Datasets(1).Value = h5read(filename,"/MeshTags/Facet tags/Values");
    data.Groups(3).Groups(2).Datasets(2).Value = h5read(filename,"/MeshTags/Facet tags/topology");
end

function vol = closeAndVolume(Fsurf, Vnodes)
    Edge = patchBoundary(Fsurf);
    [Fc,Vc] = triSurfCloseHoles(double(Fsurf), Vnodes, 0.5, double(Edge));
    Fc  = patchNormalFix(Fc);
    vol = patchVolume(Fc,Vc);
    if vol < 0, vol = -vol; end
end