clear all
clc
close all
%% ========================================================================
%  ACTIVE CONTRACTION MODEL -- BUILT-IN TRANSVERSE-ISOTROPIC ACTIVE MATERIAL
%
%  FIXED VERSION:
%   1) Bulk modulus k_mat moved from the nested Holzapfel_Ogden solid onto
%      the outer 'uncoupled solid mixture' material for all three regions.
%      Per the Holzapfel_Ogden plugin's own manual (already noted in
%      FiberStressTest.m): "as any other uncoupled material, the
%      Holzapfel-Ogden model can be included within an uncoupled solid
%      mixture. In that case, the bulk modulus should be outside any solid
%      domain." Leaving k nested inside the passive solid, as the original
%      IsotropicTest.m did, is inconsistent with that guidance and is a
%      likely source of incorrect near-incompressibility behavior, which
%      directly affects cavity volumes (EDV/ESV).
%   2) Added a calibration diagnostic block after the volume-trajectory
%      loop: prints pressure fraction, calcium fraction, and both cavity
%      volumes per step, checks step 1 against the analytic ED target, and
%      checks whether the step of minimum volume (simulated ES) lines up
%      with the step where pressure reaches its target level.
% ========================================================================
addpath(genpath('C:\Program Files\MATLAB\R2025a'))
gibbonRoot = 'C:\Users\shreygandhi\Desktop\GIBBON';
addpath(genpath(gibbonRoot));
savepath;
febioExe = 'C:\Program Files\FEBioStudio\bin\febio4.exe';
febioAnalysis.febioPath = febioExe;
febioAnalysis.runMode   = 'internal';
savePath = 'C:\Users\shreygandhi\Desktop\ActiveModel\FeBio';
ratname  = 'Z210W0';

%% ---- Solver settings ----
% Same settings established as stable for the custom-plugin active run.
% startMode:
%   'ED'  = start from the passive script's saved ED_sim geometry (V_ED_final).
%           Needs NO change to the passive script. ED is already reached at
%           t=0, so only the pressure ABOVE end-diastolic pressure is applied.
%   'UZP' = two-phase run from the saved stress-free geometry (V_UZP_final):
%           fill to EDP, then systole. Only valid if V_UZP_final is smaller
%           than the ED target (check the printed 'Start volume' line).
startMode = 'ED';
if strcmp(startMode,'ED')
    fastMode     = true;   % true = about 3 to 4x faster; false = v4 settings (slow, finer)
    if fastMode
        simDuration  = 0.17;   % pressure is back near EDP by t = 0.16, nothing useful after
        numTimeSteps = 34;     % dt = 0.005
    else
        simDuration  = 0.2;
        numTimeSteps = 80;     % dt = 0.0025
    end
    tED          = 0;      % ED is the starting state
else
    simDuration  = 0.25;
    numTimeSteps = 50;     % dt = 0.005, t = 0.05 lands on a step
    tED          = 0.05;   % end of passive-fill phase (ED)
end
fineDt       = simDuration/numTimeSteps;
baseDt       = fineDt;
max_refs     = 25;
max_ups      = 10;
opt_iter     = 25;
if exist('fastMode','var') && fastMode
    max_retries = 6;           % fail fast instead of burning an hour on cut steps
    dtmin       = fineDt/20;
else
    max_retries = 12;
    dtmin       = fineDt/100;
end
dtmax        = fineDt;
symmetric_stiffness = 0;
pressureScale = 1;

%% ---- Holzapfel-Ogden passive parameters (per region) ----

a_mat_LV  = 0.186973;  b_mat_LV  = 4.80439;  af_mat_LV = 0.136467;  bf_mat_LV = 4.15537;  as_mat_LV = 1.33191;  bs_mat_LV = 3.334;
a_mat_RV  = 0.336956;  b_mat_RV  = 5.68276;  af_mat_RV = 0.352679;  bf_mat_RV = 6;        as_mat_RV = 1.72124;  bs_mat_RV = 5.85127;
a_mat_S   = 0.127484;  b_mat_S   = 8.49002;  af_mat_S  = 0.642743;  bf_mat_S  = 3.33202;  as_mat_S  = 1.41397;  bs_mat_S  = 5.2872;
k_mat     = 2000;

%% ---- Active contraction parameter (built-in material) ----
% T0: the prescribed transverse active stress magnitude [P]. Using the
% same 98 kPa peak value as Tmax in the custom-plugin run, uniformly
% across LV/RV/S, so the two approaches are comparable on a like-for-like
% basis -- adjust here if you want a different scale for this material.
% Driven by LC_calcium (id=2, real fluorescence data) below -- NO
% length-dependence, unlike the custom plugin's Ca0/beta/l0/refl terms.
% activeType:
%   'fiberstress' = FEBio 'uncoupled active fiber stress' (has force-length and
%                   force-velocity curves; the material your advisor suggested).
%                   NOTE: field names follow the lab's FiberStressTest.m, whose own
%                   header says the syntax is UNVERIFIED. If FEBio rejects the file,
%                   send me the first error lines.
%   'transiso'    = 'uncoupled prescribed trans iso active contraction' (T0 only).
%                   Run B showed this one INFLATES the cavity at zero pressure.
activeType = 'fiberstress';
tES        = 0.06;    % time pressure reaches ESP in 'ED' mode (0.10 = original timing)
smax       = 100;     % peak active fiber stress (kPa), used by 'fiberstress'
T0_peak    = 50;      % used by 'transiso'
zeroPressureTest = false;   % true = active stress only, no pressure (direction test). Set false for real runs.

%% ---- Pressure (mmHg -> kPa) ----
lowPressureTest = false;   % true = diagnostic run at about half pressure; set false for real ESP
if lowPressureTest
    Pressure_LVRV = [50.0  11.0];
else
    Pressure_LVRV = [100.0  21.9];   % [LV, RV] mmHg, END-SYSTOLIC
end
P_LV = Pressure_LVRV(1) * 0.133;
P_RV = Pressure_LVRV(2) * 0.133;

%% ---- File names ----
% Distinct "_TransIso" suffix so this never overwrites the custom-plugin
% active run's .feb/.xplt/log files in the same savePath directory.
febioFebFileNamePart = [ratname '_Active_TransIso'];
febioFebFileName     = fullfile(savePath, [febioFebFileNamePart '.feb']);
febioLogFileName     = [febioFebFileNamePart '.txt'];
febioLogFileName_disp = [febioFebFileNamePart '_disp_out.txt'];
febioLogFileName_stress = [febioFebFileNamePart '_stress_out.txt'];
febioLogFileName_principalstress = [febioFebFileNamePart '_prinstress_out.txt'];

%% ---- Load mesh ----
parentdir = 'C:\Users\shreygandhi\Desktop\Volume_Meshes';
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
LV = sort(unique(readmatrix(fullfile(baseDir, 'LV_geo_clean.csv'))));
RV = sort(unique(readmatrix(fullfile(baseDir, 'RV_geo_clean.csv'))));
S  = sort(unique(readmatrix(fullfile(baseDir,  'S_geo_clean.csv'))));
fprintf('Regions: LV=%d  RV=%d  S=%d\n', numel(LV), numel(RV), numel(S));

%% ---- Boundary surfaces ----
Fb = patchNormalFix(Fb);
F_base_BC     = Fb(Cb==1,:);
bcSupportList = unique(F_base_BC(:));
% Same flipped (non-fliplr'd) pressure surfaces as the custom-plugin
% active script -- genuinely inflating convention, verified numerically.
F_LV_pressure = patchNormalFix(Fb(Cb==3,:));
F_RV_pressure = patchNormalFix(Fb(Cb==4,:));

%% ---- Load ED (pressurized) reference geometry ----
% Same scope as the custom-plugin script: contraction only, starting from
% the passive script's own properly-computed ED baseline.
geomFile = fullfile(savePath, [ratname '_geomStates.mat']);
if ~exist(geomFile,'file')
    error('Geometry file not found: %s\nRun the passive script first.', geomFile);
end
geomData = load(geomFile);
if strcmp(startMode,'ED')
    V_def = geomData.V_ED_final;    % passive script's loaded ED_sim state
else
    V_def = geomData.V_UZP_final;   % stress-free start; phase 1 inflates to ED
end
fprintf('Start geometry (%s): LV=%.1f  RV=%.1f uL\n', startMode, ...
    closeAndVolume(F_LV_pressure, V_def), closeAndVolume(F_RV_pressure, V_def));
fprintf('Passive-script volumes: UZP LV=%.1f RV=%.1f | ED_sim LV=%.1f RV=%.1f | MRI ED LV=%.1f RV=%.1f\n', ...
    closeAndVolume(F_LV_pressure, geomData.V_UZP_final), closeAndVolume(F_RV_pressure, geomData.V_UZP_final), ...
    closeAndVolume(F_LV_pressure, geomData.V_ED_final),  closeAndVolume(F_RV_pressure, geomData.V_ED_final), ...
    closeAndVolume(F_LV_pressure, geomData.V_MRI_ED),    closeAndVolume(F_RV_pressure, geomData.V_MRI_ED));
fprintf('Loaded ED (pressurized) reference geometry: %d nodes\n', size(V_def,1));

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
febio_spec.Control.solver.rtol  = 0.001;
febio_spec.Control.solver.etol  = 0.01;
febio_spec.Control.solver.dtol  = 0.001;
febio_spec.Control.solver.lstol = 0.9;
febio_spec.Control.time_stepper.dtmin       = dtmin;
febio_spec.Control.time_stepper.dtmax       = dtmax;
febio_spec.Control.time_stepper.max_retries = max_retries;
febio_spec.Control.time_stepper.opt_iter    = opt_iter;

%% ---- Materials: 'uncoupled solid mixture' (Holzapfel_Ogden + built-in active) ----
% Each region gets its own mixture material with its own passive
% constants (matching the passive script exactly) plus the shared T0.
%
% FIX: k_mat now lives on the OUTER uncoupled solid mixture material, not
% nested inside the Holzapfel_Ogden solid block, per the plugin's own
% documented guidance (see header note above and FiberStressTest.m).
materialName1 = 'Material1';
febio_spec.Material.material{1}.ATTR.name = materialName1;
febio_spec.Material.material{1}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{1}.ATTR.id   = 1;
febio_spec.Material.material{1}.k = k_mat;   % <-- moved here (was nested in solid{1})
febio_spec.Material.material{1}.solid{1}.ATTR.type = 'Holzapfel_Ogden';
febio_spec.Material.material{1}.solid{1}.a    = a_mat_LV;
febio_spec.Material.material{1}.solid{1}.b    = b_mat_LV;
febio_spec.Material.material{1}.solid{1}.af   = af_mat_LV;
febio_spec.Material.material{1}.solid{1}.bf   = bf_mat_LV;
febio_spec.Material.material{1}.solid{1}.as   = as_mat_LV;
febio_spec.Material.material{1}.solid{1}.bs   = bs_mat_LV;
febio_spec.Material.material{1}.solid{1}.afs  = 0.0;
febio_spec.Material.material{1}.solid{1}.bfs  = 0.0;
febio_spec.Material.material{1}.solid{1}.asn  = 0.0;
febio_spec.Material.material{1}.solid{1}.bsn  = 0.0;
febio_spec.Material.material{1}.solid{1}.anf  = 0.0;
febio_spec.Material.material{1}.solid{1}.bnf  = 0.0;
febio_spec.Material.material{1}.solid{2}.ATTR.type = 'uncoupled prescribed trans iso active contraction';
febio_spec.Material.material{1}.solid{2}.T0.ATTR.lc = 2;
febio_spec.Material.material{1}.solid{2}.T0.VAL     = T0_peak;

materialName2 = 'Material2';
febio_spec.Material.material{2}.ATTR.name = materialName2;
febio_spec.Material.material{2}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{2}.ATTR.id   = 2;
febio_spec.Material.material{2}.k = k_mat;   % <-- moved here (was nested in solid{1})
febio_spec.Material.material{2}.solid{1}.ATTR.type = 'Holzapfel_Ogden';
febio_spec.Material.material{2}.solid{1}.a    = a_mat_RV;
febio_spec.Material.material{2}.solid{1}.b    = b_mat_RV;
febio_spec.Material.material{2}.solid{1}.af   = af_mat_RV;
febio_spec.Material.material{2}.solid{1}.bf   = bf_mat_RV;
febio_spec.Material.material{2}.solid{1}.as   = as_mat_RV;
febio_spec.Material.material{2}.solid{1}.bs   = bs_mat_RV;
febio_spec.Material.material{2}.solid{1}.afs  = 0.0;
febio_spec.Material.material{2}.solid{1}.bfs  = 0.0;
febio_spec.Material.material{2}.solid{1}.asn  = 0.0;
febio_spec.Material.material{2}.solid{1}.bsn  = 0.0;
febio_spec.Material.material{2}.solid{1}.anf  = 0.0;
febio_spec.Material.material{2}.solid{1}.bnf  = 0.0;
febio_spec.Material.material{2}.solid{2}.ATTR.type = 'uncoupled prescribed trans iso active contraction';
febio_spec.Material.material{2}.solid{2}.T0.ATTR.lc = 2;
febio_spec.Material.material{2}.solid{2}.T0.VAL     = T0_peak;

materialName3 = 'Material3';
febio_spec.Material.material{3}.ATTR.name = materialName3;
febio_spec.Material.material{3}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{3}.ATTR.id   = 3;
febio_spec.Material.material{3}.k = k_mat;   % <-- moved here (was nested in solid{1})
febio_spec.Material.material{3}.solid{1}.ATTR.type = 'Holzapfel_Ogden';
febio_spec.Material.material{3}.solid{1}.a    = a_mat_S;
febio_spec.Material.material{3}.solid{1}.b    = b_mat_S;
febio_spec.Material.material{3}.solid{1}.af   = af_mat_S;
febio_spec.Material.material{3}.solid{1}.bf   = bf_mat_S;
febio_spec.Material.material{3}.solid{1}.as   = as_mat_S;
febio_spec.Material.material{3}.solid{1}.bs   = bs_mat_S;
febio_spec.Material.material{3}.solid{1}.afs  = 0.0;
febio_spec.Material.material{3}.solid{1}.bfs  = 0.0;
febio_spec.Material.material{3}.solid{1}.asn  = 0.0;
febio_spec.Material.material{3}.solid{1}.bsn  = 0.0;
febio_spec.Material.material{3}.solid{1}.anf  = 0.0;
febio_spec.Material.material{3}.solid{1}.bnf  = 0.0;
febio_spec.Material.material{3}.solid{2}.ATTR.type = 'uncoupled prescribed trans iso active contraction';
febio_spec.Material.material{3}.solid{2}.T0.ATTR.lc = 2;
febio_spec.Material.material{3}.solid{2}.T0.VAL     = T0_peak;

%% ---- Optional: swap the active component for the fiber-stress material ----
if strcmp(activeType,'fiberstress')
    for mi = 1:3
        febio_spec.Material.material{mi}.solid{2} = struct();
        febio_spec.Material.material{mi}.solid{2}.ATTR.type = 'uncoupled active fiber stress';
        febio_spec.Material.material{mi}.solid{2}.smax = smax;
        febio_spec.Material.material{mi}.solid{2}.activation.ATTR.lc = 2;
        febio_spec.Material.material{mi}.solid{2}.activation.VAL     = 1;
        febio_spec.Material.material{mi}.solid{2}.stl.ATTR.type      = 'point';
        febio_spec.Material.material{mi}.solid{2}.stl.interpolate    = 'linear';
        febio_spec.Material.material{mi}.solid{2}.stl.points.pt.VAL  = [0.7 0.3; 0.85 0.7; 1.0 1.0; 1.15 0.7; 1.3 0.3];
        febio_spec.Material.material{mi}.solid{2}.stv.ATTR.type      = 'point';
        febio_spec.Material.material{mi}.solid{2}.stv.interpolate    = 'linear';
        febio_spec.Material.material{mi}.solid{2}.stv.points.pt.VAL  = [-1 0.5; 0 1; 1 1.5];
    end
end

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
% Same local coordinate system used everywhere -- the manual confirms this
% built-in material takes e1 from mat_axis (or the global system if none
% is specified), so no special handling is needed for the new material.
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
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc  = 3;
febio_spec.Loads.surface_load{2}.pressure.VAL       = P_RV * pressureScale;
febio_spec.Loads.surface_load{2}.symmetric_stiffness = 0;

%% ---- Load curves ----
% Same pressure ramp shapes (retimed to track calcium, peak-then-decline)
% as the custom-plugin active script.
fED_LV = 0.0583;   % EDP/ESP for LV (5.83/100 mmHg)
fED_RV = 0.0680;   % EDP/ESP for RV (1.49/21.9 mmHg)
Ca_orig = [
    0.0000 0.0000; 0.0240 0.0790; 0.0440 0.4415; 0.0640 0.8141;
    0.0840 0.9814; 0.0960 1.0000; 0.1240 0.9619; 0.1640 0.9022;
    0.2040 0.8331; 0.2640 0.6752; 0.3240 0.4902; 0.3840 0.3339;
    0.4640 0.2006; 0.5440 0.1220; 0.6240 0.0755; 0.7040 0.0471;
    0.7840 0.0267; 0.8640 0.0105];
if strcmp(startMode,'ED')
    % Already at ED: apply only the pressure above EDP; calcium starts at t=0.
    % tES = time pressure reaches ESP (original curve used 0.10). Earlier = ESP lines up
    % with the contraction-driven volume minimum instead of arriving after it.
    % Pressure starts rising at t=0, in step with calcium (no 0.01 hold). In v2 the hold let
    % calcium shrink the unloaded cavity 508 -> 394 before any pressure existed.
    tp = [0.0 0.044 0.07 0.10];
    tp = tp * tES / 0.10;            % compress the whole pressure rise so ESP lands at tES
    % After peak, pressure relaxes toward EDP (real hearts do this). In v3 it held at 85% of ESP
    % while calcium decayed, so the weakened wall ballooned (negative jacobians at t~0.14 to 0.17).
    tEnd = tp(4) + 0.04;
    tRel = tp(4) + 0.10;
    PLV = [0 0; tp(2) 0.25-fED_LV; tp(3) 0.60-fED_LV; tp(4) 1.0-fED_LV; tEnd 0.85-fED_LV; tRel 0.25-fED_LV; 8.45 0.25-fED_LV];
    PRV = [0 0; tp(2) 0.25-fED_RV; tp(3) 0.60-fED_RV; tp(4) 1.0-fED_RV; tEnd 0.85-fED_RV; tRel 0.25-fED_RV; 8.45 0.25-fED_RV];
    % Calcium time axis compressed so its peak (orig t=0.096) lands at tES. In v2 calcium peaked at
    % 0.096 while pressure peaked at 0.06, so volume kept falling until step 42 (misaligned).
    % Rise follows the original shape (compressed to peak at tES). The decay is tied to the pressure
    % relaxation: calcium about 5 to 10% above the pressure fraction. In v5 calcium was still 0.69
    % at t=0.16 while pressure was 0.19, so the heart kept emptying (volume fell 256 -> 227).
    iRise = find(Ca_orig(:,1) * tES / 0.096 <= tES + 1e-9);
    CaRise = [Ca_orig(iRise,1) * tES / 0.096, Ca_orig(iRise,2)];
    CaDecay = [tES + 0.01, 0.97; tEnd, 0.90; tp(4) + 0.07, 0.55; tRel, 0.25; tRel + 0.04, 0.05; 8.45, 0.0];
    CaC = [CaRise; CaDecay];
else
    % Two-phase: fill to EDP by t=0.05 with no calcium, then systole.
    PLV = [0 0; 0.05 fED_LV; 0.084 0.25; 0.11 0.60; 0.14 1.0; 0.24 0.85; 8.45 0.85];
    PRV = [0 0; 0.05 fED_RV; 0.084 0.25; 0.11 0.60; 0.14 1.0; 0.24 0.85; 8.45 0.85];
    CaC = [0 0; Ca_orig(:,1)+0.05, Ca_orig(:,2)];
end
if zeroPressureTest
    PLV(:,2) = 0; PRV(:,2) = 0;   % diagnostic: no pressure at all
end
febio_spec.LoadData.load_controller{1}.ATTR.name   = 'LC_pressure';
febio_spec.LoadData.load_controller{1}.ATTR.id     = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type   = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
febio_spec.LoadData.load_controller{1}.points.pt.VAL = PLV;
febio_spec.LoadData.load_controller{3}.ATTR.name   = 'LC_pressure_RV';
febio_spec.LoadData.load_controller{3}.ATTR.id     = 3;
febio_spec.LoadData.load_controller{3}.ATTR.type   = 'loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate = 'LINEAR';
febio_spec.LoadData.load_controller{3}.points.pt.VAL = PRV;
febio_spec.LoadData.load_controller{2}.ATTR.name   = 'LC_calcium';
febio_spec.LoadData.load_controller{2}.ATTR.id     = 2;
febio_spec.LoadData.load_controller{2}.ATTR.type   = 'loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate = 'LINEAR';
febio_spec.LoadData.load_controller{2}.points.pt.VAL = CaC;

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

    if isfield(dataDisp, 'time') && numel(dataDisp.time) == nSteps
        stepTimes = dataDisp.time(:);
    else
        stepTimes = simDuration * (1:nSteps)' / nSteps;
    end
    [~, iED] = min(abs(stepTimes - tED));
    Vol_LV_ED = Vol_LV_all(iED);
    Vol_RV_ED = Vol_RV_all(iED);
    fprintf('ED taken at step %d (t=%.4f)\n', iED, stepTimes(iED));

    % ES = step where pressure peaks (ESP), per the calibration goal
    interpLC = @(pts,t) interp1(pts(:,1), pts(:,2), t, 'linear', 'extrap');
    [~, k_LV] = max(interpLC(febio_spec.LoadData.load_controller{1}.points.pt.VAL, stepTimes));
    [~, k_RV] = max(interpLC(febio_spec.LoadData.load_controller{3}.points.pt.VAL, stepTimes));
    Vol_LV_ES = Vol_LV_all(k_LV);
    Vol_RV_ES = Vol_RV_all(k_RV);

    SV_LV = Vol_LV_ED - Vol_LV_ES;
    SV_RV = Vol_RV_ED - Vol_RV_ES;

    fprintf('\n=== Results ===\n');
    fprintf('LV: EDV=%.1f  ESV=%.1f (at step %d)  SV=%.1f  EF=%.1f%%\n', ...
        Vol_LV_ED, Vol_LV_ES, k_LV, SV_LV, SV_LV/Vol_LV_ED*100);
    fprintf('RV: EDV=%.1f  ESV=%.1f (at step %d)  SV=%.1f  EF=%.1f%%\n', ...
        Vol_RV_ED, Vol_RV_ES, k_RV, SV_RV, SV_RV/Vol_RV_ED*100);


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

    %% ---- CALIBRATION DIAGNOSTIC (new) ----
    % Checks: (1) does step 1 land close to the analytic ED target,
    % (2) does the step of minimum volume (simulated ES) coincide with
    % the step where pressure reaches its target level.
    ESP_target_LV_mmHg = Pressure_LVRV(1);   % reuse the same targets already set above
    ESP_target_RV_mmHg = Pressure_LVRV(2);

    interpLC = @(pts,t) interp1(pts(:,1), pts(:,2), t, 'linear', 'extrap');
    pressLV_pts = febio_spec.LoadData.load_controller{1}.points.pt.VAL;
    pressRV_pts = febio_spec.LoadData.load_controller{3}.points.pt.VAL;
    calcium_pts = febio_spec.LoadData.load_controller{2}.points.pt.VAL;

    fracLV = interpLC(pressLV_pts, stepTimes);
    fracRV = interpLC(pressRV_pts, stepTimes);
    fracCa = interpLC(calcium_pts, stepTimes);

    iP_LV = find(fracLV >= 0.95*max(fracLV), 1, 'first');
    iP_RV = find(fracRV >= 0.95*max(fracRV), 1, 'first');
    % ESV is searched only up to a few steps past peak pressure. Searching the whole run picked
    % the relaxation tail (v5: step 46), where pressure is already low, which is not end systole.
    winEnd_LV = nSteps; if ~isempty(iP_LV), winEnd_LV = min(nSteps, iP_LV + 3); end
    winEnd_RV = nSteps; if ~isempty(iP_RV), winEnd_RV = min(nSteps, iP_RV + 3); end
    [minVol_LV, iES_LV] = min(Vol_LV_all(iED+1:winEnd_LV)); iES_LV = iES_LV + iED;
    [minVol_RV, iES_RV] = min(Vol_RV_all(iED+1:winEnd_RV)); iES_RV = iES_RV + iED;

    fprintf('\n=== CALIBRATION DIAGNOSTIC ===\n');
    fprintf('%4s  %8s  %8s  %8s  %8s  %10s  %10s\n', ...
        'step','t','fracP_LV','fracP_RV','fracCa','Vol_LV','Vol_RV');
    for k = 1:nSteps
        fprintf('%4d  %8.4f  %8.3f  %8.3f  %8.3f  %10.1f  %10.1f\n', ...
            k, stepTimes(k), fracLV(k), fracRV(k), fracCa(k), Vol_LV_all(k), Vol_RV_all(k));
    end

    fprintf('\n--- ED (end of fill) vs analytic ED target ---\n');
    fprintf('LV: sim=%.1f  target=%.1f  (delta=%.1f, %.2f%%)\n', ...
        Vol_LV_ED, ED_target_LV, Vol_LV_ED-ED_target_LV, 100*(Vol_LV_ED-ED_target_LV)/ED_target_LV);
    fprintf('RV: sim=%.1f  target=%.1f  (delta=%.1f, %.2f%%)\n', ...
        Vol_RV_ED, ED_target_RV, Vol_RV_ED-ED_target_RV, 100*(Vol_RV_ED-ED_target_RV)/ED_target_RV);

    fprintf('\n--- Alignment: step of min volume (ESV) vs step pressure reaches target ---\n');
    if isempty(iP_LV)
        fprintf('LV: ESV at step %d (t=%.4f)   pressure NEVER reaches target fraction\n', iES_LV, stepTimes(iES_LV));
    else
        fprintf('LV: ESV at step %d (t=%.4f)   pressure target reached at step %d (t=%.4f)%s\n', ...
            iES_LV, stepTimes(iES_LV), iP_LV, stepTimes(iP_LV), ...
            merge(abs(iES_LV-iP_LV)<=3, '  [ALIGNED within 3 steps]', '  [MISALIGNED]'));
    end
    if isempty(iP_RV)
        fprintf('RV: ESV at step %d (t=%.4f)   pressure NEVER reaches target fraction\n', iES_RV, stepTimes(iES_RV));
    else
        fprintf('RV: ESV at step %d (t=%.4f)   pressure target reached at step %d (t=%.4f)%s\n', ...
            iES_RV, stepTimes(iES_RV), iP_RV, stepTimes(iP_RV), ...
            merge(abs(iES_RV-iP_RV)<=3, '  [ALIGNED within 3 steps]', '  [MISALIGNED]'));
    end

    if ~isempty(iP_LV) && ~isempty(iP_RV)
        fprintf('\nEF at peak pressure step:  LV = %.1f%%  (EDV %.1f -> V %.1f)   RV = %.1f%%  (EDV %.1f -> V %.1f)\n', ...
            100*(Vol_LV_ED-Vol_LV_all(iP_LV))/Vol_LV_ED, Vol_LV_ED, Vol_LV_all(iP_LV), ...
            100*(Vol_RV_ED-Vol_RV_all(iP_RV))/Vol_RV_ED, Vol_RV_ED, Vol_RV_all(iP_RV));
    end

    figure('Name','Volume trajectory (trans-iso active material) vs. true ED target');
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

function out = merge(cond, a, b)
    if cond, out = a; else, out = b; end
end
