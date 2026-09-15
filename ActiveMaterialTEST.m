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
simDuration  = 0.2;
numTimeSteps = 20;                 % matches the passive model's step count
fineDt       = simDuration/numTimeSteps;   % = 0.01, finer resolution near ES
baseDt       = fineDt;
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
% geometry 
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
Tmax  = 136;
Ca0   = 4.35;
B_LV  = 4.75;
B_RV  = 11;
B_S   = 4.75;
l0    = 1.8;
refl  = 2.20;

%% ---- Pressure (mmHg -> kPa) ----
Pressure_LVRV = [91  21.9];   % [LV, RV] mmHg, END-SYSTOLIC
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
F_LV_pressure = patchNormalFix(Fb(Cb==3,:));
F_RV_pressure = patchNormalFix(Fb(Cb==4,:));

%% ---- Load ED (pressurized) reference geometry ----
geomFile = fullfile(savePath, [ratname '_geomStates.mat']);
if ~exist(geomFile,'file')
    error('Geometry file not found: %s\nRun the passive script first.', geomFile);
end
geomData = load(geomFile);
V_def = geomData.V_ED_final;
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

%% ---- Materials: Holzapfel-Ogden_ACTIVE (custom plugin, single class) ----
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
% LC 3: RV's own ED->systolic ramp 
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
% LC 2: calcium activation curve -- built directly from real Ca2+ transient data 
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