clear all
clc
close all
%% ========================================================================
%  ACTIVE CONTRACTION MODEL -- FEBio built-in "uncoupled active fiber
%  stress" variant
%
%  Everything in this script (solver settings, geometry loading, pressure
%  load, load curves, output/diagnostics, results parsing) is copied
%  UNCHANGED from the working ACTIVE_FEA_Model_Only.m (custom
%  Holzapfel-Ogden_ACTIVE plugin version). The ONLY change is the
%  Material section: instead of a single custom plugin class that hand-
%  computes passive HO + Hill/calcium active stress together, this
%  version combines two SEPARATE materials via FEBio's built-in
%  "uncoupled solid mixture" container, per region:
%    1) 'Holzapfel_Ogden'            -- the ORIGINAL, untouched passive
%                                        plugin (same one validated by
%                                        the passive model), unmodified.
%    2) 'uncoupled active fiber stress' -- FEBio BUILT-IN material
%                                        (Section 4.12.1.4 of the FEBio
%                                        user manual), Hill-type active
%                                        fiber stress WITH force-velocity
%                                        (rate) dependence via sTV(l_dot)
%                                        -- the physics missing from the
%                                        custom plugin's Ca0/length-only
%                                        formulation.
%
%  Per the Holzapfel_Ogden plugin's own theory manual: "as any other
%  uncoupled material, the Holzapfel-Ogden model can be included within
%  an uncoupled solid mixture. In that case, the bulk modulus should be
%  outside any solid domain." -- i.e. k belongs on the OUTER 'uncoupled
%  solid mixture' material, NOT inside the nested Holzapfel_Ogden solid
%  block. Implemented that way below.
%
%  *** UNVERIFIED / TO CONFIRM BEFORE TRUSTING RESULTS ***
%  The FEBio manual text this was drafted from lost its XML tag
%  structure when copied out of the PDF (parameter VALUES came through,
%  but the literal element names did not, e.g. it's unclear whether the
%  activation-level parameter is tagged <activation> or something else,
%  and the exact tags/structure for the sTL/sTV functions -- math
%  expression vs. point-list form). The field names used below
%  (smax, activation, fiber, sTL, sTV) are best-effort guesses following
%  standard FEBio4 material parameter naming conventions, NOT confirmed
%  against a working example. Before trusting any run:
%    1) Open the generated .feb in FEBioStudio (or grep it directly) and
%       check the 'uncoupled active fiber stress' block against
%       FEBioStudio's material parameter list / the FEBio user manual
%       (Section 4.12.1.4) for that material.
%    2) Confirm whether an explicit <fiber> direction is needed inside
%       the active material block, or whether (like Holzapfel_Ogden) it
%       correctly inherits orientation from the mat_axis MeshData already
%       defined below -- left OUT here on the assumption it inherits,
%       matching how every other material in this project gets its fiber
%       direction.
%    3) sTL/sTV below are PLACEHOLDER functions (simple linear ramps),
%       NOT literature-derived force-length/force-velocity curves --
%       replace with real data (e.g. de Tombe & Stienen, or the
%       Krueger/Tsujioka force-velocity data already used elsewhere in
%       this project) before treating results as meaningful.
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

%% ---- Solver settings (unchanged from ACTIVE_FEA_Model_Only.m) ----
simDuration  = 0.45;
numTimeSteps = 20;
fineDt       = simDuration/numTimeSteps;
baseDt       = fineDt;
max_refs     = 25;
max_ups      = 10;
opt_iter     = 25;
max_retries  = 16;
dtmin        = fineDt/100;
dtmax        = fineDt;
symmetric_stiffness = 0;
pressureScale = 1;

%% ---- Holzapfel-Ogden passive parameters (per region, unchanged) ----
a_mat_LV  = 0.186973;  b_mat_LV  = 4.80439;  af_mat_LV = 0.136467;  bf_mat_LV = 4.15537;  as_mat_LV = 1.33191;  bs_mat_LV = 3.334;
a_mat_RV  = 0.336956;  b_mat_RV  = 5.68276;  af_mat_RV = 0.352679;  bf_mat_RV = 6;        as_mat_RV = 1.72124;  bs_mat_RV = 5.85127;
a_mat_S   = 0.127484;  b_mat_S   = 8.49002;  af_mat_S  = 0.642743;  bf_mat_S  = 3.33202;  as_mat_S  = 1.41397;  bs_mat_S  = 5.2872;
k_mat     = 2000;   % now lives on the OUTER 'uncoupled solid mixture', per the HO plugin's own guidance

%% ---- Active fiber stress parameters (per region) ----
% smax: scale factor for maximum active stress [kPa]. Starting from the
% same Tmax=98 confirmed-stable ceiling used in the custom-plugin model,
% since that's the only empirically-tested magnitude we have so far for
% this mesh/solver combination -- NOT assumed to still be the right
% ceiling with a different material (may need its own stability sweep).
smax_LV = 98;
smax_RV = 98;
smax_S  = 98;

%% ---- Pressure (mmHg -> kPa), unchanged ----
Pressure_LVRV = [100.0  21.9];
P_LV = Pressure_LVRV(1) * 0.133;
P_RV = Pressure_LVRV(2) * 0.133;

%% ---- File names ----
febioFebFileNamePart = [ratname '_ActiveFiberStress'];
febioFebFileName     = fullfile(savePath, [febioFebFileNamePart '.feb']);
febioLogFileName     = [febioFebFileNamePart '.txt'];
febioLogFileName_disp = [febioFebFileNamePart '_disp_out.txt'];
febioLogFileName_stress = [febioFebFileNamePart '_stress_out.txt'];
febioLogFileName_principalstress = [febioFebFileNamePart '_prinstress_out.txt'];

%% ---- Load mesh (unchanged) ----
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

%% ---- Load region CSVs (unchanged) ----
LV = sort(unique(readmatrix(fullfile(baseDir, 'LV_freewall.csv'))));
RV = sort(unique(readmatrix(fullfile(baseDir, 'RV_freewall.csv'))));
S  = sort(unique(readmatrix(fullfile(baseDir,  'S_freewall.csv'))));
fprintf('Regions: LV=%d  RV=%d  S=%d\n', numel(LV), numel(RV), numel(S));

%% ---- Boundary surfaces (unchanged) ----
Fb = patchNormalFix(Fb);
F_base_BC     = Fb(Cb==1,:);
bcSupportList = unique(F_base_BC(:));
F_LV_pressure = patchNormalFix(fliplr(Fb(Cb==3,:)));
F_RV_pressure = patchNormalFix(fliplr(Fb(Cb==4,:)));

%% ---- Load ED (pressurized) reference geometry (unchanged) ----
geomFile = fullfile(savePath, [ratname '_geomStates.mat']);
if ~exist(geomFile,'file')
    error('Geometry file not found: %s\nRun the passive script first.', geomFile);
end
geomData = load(geomFile);
V_def = geomData.V_ED_final;
fprintf('Loaded ED (pressurized) reference geometry: %d nodes\n', size(V_def,1));

%% ---- FEBio spec (unchanged solver block) ----
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

%% ---- Materials: 'uncoupled solid mixture' = Holzapfel_Ogden + uncoupled active fiber stress ----
% *** SEE UNVERIFIED-FIELDS WARNING AT TOP OF FILE BEFORE TRUSTING RESULTS ***
materialName1 = 'Material1';
febio_spec.Material.material{1}.ATTR.name = materialName1;
febio_spec.Material.material{1}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{1}.ATTR.id   = 1;
febio_spec.Material.material{1}.k = k_mat;   % bulk modulus lives on the mixture, not the nested HO solid
% -- passive component --
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
% -- active component (built-in, velocity-dependent) --
% Bisection test confirmed passive-only (Holzapfel_Ogden inside uncoupled
% solid mixture) runs stably -- the wrapper is NOT the problem. Re-enabling
% the active material with the corrected tent-shaped stl (peaked at l=1,
% falls off both sides) instead of the original backwards ramp.
febio_spec.Material.material{1}.solid{2}.ATTR.type = 'uncoupled active fiber stress';
febio_spec.Material.material{1}.solid{2}.smax = smax_LV;
febio_spec.Material.material{1}.solid{2}.activation.ATTR.lc = 2;
febio_spec.Material.material{1}.solid{2}.activation.VAL     = 1;
febio_spec.Material.material{1}.solid{2}.stl.ATTR.type      = 'point';
febio_spec.Material.material{1}.solid{2}.stl.interpolate    = 'linear';
febio_spec.Material.material{1}.solid{2}.stl.points.pt.VAL  = [0.7 0.3; 0.85 0.7; 1.0 1.0; 1.15 0.7; 1.3 0.3];
febio_spec.Material.material{1}.solid{2}.stv.ATTR.type      = 'point';
febio_spec.Material.material{1}.solid{2}.stv.interpolate    = 'linear';
febio_spec.Material.material{1}.solid{2}.stv.points.pt.VAL  = [-1 0.5; 0 1; 1 1.5];

materialName2 = 'Material2';
febio_spec.Material.material{2}.ATTR.name = materialName2;
febio_spec.Material.material{2}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{2}.ATTR.id   = 2;
febio_spec.Material.material{2}.k = k_mat;
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
febio_spec.Material.material{2}.solid{2}.ATTR.type = 'uncoupled active fiber stress';
febio_spec.Material.material{2}.solid{2}.smax = smax_RV;
febio_spec.Material.material{2}.solid{2}.activation.ATTR.lc = 2;
febio_spec.Material.material{2}.solid{2}.activation.VAL     = 1;
febio_spec.Material.material{2}.solid{2}.stl.ATTR.type      = 'point';
febio_spec.Material.material{2}.solid{2}.stl.interpolate    = 'linear';
febio_spec.Material.material{2}.solid{2}.stl.points.pt.VAL  = [0.7 0.3; 0.85 0.7; 1.0 1.0; 1.15 0.7; 1.3 0.3];
febio_spec.Material.material{2}.solid{2}.stv.ATTR.type      = 'point';
febio_spec.Material.material{2}.solid{2}.stv.interpolate    = 'linear';
febio_spec.Material.material{2}.solid{2}.stv.points.pt.VAL  = [-1 0.5; 0 1; 1 1.5];

materialName3 = 'Material3';
febio_spec.Material.material{3}.ATTR.name = materialName3;
febio_spec.Material.material{3}.ATTR.type = 'uncoupled solid mixture';
febio_spec.Material.material{3}.ATTR.id   = 3;
febio_spec.Material.material{3}.k = k_mat;
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
febio_spec.Material.material{3}.solid{2}.ATTR.type = 'uncoupled active fiber stress';
febio_spec.Material.material{3}.solid{2}.smax = smax_S;
febio_spec.Material.material{3}.solid{2}.activation.ATTR.lc = 2;
febio_spec.Material.material{3}.solid{2}.activation.VAL     = 1;
febio_spec.Material.material{3}.solid{2}.stl.ATTR.type      = 'point';
febio_spec.Material.material{3}.solid{2}.stl.interpolate    = 'linear';
febio_spec.Material.material{3}.solid{2}.stl.points.pt.VAL  = [0.7 0.3; 0.85 0.7; 1.0 1.0; 1.15 0.7; 1.3 0.3];
febio_spec.Material.material{3}.solid{2}.stv.ATTR.type      = 'point';
febio_spec.Material.material{3}.solid{2}.stv.interpolate    = 'linear';
febio_spec.Material.material{3}.solid{2}.stv.points.pt.VAL  = [-1 0.5; 0 1; 1 1.5];

%% ---- Elements (unchanged) ----
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

%% ---- MeshData: mat_axis (unchanged -- supplies fiber/sheet direction) ----
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

%% ---- MeshDomains (unchanged) ----
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat  = materialName1;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat  = materialName2;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.name = partName3;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat  = materialName3;

%% ---- Nodes (unchanged) ----
febio_spec.Mesh.Nodes{1}.ATTR.name    = 'Object1';
febio_spec.Mesh.Nodes{1}.node.ATTR.id = (1:size(V_def,1))';
febio_spec.Mesh.Nodes{1}.node.VAL     = V_def;

%% ---- Surfaces (unchanged) ----
febio_spec.Mesh.Surface{1}.ATTR.name    = 'LVPressure';
febio_spec.Mesh.Surface{1}.tri3.ATTR.id = (1:size(F_LV_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL     = F_LV_pressure;
febio_spec.Mesh.Surface{2}.ATTR.name    = 'RVPressure';
febio_spec.Mesh.Surface{2}.tri3.ATTR.id = (1:size(F_RV_pressure,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL     = F_RV_pressure;

%% ---- NodeSets + BCs (unchanged) ----
febio_spec.Mesh.NodeSet{1}.ATTR.name = 'bcSupportList';
febio_spec.Mesh.NodeSet{1}.VAL       = mrow(bcSupportList);
febio_spec.Boundary.bc{1}.ATTR.name     = 'FixedBase';
febio_spec.Boundary.bc{1}.ATTR.type     = 'zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set = 'bcSupportList';
febio_spec.Boundary.bc{1}.x_dof = 1;
febio_spec.Boundary.bc{1}.y_dof = 1;
febio_spec.Boundary.bc{1}.z_dof = 1;

%% ---- Pressure loads (unchanged) ----
febio_spec.Loads.surface_load{1}.ATTR.type    = 'pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface = 'LVPressure';
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc  = 1;
febio_spec.Loads.surface_load{1}.pressure.VAL       = P_LV * pressureScale;
febio_spec.Loads.surface_load{1}.symmetric_stiffness = 0;
febio_spec.Loads.surface_load{2}.ATTR.type    = 'pressure';
febio_spec.Loads.surface_load{2}.ATTR.surface = 'RVPressure';
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc  = 1;
febio_spec.Loads.surface_load{2}.pressure.VAL       = P_RV * pressureScale;
febio_spec.Loads.surface_load{2}.symmetric_stiffness = 0;

%% ---- Load curves (unchanged) ----
febio_spec.LoadData.load_controller{1}.ATTR.name       = 'LC_pressure';
febio_spec.LoadData.load_controller{1}.ATTR.id         = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type       = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate     = 'LINEAR';
febio_spec.LoadData.load_controller{1}.points.pt.VAL   = [
    0.0000   0.0000;
    0.0750   0.1563;
    0.1500   0.5000;
    0.2250   0.8438;
    0.3000   1.0000;
    0.4500   1.0000;
    8.4500   1.0000;
];
febio_spec.LoadData.load_controller{2}.ATTR.name       = 'LC_calcium';
febio_spec.LoadData.load_controller{2}.ATTR.id         = 2;
febio_spec.LoadData.load_controller{2}.ATTR.type       = 'loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate     = 'LINEAR';
febio_spec.LoadData.load_controller{2}.points.pt.VAL   = [
    0.0000   0.0000;
    0.0240   0.0790;
    0.0440   0.4415;
    0.0640   0.8141;
    0.0840   0.9814;
    0.0960   1.0000;
    0.1240   0.9619;
    0.1640   0.9022;
    0.2040   0.8331;
    0.2640   0.6752;
    0.3240   0.4902;
    0.3840   0.3339;
    0.4640   0.2006;
    0.5440   0.1220;
    0.6240   0.0755;
    0.7040   0.0471;
    0.7840   0.0267;
    0.8640   0.0105;
];

%% ---- Output (unchanged) ----
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

%% ---- Write and run (unchanged) ----
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

%% ---- Diagnostic: parse console output for negative-jacobian failures (unchanged) ----
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

%% ---- Results: scan all steps for true end-systolic (minimum) volume (unchanged) ----
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