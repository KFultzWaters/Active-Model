clear all
clc
close all
addpath(genpath('C:\Program Files\MATLAB\R2025a'))
gibbonRoot = 'C:\Users\kmagi\Downloads\GIBBON-Master';
addpath(genpath(gibbonRoot));   % add ALL subfolders
savepath;  

% Point directly to your verified FEBio executable
febioExe ='C:\Program Files\FEBioStudio\bin\febio4.exe';
febioAnalysis.febioPath = febioExe;

% Run mode: choose 'internal' (via MATLAB) or 'external' (calls FEBio exe)
febioAnalysis.runMode = 'internal';

% Add this to your FEBio analysis struct (GIBBON or your own)
febioAnalysis.febioPath = febioExe;
% % Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
markerSize=40;
lineWidth=3;

% Control parameters
% Path names
savePath=fullfile(pwd,"FeBio");


ratname = 'Z210W0';

pressureScale = 1;
%materialScale = 50;
materialScale = 5;
k = 2000;

%Material parameter set
%material_LV = [1.204551483,36.69618646,2.403175723,30,11.22519682,8.139261935,0,k];
%material_LV = [0.926163882,	43.48870923,	5.172917043,	8,	16.12644498,	60,0,k];
material_LV = [0.9348671438,	24.02194001,	0.6823358098,	20.77682741,	6.659563712,	16.66999221,0,k];
material_RV = [1.684779022,28.41380226,1.76339444,30,8.606194526,29.25633166,0,k];
material_S  = [0.637422474,42.45007529,3.213712742,16.66010499,7.069844357,26.43599904,0,k];


% Applied Pressure in mmHg
Pressure_LVRV = [5.83	1.49];
%Pressure_LVRV = [1.493674971	1.493674971];
%Pressure_LVRV = [5 25];
%Turning mmHg into kPa
P_LV = Pressure_LVRV(1)*0.133;
P_RV = Pressure_LVRV(2)*0.133;

a_mat_LV=material_LV(1);
b_mat_LV=material_LV(2);
af_mat_LV=material_LV(3);
bf_mat_LV=material_LV(4);
as_mat_LV=material_LV(5);
bs_mat_LV=material_LV(6);
k_mat_LV=material_LV(8);

a_mat_RV=material_RV(1);
b_mat_RV=material_RV(2);
af_mat_RV=material_RV(3);
bf_mat_RV=material_RV(4);
as_mat_RV=material_RV(5);
bs_mat_RV=material_RV(6);
k_mat_RV=material_RV(8);

a_mat_S=material_S(1);
b_mat_S=material_S(2);
af_mat_S=material_S(3);
bf_mat_S=material_S(4);
as_mat_S=material_S(5);
bs_mat_S=material_S(6);
k_mat_S=material_S(8);

%%
% Defining file names
febioFebFileNamePart=ratname;
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_principalstress=[febioFebFileNamePart,'_prinstress_out.txt']; %Log file name for exporting stress
febioLogFileName_strain        = [febioFebFileNamePart,'_strain_out.txt'];
febioLogFileName_principalstrain = [febioFebFileNamePart,'_prinstrain_out.txt'];
febioLogFileName_sed = [febioFebFileNamePart,'_sed_out.txt'];



% FEA control settings
analysisType='DYNAMIC';
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=10; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=6; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-3; %1e-10;
symmetric_stiffness=0;
runMode='internal';% 'internal' or 'external'

% Load Geometry
parentdir = 'C:\Users\kmagi\Downloads\PAH_Imaging_FEA\Volume_Meshes';
baseDir = fullfile(parentdir, ratname);   % adds the W282W0 subfolder
filename = fullfile(baseDir, [ratname '_With_Fibers.h5']);
data = readXDMF(filename);

V = data.Groups(2).Groups(1).Datasets(1).Value';
E = data.Groups(2).Groups(1).Datasets(2).Value';
E = double(E)+1;
Fb = data.Groups(3).Groups(2).Datasets(2).Value';
Fb = Fb+1;
Cb = data.Groups(3).Groups(2).Datasets(1).Value';
CE = data.Groups(3).Groups(1).Datasets(1).Value';

% Load Fibers
f0 = data.Groups(1).Groups(1).Datasets(1).Value';
n0 = data.Groups(1).Groups(2).Datasets(1).Value';
s0 = data.Groups(1).Groups(3).Datasets(1).Value';

meshOutput = struct();
meshOutput.nodes = V;
meshOutput.facesBoundary = Fb;
meshOutput.boundaryMarker = Cb;
meshOutput.elements = E;
meshOutput.elementMaterialID = ones(size(CE));

%meshView(meshOutput);

figure;
gpatch(Fb, V, Cb, 'k', 0.5);
axisGeom;
camlight headlight;
title('Mesh View');

%%
LV_RV_S = fullfile(parentdir, ratname);
LV_raw = readmatrix(fullfile(LV_RV_S, 'LV_freewall.csv'));
RV_raw = readmatrix(fullfile(LV_RV_S, 'RV_freewall.csv'));
S_raw  = readmatrix(fullfile(LV_RV_S, 'S_freewall.csv'));

LV = sort(unique(LV_raw));
RV = sort(unique(RV_raw));
S  = sort(unique(S_raw));

% Compute centroids
centroids = zeros(size(E,1), 3);
for q = 1:size(E,1)
    centroids(q,:) = mean(V(E(q,:),:), 1);
end

fprintf('CSV-labeled elements: LV=%d RV=%d S=%d (total %d of %d = %.1f%%)\n', ...
    numel(LV), numel(RV), numel(S), numel(LV)+numel(RV)+numel(S), size(E,1), ...
    100*(numel(LV)+numel(RV)+numel(S))/size(E,1));

% --- Region growing via true mesh connectivity (distance-weighted multi-source Dijkstra) ---
% REPLACES the old per-element single-axis-projection fallback, which
% classified each unlabeled element independently based only on its
% centroid's position along one straight LV-RV axis. That approach was
% confirmed (via visualization) to produce widespread, non-anatomical
% speckling across a large fraction of the mesh, since a single straight
% axis is a poor approximation of real curved anatomy away from mid-cavity.
%
% Instead: CSV-labeled elements are trusted seeds. Every unlabeled element
% is assigned the label of whichever seed region is truly CLOSEST via
% shortest weighted path (real Euclidean distance) through actual
% shared-face mesh connectivity. This guarantees spatially contiguous,
% smooth regions by construction, and resolves ties by genuine geometric
% distance rather than element-index order (an earlier unweighted-BFS
% attempt was order-biased and is NOT used here).
nElemTot = size(E,1);
tetFaces  = [E(:,[1 2 3]); E(:,[1 2 4]); E(:,[1 3 4]); E(:,[2 3 4])];
elemOwner = repmat((1:nElemTot)', 4, 1);
faceSorted = sort(tetFaces, 2);
[~, ~, ic] = unique(faceSorted, 'rows');
faceCounts = accumarray(ic, 1);
interiorFaceIdx = find(faceCounts == 2);
isInteriorRow = ismember(ic, interiorFaceIdx);
ic_int    = ic(isInteriorRow);
owner_int = elemOwner(isInteriorRow);
[ic_int_sorted, sIdx] = sort(ic_int);
owner_int_sorted = owner_int(sIdx);
adjPairs = reshape(owner_int_sorted, 2, [])';   % [nInteriorFaces x 2] element adjacency

edgeWeights = vecnorm(centroids(adjPairs(:,1),:) - centroids(adjPairs(:,2),:), 2, 2);

nVirtual = 3;
srcLVidx = nElemTot + 1;
srcRVidx = nElemTot + 2;
srcSidx  = nElemTot + 3;

virtEdgesFrom = [repmat(srcLVidx, numel(LV),1); repmat(srcRVidx, numel(RV),1); repmat(srcSidx, numel(S),1)];
virtEdgesTo   = [LV(:); RV(:); S(:)];
virtWeights   = zeros(numel(virtEdgesFrom),1);

allEdgesFrom = [adjPairs(:,1); virtEdgesFrom];
allEdgesTo   = [adjPairs(:,2); virtEdgesTo];
allWeights   = [edgeWeights;   virtWeights];

G = graph(allEdgesFrom, allEdgesTo, allWeights, nElemTot + nVirtual);

distToLV = distances(G, srcLVidx, 1:nElemTot);
distToRV = distances(G, srcRVidx, 1:nElemTot);
distToS  = distances(G, srcSidx,  1:nElemTot);

distAll = [distToLV(:), distToRV(:), distToS(:)];
[~, nearestRegion] = min(distAll, [], 2);

label = nearestRegion;
label(LV) = 1;   % trusted CSV seeds keep their original label
label(RV) = 2;
label(S)  = 3;

nUnreached = sum(~isfinite(min(distAll,[],2)));
if nUnreached > 0
    warning(['%d elements are disconnected from ALL seed regions via mesh connectivity ', ...
             '(infinite distance to every region). Check mesh connectivity/manifoldness.'], ...
             nUnreached);
end

LV = find(label == 1);
RV = find(label == 2);
S  = find(label == 3);

fprintf('After distance-weighted region growing - LV: %d  RV: %d  S: %d  Total: %d\n', ...
    numel(LV), numel(RV), numel(S), numel(LV)+numel(RV)+numel(S));




meshOutput.elementMaterialID(LV) = ones([length(LV),1])*1;
meshOutput.elementMaterialID(RV) = ones([length(RV),1])*2;
meshOutput.elementMaterialID(S) = ones([length(S),1])*3;

%%
% -> Fix Faces
figure; 
subplot(1,2,1); hold on;
gpatch(Fb,V,Cb,'k',1);
patchNormPlot(Fb,V);
axisGeom; camlight headlight;
drawnow;

Fb_fix = patchNormalFix(Fb);

subplot(1,2,2); hold on;
gpatch(Fb_fix,V,Cb,'k',1);
patchNormPlot(Fb_fix,V);
axisGeom; camlight headlight;
drawnow;

Fb = Fb_fix;

%%

%BCs
F_base_BC = Fb(Cb==1,:);
bcSupportList = unique(F_base_BC(:));

%PressureSurfaces
F_LV_pressure=Fb(Cb==3,:);

F_RV_pressure=Fb(Cb==4,:);
F_LV_pressure = fliplr(F_LV_pressure);
F_RV_pressure = fliplr(F_RV_pressure);

F_LV_pressure = patchNormalFix(F_LV_pressure);
F_RV_pressure = patchNormalFix(F_RV_pressure);

% Visualize BC and Load
visualizeBC(Fb,V,bcSupportList,F_LV_pressure,F_RV_pressure)

%%

Edge_LV=patchBoundary(F_LV_pressure);
[F_LV,V_LV,C_LV] = triSurfCloseHoles(double(F_LV_pressure),V,0.5,double(Edge_LV));
F_LV = patchNormalFix(F_LV);
[volEst_LV]  = patchVolume(F_LV,V_LV);

% original: if volEst_LV<0, F_LV = fliplr(F_LV); volEst_LV = -volEst_LV; end
% Vol_LV = volEst_LV;
if volEst_LV<0, volEst_LV = -volEst_LV; end
Vol_LV = volEst_LV;

Edge_RV=patchBoundary(F_RV_pressure);
[F_RV,V_RV,C_RV] = triSurfCloseHoles(double(F_RV_pressure),V,0.5,double(Edge_RV));
F_RV = patchNormalFix(F_RV);
[volEst_RV]  = patchVolume(F_RV,V_RV);

% if volEst_RV<0, F_RV = fliplr(F_RV); volEst_RV = -volEst_RV; end
% Vol_RV = volEst_RV;
if volEst_RV<0, volEst_RV = -volEst_RV; end
Vol_RV = volEst_RV;


figure; 
hold on;
gpatch(F_LV,V_LV,C_LV);
patchNormPlot(F_LV,V_LV);
gpatch(F_RV,V_RV,C_RV);
patchNormPlot(F_RV,V_RV);
axisGeom; camlight headlight;
drawnow;


%%
%Get a template with default settings
[febio_spec]=febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version='4.0';

%Module section
febio_spec.Module.ATTR.type='solid';

%Control section
febio_spec.Control.analysis='STATIC';%'DYNAMIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.ATTR.type='Broyden';
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;


%Material section

materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Holzapfel_Ogden';
%febio_spec.Material.material{1}.ATTR.type= 'Holzapfel-Gasser-Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.a= a_mat_LV /materialScale;
febio_spec.Material.material{1}.b= b_mat_LV /materialScale;
febio_spec.Material.material{1}.af= af_mat_LV /materialScale;
febio_spec.Material.material{1}.bf= bf_mat_LV /materialScale;
febio_spec.Material.material{1}.as= as_mat_LV /materialScale;
febio_spec.Material.material{1}.bs= bs_mat_LV /materialScale;
febio_spec.Material.material{1}.afs=0.0;
febio_spec.Material.material{1}.bfs=0.0;
febio_spec.Material.material{1}.asn=0.0;
febio_spec.Material.material{1}.bsn=0.0;
febio_spec.Material.material{1}.anf=0.0;
febio_spec.Material.material{1}.bnf=0.0;
febio_spec.Material.material{1}.k= k_mat_LV;%10000.0;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Holzapfel_Ogden';
%febio_spec.Material.material{2}.ATTR.type= 'Holzapfel-Gasser-Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.a= a_mat_RV /materialScale;
febio_spec.Material.material{2}.b= b_mat_RV /materialScale;
febio_spec.Material.material{2}.af= af_mat_RV /materialScale;
febio_spec.Material.material{2}.bf= bf_mat_RV /materialScale;
febio_spec.Material.material{2}.as= as_mat_RV /materialScale;
febio_spec.Material.material{2}.bs= bs_mat_RV /materialScale;
febio_spec.Material.material{2}.afs=0.0;
febio_spec.Material.material{2}.bfs=0.0;
febio_spec.Material.material{2}.asn=0.0;
febio_spec.Material.material{2}.bsn=0.0;
febio_spec.Material.material{2}.anf=0.0;
febio_spec.Material.material{2}.bnf=0.0;
febio_spec.Material.material{2}.k= k_mat_RV;%10000.0;

materialName3='Material3';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='Holzapfel_Ogden';
%febio_spec.Material.material{3}.ATTR.type= 'Holzapfel-Gasser-Ogden';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.a= a_mat_S /materialScale;
febio_spec.Material.material{3}.b= b_mat_S /materialScale;
febio_spec.Material.material{3}.af= af_mat_S /materialScale;
febio_spec.Material.material{3}.bf= bf_mat_S /materialScale;
febio_spec.Material.material{3}.as= as_mat_S /materialScale;
febio_spec.Material.material{3}.bs= bs_mat_S /materialScale;
febio_spec.Material.material{3}.afs=0.0;
febio_spec.Material.material{3}.bfs=0.0;
febio_spec.Material.material{3}.asn=0.0;
febio_spec.Material.material{3}.bsn=0.0;
febio_spec.Material.material{3}.anf=0.0;
febio_spec.Material.material{3}.bnf=0.0;
febio_spec.Material.material{3}.k= k_mat_S;%10000.0;


% -> Elements
partName1='LV';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=LV; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E(LV,:); %The element matrix

partName2='RV';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{2}.elem.ATTR.id=RV; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E(RV,:); %The element matrix

partName3='S';
febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of this part
febio_spec.Mesh.Elements{3}.ATTR.type='tet4'; %Element type
febio_spec.Mesh.Elements{3}.elem.ATTR.id=S; %Element id's
febio_spec.Mesh.Elements{3}.elem.VAL=E(S,:); %The element matrix

% % % % % % % % % % % % %MeshData section
% -> Fibers
% -> ElementData
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

febio_spec.MeshData.ElementData{2}.ATTR.elem_set=partName2;
febio_spec.MeshData.ElementData{2}.ATTR.type='mat_axis';

febio_spec.MeshData.ElementData{3}.ATTR.elem_set=partName3;
febio_spec.MeshData.ElementData{3}.ATTR.type='mat_axis';

lid_LV = zeros(size(E,1),1); lid_LV(LV) = 1:numel(LV);
lid_RV = zeros(size(E,1),1); lid_RV(RV) = 1:numel(RV);
lid_S = zeros(size(E,1),1); lid_S(S) = 1:numel(S);

elem_count_LV = 1;
elem_count_RV = 1;
elem_count_S = 1;

for q=1:size(E,1)
    if ismember(q,LV)
        LID = lid_LV(q);
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.ATTR.lid=LID;
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.a=f0(q,:);
        febio_spec.MeshData.ElementData{1}.elem{elem_count_LV}.d=s0(q,:);
        elem_count_LV = elem_count_LV +1;
    elseif ismember(q,RV)
        LID = lid_RV(q);
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.ATTR.lid=LID;
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.a=f0(q,:);
        febio_spec.MeshData.ElementData{2}.elem{elem_count_RV}.d=s0(q,:);
        elem_count_RV = elem_count_RV +1;
    else
        LID = lid_S(q);
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.ATTR.lid=LID;
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.a=f0(q,:);
        febio_spec.MeshData.ElementData{3}.elem{elem_count_S}.d=s0(q,:);
        elem_count_S = elem_count_S +1;
    end
end


    %%

% -> Surfaces
LV_surfaceName='LVPressure';
febio_spec.Mesh.Surface{1}.ATTR.name=LV_surfaceName;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_LV_pressure,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_LV_pressure;

RV_surfaceName='RVPressure';
febio_spec.Mesh.Surface{2}.ATTR.name=RV_surfaceName;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(F_RV_pressure,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL=F_RV_pressure;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

febio_spec.MeshDomains.SolidDomain{3}.ATTR.name=partName3;
febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat=materialName3;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='FixedDisplacement01';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

%Loads section
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=LV_surfaceName;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=P_LV * pressureScale;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=0;

febio_spec.Loads.surface_load{2}.ATTR.type='pressure';
febio_spec.Loads.surface_load{2}.ATTR.surface=RV_surfaceName;
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{2}.pressure.VAL=P_RV * pressureScale; 
febio_spec.Loads.surface_load{2}.symmetric_stiffness=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sx;sy;sz;sxy;syz;sxz';
%febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
%febio_spec.Output.logfile.element_data{1}.ATTR.data='s';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_principalstress;
febio_spec.Output.logfile.element_data{2}.ATTR.data='s1;s2;s3';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';

% NEW: full Green–Lagrange strain tensor (Ex, Ey, Ez, Exy, Eyz, Exz)
febio_spec.Output.logfile.element_data{3}.ATTR.file  = febioLogFileName_strain;
febio_spec.Output.logfile.element_data{3}.ATTR.data  = 'Ex;Ey;Ez;Exy;Eyz;Exz';
febio_spec.Output.logfile.element_data{3}.ATTR.delim = ',';

% NEW: principal Green–Lagrange strains (E1, E2, E3)
febio_spec.Output.logfile.element_data{4}.ATTR.file  = febioLogFileName_principalstrain;
febio_spec.Output.logfile.element_data{4}.ATTR.data  = 'E1;E2;E3';
febio_spec.Output.logfile.element_data{4}.ATTR.delim = ',';

% NEW: strain energy density (sed) per element
febio_spec.Output.logfile.element_data{5}.ATTR.file  = febioLogFileName_sed;
febio_spec.Output.logfile.element_data{5}.ATTR.data  = 'sed';
febio_spec.Output.logfile.element_data{5}.ATTR.delim = ','; 

febio_spec.Output.plotfile.compression=0;


%% --- Initialize state and ED cavity volumes ---
% Original: (11/12)
V_original = V;              % ED geometry you loaded
V_current  = V_original;     % geometry we will update iteratively
mesh_update = 1;
V_all = V_current;           % store geometries over iterations

% Turn ED scalar volumes into arrays we can append to
Vol_LV = Vol_LV(:).';  % [ED]
Vol_RV = Vol_RV(:).';  % [ED]
% --- Target ED volumes from the original (imaged) geometry ---
ED_target_LV = closeAndVolume(F_LV_pressure, V_original);
ED_target_RV = closeAndVolume(F_RV_pressure, V_original);
fprintf('[Target ED volumes] LV=%.1f  RV=%.1f (µL)\n', ED_target_LV, ED_target_RV);


%% --- FEBio spec (parts already defined above) ---
febio_spec.Control.analysis = 'STATIC';

%% --- Iterative pressure solve + ES volume measure + inverse update ---
% Convenience: positive tet check
tetOK = @(Vtry,E) all( dot(cross(Vtry(E(:,2),:)-Vtry(E(:,1),:), ...
                                 Vtry(E(:,3),:)-Vtry(E(:,1),:),2), ...
                                Vtry(E(:,4),:)-Vtry(E(:,1),:),2) / 6 > 0 );


% === Inverse iteration controls (place before while true) ===
er_num = 0; terminate = 0;

maxOuter   = 10;     % maximum inverse iterations
%volTol     = 70;     % µL tolerance per cavity for convergence (tune) --> try to get it to 20
volTol_pct = 0.01; % 20% tolerance on ED volume
% original: beta/min/max .5/.1/.8
betaInit   = 0.25;    % relaxation for inverse update (0<beta<=1)
betaMin    = 0.1;    % lower bound on beta
betaMax    = 0.8;    % upper bound on beta
j          = 0;      % iteration counter (starts at 0; increment inside loop)
beta       = betaInit;

Vol_LV_ES_hist = [];     % ES volumes at each solve
Vol_RV_ES_hist = [];

Vdef_all   = [];         % deformed (ES) geometry per solve
Vsolve_all = [];         % geometry you solved from (unloaded guess) per solve

% (Optional) make sure Vol_LV / Vol_RV exist and have ED at index 1
assert(exist('Vol_LV','var')==1 && ~isempty(Vol_LV), 'Vol_LV not set before loop.');
assert(exist('Vol_RV','var')==1 && ~isempty(Vol_RV), 'Vol_RV not set before loop.');
solveSucceeded = false;
V_def = [];


while true
    j = j + 1;
    fprintf('\n=== Inverse iteration %d/%d ===\n', j, maxOuter);

    % --- Write current nodes to FEB file ---
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1';
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_current,1))';
    febio_spec.Mesh.Nodes{1}.node.VAL=V_current;

    % Make sure current mesh is valid (if not, backtrack toward original)
    if ~tetOK(V_current, E)
        warning('Invalid tets BEFORE run; backtracking toward original.');
        alpha_bt = 0.5;
        while ~tetOK(V_current, E) && alpha_bt > 1e-3
            V_current = V_current - alpha_bt*(V_current - V_original);
            alpha_bt  = alpha_bt*0.5;
        end
        if ~tetOK(V_current, E)
            error('Could not recover a valid mesh. Check step size / labels.');
        end
        % Rewrite nodes after backtrack
        febio_spec.Mesh.Nodes{1}.node.VAL=V_current;
    end
    % --- Right before febioStruct2xml(...) and runMonitorFEBio(...)
    V_atSolve = V_current;              % this is the geometry we actually solved from

 if ~exist(fileparts(febioFebFileName), 'dir')
    mkdir(fileparts(febioFebFileName));
end
febioStruct2xml(febio_spec, febioFebFileName);
if exist(febioFebFileName, 'file')
    fprintf('FEB file written: %d bytes\n', dir(febioFebFileName).bytes);
else
    error('febioStruct2xml failed to create file at: %s', febioFebFileName);
end

    % Write FEB file
    %febioStruct2xml(febio_spec, febioFebFileName);



    % --- Patch/normalize the run struct (ensure CHARs, set output path) ---
    outDir = fileparts(febioFebFileName);   % <-- folder of the .feb file
    febioAnalysis.run_filename    = char(febioFebFileName);
    febioAnalysis.run_output_path = char(outDir);
    [~,base,~]                    = fileparts(char(febioAnalysis.run_filename));
    febioAnalysis.run_output_name = char(base);
    febioAnalysis.fileName_plot   = char(fullfile(outDir, [base '__plot.mat']));
    febioAnalysis.run_logname     = char([base '.txt']);
    febioAnalysis.febioPath       = char(febioAnalysis.febioPath);
    febioAnalysis.runMode         = char(febioAnalysis.runMode);
    checkFEBioRunStruct(febioAnalysis);


    febioAnalysis.run_string = sprintf('"%s" -i "%s"', ...
    char(febioExe), char(febioFebFileName));

    % --- Run FEBio and handle failure with geometric backtrack ---
    [runFlag] = runMonitorFEBio(febioAnalysis);
    if runFlag ~= 1
        fprintf('\nFEBio solve failed → geometric backtrack\n');
        N_disp_from_original = V_current - V_original;
        if     er_num == 0, alpha = 0.50;
        elseif er_num <  5, alpha = 0.10;
        else                 alpha = 0.05; terminate = 1;
        end
        V_current = V_current - alpha * N_disp_from_original;
        er_num = er_num + 1;
        if terminate, warning('Stopping after repeated failures.'); break; end
        continue; % retry loop with backtracked geometry
    end
    er_num = 0;
    V_atSolve_last = V_atSolve;         % keep the last-solve reference for visuals/metrics


    % Read last-step nodal displacements; build pressurized geometry
    dataDisp   = importFEBio_logfile(fullfile(outDir, febioLogFileName_disp), 0, 1);   % <— outDir
    N_disp_mat = dataDisp.data;                       % [nNodes x 3 x nSteps]
    U_last     = N_disp_mat(:,:,end);                 % final increment
    V_def      = V_current + U_last;                  % pressurized state ("ED_sim")


   if j > 1
        [F_LV_pressure, ~, fracLV] = ensureMinusN(F_LV_pressure, V_current, U_last, 0.60);
        [F_RV_pressure, ~, fracRV] = ensureMinusN(F_RV_pressure, V_current, U_last, 0.60);
        fprintf('U·n<0 agreement (expect high %% under internal pressure):  LV=%.1f%%  RV=%.1f%%\n', ...
                100*fracLV, 100*fracRV);
    
        febio_spec.Mesh.Surface{1}.tri3.VAL = F_LV_pressure;
        febio_spec.Mesh.Surface{2}.tri3.VAL = F_RV_pressure;
    end

    % % Read last-step element stresses 
% --- Read element stresses log ---
dataStrs     = importFEBio_logfile(fullfile(outDir, febioLogFileName_stress), 0, 1);
S_elem       = dataStrs.data;                  % [nElem x nComp x nSteps]
[nElem_log, nComp, nSteps] = size(S_elem);

nElem  = size(E,1);
nNodes = size(V_current,1);
assert(nElem_log==nElem, 'Mismatch: stress log nElem != mesh nElem');

% Use the LAST step (ED_sim)
S_last = S_elem(:,:,end);                      % [nElem x nComp]

% Always have at least sx, sy, sz:
sx = S_last(:,1);
sy = S_last(:,2);
sz = S_last(:,3);

% If shear components exist, use them. Otherwise set them to zero.
if nComp >= 4
    sxy = S_last(:,4);
    syz = S_last(:,5);
    sxz = S_last(:,6);
else
    sxy = zeros(nElem,1);
    syz = zeros(nElem,1);
    sxz = zeros(nElem,1);
end

% --- Fiber-direction Cauchy stress σ_f = fᵀ σ f ---
if exist('f0','var')
    % Unify to per-element fibers: average from nodes if needed
    if size(f0,1) == nElem
        f_elem = f0;
    elseif size(f0,1) == nNodes
        f_elem = (f0(E(:,1),:) + f0(E(:,2),:) + f0(E(:,3),:) + f0(E(:,4),:)) / 4;
    else
        error('f0 has incompatible dimensions (need nElem×3 or nNodes×3).');
    end

    % Normalize fiber directions
    f_elem = f_elem ./ max(eps, vecnorm(f_elem,2,2));
    fx = f_elem(:,1); fy = f_elem(:,2); fz = f_elem(:,3);

    % σ_f per element (scalar)
    sf_elem = sx.*fx.^2 + sy.*fy.^2 + sz.*fz.^2 ...
            + 2*( sxy.*fx.*fy + sxz.*fx.*fz + syz.*fy.*fz );

    % Map element values to nodes
    fib_nodes = elemToNodeMean(E, sf_elem, nNodes);   % [nNodes x 1]

    % Store history across solves
    if ~exist('FIB_nodes_all','var')
        FIB_nodes_all = fib_nodes;          % [nNodes x 1] for first solve
    else
        FIB_nodes_all(:,end+1) = fib_nodes; % [nNodes x nSolves]
    end
end

    % 
    Vsolve_all(:,:,end+1) = V_atSolve;
    Vdef_all(:,:,end+1)   = V_def;
    solveSucceeded = true;% --- Read last-step element stresses ---

    % Volumes on correct configurations 
    % UZP (unloaded) on V_current
    Vol_LV_UZP = closeAndVolume(F_LV_pressure, V_current);
    Vol_RV_UZP = closeAndVolume(F_RV_pressure, V_current);
    % ED_sim (pressurized) on V_def
    Vol_LV_ED  = closeAndVolume(F_LV_pressure, V_def);
    Vol_RV_ED  = closeAndVolume(F_RV_pressure, V_def);

    % Report vs target ED (from imaged geometry) 
    fprintf('Cavity volumes (µL)  |  UZP → ED_sim   (Δ)\n');
    fprintf('LV: %8.1f → %8.1f  (%+8.1f)\n', Vol_LV_UZP, Vol_LV_ED, Vol_LV_ED-Vol_LV_UZP);
    fprintf('RV: %8.1f → %8.1f  (%+8.1f)\n', Vol_RV_UZP, Vol_RV_ED, Vol_RV_ED-Vol_RV_UZP);

    % Errors vs target ED (in ul from the imaged geometry)
    eLV = abs(Vol_LV_ED - ED_target_LV); % was named errLV_vol
    eRV = abs(Vol_RV_ED - ED_target_RV); % was named errRV_vol

    % Relative Error (fraction of target) 
    relLV = eLV / max(ED_target_LV,1);
    relRV = eRV / max(ED_target_RV,1);

    % Print error: 
    fprintf('ED volume errors vs target → LV: %.1f µL (%.1f%%)  RV: %.1f µL (%.1f%%)\n', ...
            eLV, 100*relLV, eRV, 100*relRV);


        %% === Import principal stresses (s1, s2, s3) ===
    
    % Path to principal stress logfile
    log_principal = fullfile(outDir, febioLogFileName_principalstress);
    
    % Import using GIBBON
    dataP = importFEBio_logfile(log_principal, 0, 1);
    
    P_elem = dataP.data;                    % [nElem x 3 x nSteps]
    [nElemP, nCompP, nStepsP] = size(P_elem);
    
    if nCompP ~= 3
        error('Expected 3 principal stress components (s1,s2,s3), got %d', nCompP);
    end
    
    % === Use principal stresses at FINAL time step ===
    s1 = P_elem(:,1,end);   % max principal
    s2 = P_elem(:,2,end);   % mid principal
    s3 = P_elem(:,3,end);   % min principal
    
    fprintf('Loaded principal stresses: s1 range = [%.3f , %.3f]\n', min(s1), max(s1));
    
    %% === Map element principal stresses → nodes ===
    % Use simple element→node averaging
    % === Map element → node ===
    nNodes = size(V_def,1);
    s1_nodes = elemToNodeMean(E, s1, nNodes);
    s2_nodes = elemToNodeMean(E, s2, nNodes);
    s3_nodes = elemToNodeMean(E, s3, nNodes);
    
    % === Smooth the nodal fields (Laplacian smoothing) ===
    s1_nodes = lapSmoothScalarField(E, s1_nodes, 5, 0.5);
    s2_nodes = lapSmoothScalarField(E, s2_nodes, 5, 0.5);
    s3_nodes = lapSmoothScalarField(E, s3_nodes, 5, 0.5);    
   

    % Assigns node masks for each RV and LV 
    nNodes = size(V_current,1);
    isLV = false(nNodes,1);  isLV(unique(F_LV_pressure(:))) = true;
    isRV = false(nNodes,1);  isRV(unique(F_RV_pressure(:))) = true;


    %  Using scalar weights per node 
    gamma = 0.20;                          % 11/11 original was 0.75 aggressiveness (0.5–0.8 is safe) 11/14 0.25
    w = ones(nNodes,1);
    w(isLV) = w(isLV) .* (1 + gamma * (eLV/max(ED_target_LV,1)));
    w(isRV) = w(isRV) .* (1 + gamma * (eRV/max(ED_target_RV,1)));
    w = max(0.2, min(3.0, w));  


    % Light smoothing of w to avoid sharp updates (3 Jacobi passes)
    A = sparse([],[],[],nNodes,nNodes,6*size(E,1));
    pairs = [E(:,[1 2]);E(:,[1 3]);E(:,[1 4]);E(:,[2 3]);E(:,[2 4]);E(:,[3 4])];
    A = A + sparse(pairs(:,1),pairs(:,2),1,nNodes,nNodes) + sparse(pairs(:,2),pairs(:,1),1,nNodes,nNodes);
    deg = sum(A,2); DinvA = spdiags(1./max(deg,1),0,nNodes,nNodes) * A;
    for it = 1:3, w = 0.5*w + 0.5*(DinvA*w); end
   
    % Inverse update: pull unloaded mesh by the *applied* pressurization 
    beta   = betaInit;
    step     = V_def - V_original;
    %step = V_def - V_current;
    stepW  = step .* w;                    % nodewise scaling of the vector field
    V_trial = V_current - beta * stepW;

% Keep tets positive (shrink step if needed)-
    while ~tetOK(V_trial, E) && beta > 1e-3
        beta    = 0.5 * beta;
        V_trial = V_current - beta * stepW;
    end
    if ~tetOK(V_trial, E)
        warning('Could not find a valid inverse step; stopping.');
        break;
    end
    V_update            = V_trial;       % accept the (possibly shrunk) step

    % --- Store history and accept update ---
    if ~exist('Vdef_all','var'),   Vdef_all   = V_def;   else, Vdef_all(:,:,end+1)   = V_def;   end
    if ~exist('Vsolve_all','var'), Vsolve_all = V_current; else, Vsolve_all(:,:,end+1) = V_current; end
    if ~exist('Vol_LV_ED_hist','var'), Vol_LV_ED_hist = Vol_LV_ED; else, Vol_LV_ED_hist(end+1,1) = Vol_LV_ED; end
    if ~exist('Vol_RV_ED_hist','var'), Vol_RV_ED_hist = Vol_RV_ED; else, Vol_RV_ED_hist(end+1,1) = Vol_RV_ED; end

    V_all(:,:,end+1) = V_update;   % generic history stack
    V_current        = V_update;   % advance unloaded geometry

    % --- Termination criteria ---
    if (relLV < volTol_pct) && (relRV < volTol_pct)
        fprintf('Both ventricles withing %.1f%% of target ED Volume. Stopping .\n', 100*volTol_pct);
        break;
    end
    % k >= maxOuter
    if j >= maxOuter
        fprintf('Reached max inverse iterations (%d). Stopping.\n', maxOuter);
        break;
    end

end
    if exist('V_atSolve_last','var') && ~isempty(V_atSolve_last)
        V_ref_for_plots = V_atSolve_last;   % true pre-solve reference
    else
        V_ref_for_plots = V_current;        % fallback
    end
    if ~solveSucceeded || isempty(V_def)
        warning('No deformed state (V_def) to plot. Skipping this visualization.');
    else
        % ... your gpatch(Fb, V_def, 'r','k',0.25) etc ...
    end

%%
% === Save geometries for external postprocessing (UZP & ED_sim) ===
% This creates e.g. FeBio\Y340W12_geomStates.mat

geomFile = fullfile(savePath, [ratname '_geomStates.mat']);

% MRI ED geometry from imaging
V_MRI_ED = V_original;

% Final unloaded (UZP) geometry used for last successful solve
if exist('V_ref_for_plots','var') && ~isempty(V_ref_for_plots)
    V_UZP_final = V_ref_for_plots;   % the geometry you actually solved from
else
    V_UZP_final = V_current;         % fallback
end

% Final pressurized ED_sim geometry
if exist('V_def','var') && ~isempty(V_def)
    V_ED_final = V_def;
else
    warning('V_def is empty – cannot save ED_sim geometry cleanly.');
    V_ED_final = [];
end

save(geomFile, ...
     'V_MRI_ED', ...      % imaged ED geometry
     'V_UZP_final', ...   % final zero-pressure estimate
     'V_ED_final', ...    % final ED_sim
     'E','Fb','Cb','f0','n0','s0');

fprintf('Saved geometries for postprocessing to:\n  %s\n', geomFile);



%% ========================================================================
%  THESIS-FRIENDLY POSTPROCESSING: GEOMETRIES, STRESSES, ANIMATIONS
%  - Displacement UZP → ED_sim (exaggerated)
%  - Displacement UZP → MRI ED
%  - ED_sim stresses: s1 and von Mises (static)
%  - Animations: UZP → ED_sim (von Mises, fiber-direction σ_f)
%  - LV/RV ED volume convergence across inverse iterations
% ========================================================================

fprintf('\n=== Thesis-friendly postprocessing figures & animations ===\n');

%% ------------------------------------------------------------------------
%  0) Resolve key geometries: MRI ED, final UZP, final ED_sim
% -------------------------------------------------------------------------

% MRI end-diastolic geometry (from imaging)
if ~exist('V_original','var') || isempty(V_original)
    error('V_original (MRI ED geometry) is missing.');
end
V_MRI_ED = V_original;

% Final unloaded (UZP) geometry used for last solve
if exist('V_ref_for_plots','var') && ~isempty(V_ref_for_plots)
    V_UZP_final = V_ref_for_plots;
elseif exist('V_current','var') && ~isempty(V_current)
    V_UZP_final = V_current;
else
    error('No final unloaded geometry (V_ref_for_plots or V_current) found.');
end

% Final pressurized ED_sim geometry
if exist('V_def_last','var') && ~isempty(V_def_last)
    V_ED_sim = V_def_last;
elseif exist('V_def','var') && ~isempty(V_def)
    V_ED_sim = V_def;
else
    warning('No final pressurized geometry (V_def/V_def_last) found. Skipping stress-related plots.');
    V_ED_sim = [];
end

if ~exist('Fb','var') || isempty(Fb)
    error('Fb (surface faces) is missing.');
end

%% ------------------------------------------------------------------------
%  1) Displacement from UZP → ED_sim (WITH visual exaggeration)
% -------------------------------------------------------------------------

if ~isempty(V_ED_sim)
    dispScale = 1;   % <-- try 5–12 for good visualization

    V_ED_sim_exag = V_UZP_final + dispScale*(V_ED_sim - V_UZP_final);

    U_uzp_to_ed    = V_ED_sim - V_UZP_final;
    U_mag_uzp_ed   = sqrt(sum(U_uzp_to_ed.^2,2));

    figure; hold on;
    title([ratname, ': UZP \rightarrow ED_{sim} displacement magnitude (exaggerated)'], 'FontSize', 16);

    hp = gpatch(Fb, V_ED_sim_exag, U_mag_uzp_ed, 'k', 1);   % EXAGGERATED vertices
    hp.FaceColor  = 'interp';
    hp.Marker     = '.';
    hp.MarkerSize = 4;

    axisGeom(gca,16);
    colormap(gjet(250));
    cb = colorbar;
    ylabel(cb, '|U| (UZP \rightarrow ED_{sim}) (mesh units)');
    caxis([0 max(U_mag_uzp_ed)]);
    axis(axisLim(V_ED_sim_exag));
    camlight headlight;
    drawnow;
else
    warning('Skipping UZP → ED_sim displacement plot: V_ED_sim is empty.');
end

%% ------------------------------------------------------------------------
%  2) Displacement from final UZP → MRI ED geometry
% -------------------------------------------------------------------------

U_uzp_to_mri  = V_MRI_ED - V_UZP_final;
U_mag_uzp_mri = sqrt(sum(U_uzp_to_mri.^2,2));

figure; hold on;
title([ratname, ': UZP \rightarrow MRI~ED displacement magnitude'], 'FontSize', 16);

hp = gpatch(Fb, V_MRI_ED, U_mag_uzp_mri, 'k', 1);
hp.FaceColor  = 'interp';
hp.Marker     = '.';
hp.MarkerSize = 4;

axisGeom(gca,16);
colormap(gjet(250));
cb = colorbar;
ylabel(cb, '|U| (UZP \rightarrow MRI~ED) (mesh units)');
caxis([0 max(U_mag_uzp_mri)]);
axis(axisLim(V_MRI_ED));
camlight headlight;
drawnow;

%% ------------------------------------------------------------------------
%  3) ED_sim stresses: principal (s1) and von Mises + ANIMATIONS
% -------------------------------------------------------------------------

vm_nodes = [];  % will fill below if possible

if ~isempty(V_ED_sim)

    % Re-import stress logs from last run
    outDir = fileparts(febioFebFileName);

    % --- 3a) Principal stresses (s1) static map ---
    try
        log_principal = fullfile(outDir, febioLogFileName_principalstress);
        dataP = importFEBio_logfile(log_principal, 0, 1);
        P_elem = dataP.data;                % [nElem x 3 x nSteps]
        [nElemP, nCompP, nStepsP] = size(P_elem);

        if nCompP ~= 3
            warning('Expected 3 principal stresses in principal-stress log, got %d. Skipping s1 plot.', nCompP);
        else
            s1 = P_elem(:,1,end);           % max principal at final step
            % Map element → node
            nNodes   = size(V_ED_sim,1);
            s1_nodes = elemToNodeMean(E, s1, nNodes);

            figure; hold on;
            title([ratname, ': ED_{sim} colored by s_1 (max principal stress)'], 'FontSize', 16);

            hp = gpatch(Fb, V_ED_sim, s1_nodes, 'k', 1);
            hp.FaceColor  = 'interp';
            hp.Marker     = '.';
            hp.MarkerSize = 4;

            axisGeom(gca,16);
            colormap(gjet(250));
            cb = colorbar;
            ylabel(cb,'Max principal stress s_1 (kPa)');
            caxis([min(vm_nodes) max(vm_nodes)]);

            axis(axisLim(V_ED_sim));
            camlight headlight;
            drawnow;
        end
    catch ME
        warning('Could not plot principal stresses (s1): %s', ME.message);
    end

    % --- 3b) Full stress tensor → von Mises static map ---
    try
        dataStruct_stress = importFEBio_logfile(fullfile(outDir, febioLogFileName_stress), 0, 1);
        N_stress_mat      = dataStruct_stress.data;      % [nElem x nComp x nSteps]
        S_last            = N_stress_mat(:,:,end);

        [nElem_log, nComp, ~] = size(S_last); 
        nElem  = size(E,1);
        if nElem_log ~= nElem
            error('Mismatch: stress log nElem (%d) != mesh nElem (%d)', nElem_log, nElem);
        end

        sx = S_last(:,1);
        sy = S_last(:,2);
        sz = S_last(:,3);
        if nComp >= 4
            sxy = S_last(:,4);
            syz = S_last(:,5);
            sxz = S_last(:,6);
        else
            sxy = zeros(nElem,1);
            syz = zeros(nElem,1);
            sxz = zeros(nElem,1);
        end

        % von Mises (3D)
        vm_elem = sqrt( ...
            0.5*((sx-sy).^2 + (sy-sz).^2 + (sz-sx).^2) + ...
            3*(sxy.^2 + syz.^2 + sxz.^2) );

        nNodes   = size(V_ED_sim,1);
        vm_nodes = elemToNodeMean(E, vm_elem, nNodes);
        fprintf('vm_nodes range: [%.4f  %.4f]\n', min(vm_nodes), max(vm_nodes));

        figure; hold on;
        title([ratname, ': ED_{sim} colored by von Mises stress'], 'FontSize', 16);

        hp = gpatch(Fb, V_ED_sim, vm_nodes, 'k', 1);
        hp.FaceColor  = 'interp';
        hp.Marker     = '.';
        hp.MarkerSize = 4;

        axisGeom(gca,16);
        colormap(gjet(250));
        cb = colorbar;
        ylabel(cb,'von Mises stress (kPa)');

        caxis([min(vm_nodes) max(vm_nodes)]);


        axis(axisLim(V_ED_sim));
        camlight headlight;
        drawnow;

    catch ME
        warning('Could not compute/plot von Mises stress: %s', ME.message);
    end

else
    warning('Skipping stress plots: V_ED_sim is empty.');
end

%% ------------------------------------------------------------------------
%  3c) Shared geometry interpolation for animations (UZP → ED_sim)
% -------------------------------------------------------------------------

V_anim = [];
if ~isempty(V_ED_sim)
    nNodes   = size(V_UZP_final,1);
    if ~isequal(size(V_UZP_final), size(V_ED_sim))
        warning('Cannot build animation: V_UZP_final and V_ED_sim have different sizes.');
    else
        nFrames   = 30;
        dispScale = 1;                          % purely visual exaggeration
        V_ED_exag = V_UZP_final + dispScale*(V_ED_sim - V_UZP_final);

        V_anim = zeros(nNodes,3,nFrames);
        for i = 1:nFrames
            t = (i-1)/(nFrames-1);              % 0 → 1
            V_anim(:,:,i) = (1-t)*V_UZP_final + t*V_ED_exag;
        end
    end
end

%% ------------------------------------------------------------------------
%  3d) Animation 1: UZP → ED_sim (von Mises wall stress)
% -------------------------------------------------------------------------

if ~isempty(V_anim) && ~isempty(vm_nodes)
    vm_lo = prctile(vm_nodes,5);
    vm_hi = prctile(vm_nodes,95);

    axStatic = axisLim(V_anim(:,:,end));

    hf_vm = figure; hold on;
    gtitle([febioFebFileNamePart, ': UZP \rightarrow ED_{sim} (von Mises stress)']);

    hp_vm = gpatch(Fb, V_anim(:,:,1), vm_nodes, 'k', 1);
    hp_vm.FaceColor  = 'interp';
    hp_vm.Marker     = '.';
    hp_vm.MarkerSize = markerSize/10;

    axisGeom(gca,fontSize);
    colormap(gjet(250));
    hc = colorbar;
    ylabel(hc,'von Mises stress (kPa)');

    caxis([vm_lo vm_hi]);
    axis(axStatic);
    camlight headlight;

    animStruct = struct;
    animStruct.Time = 1:size(V_anim,3);
    for i = 1:size(V_anim,3)
        animStruct.Handles{i} = hp_vm;
        animStruct.Props{i}   = {'Vertices'};
        animStruct.Set{i}     = {V_anim(:,:,i)};
    end

    anim8(hf_vm,animStruct); drawnow;
else
    warning('Skipping von Mises animation: missing V_anim or vm_nodes.');
end

%% ------------------------------------------------------------------------
%  3e) Animation 2: UZP → ED_sim, colored by fiber-direction stress σ_f
% -------------------------------------------------------------------------

if ~isempty(V_anim) && exist('fib_nodes','var') && ~isempty(fib_nodes)

    if numel(fib_nodes) ~= size(V_anim,1)
        warning('fib_nodes length (%d) != nNodes (%d). Skipping fiber-direction animation.', ...
                numel(fib_nodes), size(V_anim,1));
    else
        cmin_fs = min(fib_nodes);
        cmax_fs = max(fib_nodes);
        axStatic = axisLim(V_anim(:,:,end));

        hf_fs = figure; hold on;
        gtitle([febioFebFileNamePart, ': UZP \rightarrow ED_{sim} (colored by \sigma_f)']);

        hp_fs = gpatch(Fb, V_anim(:,:,1), fib_nodes, 'k', 1);
        hp_fs.FaceColor  = 'interp';
        hp_fs.Marker     = '.';
        hp_fs.MarkerSize = markerSize/10;

        axisGeom(gca,fontSize);
        colormap(gjet(250));
        hc = colorbar;
        ylabel(hc,'Fiber-direction Cauchy stress \sigma_f (kPa)');

        caxis([cmin_fs cmax_fs]);
        axis(axStatic);
        camlight headlight;

        animStruct = struct;
        animStruct.Time = 1:size(V_anim,3);

        for i = 1:size(V_anim,3)
            animStruct.Handles{i} = [hp_fs hp_fs];
            animStruct.Props{i}   = {'Vertices','CData'};
            animStruct.Set{i}     = {V_anim(:,:,i), fib_nodes};
        end

        anim8(hf_fs,animStruct); drawnow;
    end
else
    warning('UZP→ED_{sim} fiber-stress animation skipped: missing V_anim or fib_nodes_last.');
end

%% ------------------------------------------------------------------------
%  4) ED volume convergence across inverse iterations (LV & RV)
% -------------------------------------------------------------------------

if exist('Vol_LV_ED_hist','var') && ~isempty(Vol_LV_ED_hist) && ...
   exist('Vol_RV_ED_hist','var') && ~isempty(Vol_RV_ED_hist) && ...
   exist('ED_target_LV','var')   && exist('ED_target_RV','var')

    LV_ED_target = ED_target_LV;
    RV_ED_target = ED_target_RV;

    figure; tiledlayout(2,1);

    nexttile; hold on;
    plot(Vol_LV_ED_hist,'o-','LineWidth',1.5);
    yline(LV_ED_target,'--','LineWidth',1.2);
    title([ratname, ': LV ED volume across inverse iterations'],'FontSize',14);
    ylabel('Volume (µL)');
    xlabel('Inverse iteration #');
    legend({'Simulated ED','Target ED'},'Location','best');
    grid on;

    nexttile; hold on;
    plot(Vol_RV_ED_hist,'o-','LineWidth',1.5);
    yline(RV_ED_target,'--','LineWidth',1.2);
    title([ratname, ': RV ED volume across inverse iterations'],'FontSize',14);
    ylabel('Volume (µL)');
    xlabel('Inverse iteration #');
    legend({'Simulated ED','Target ED'},'Location','best');
    grid on;

    drawnow;

    % Numeric summary in command window
    Vol_LV_compare = [LV_ED_target, Vol_LV_ED_hist(end)];
    Vol_RV_compare = [RV_ED_target, Vol_RV_ED_hist(end)];
    fprintf('LV ED target vs final ED_{sim}: [%.1f  %.1f] µL\n', Vol_LV_compare(1), Vol_LV_compare(2));
    fprintf('RV ED target vs final ED_{sim}: [%.1f  %.1f] µL\n', Vol_RV_compare(1), Vol_RV_compare(2));
else
    warning('Volume history or targets missing – skipping volume convergence figure.');
end

fprintf('=== Postprocessing complete. ===\n');

%%
% Try to load fiber stress log (if you had it set earlier)
log_fiber = fullfile(outDir, 'fiberStress_out.txt');   % ← or your filename
if exist(log_fiber,'file')
    dataF = importFEBio_logfile(log_fiber,0,1);
    fiber_elem_last = dataF.data(:,1,end);   % [nElem x 1]
else
    warning('No fiber stress logfile found. Fiber animation will not run.');
end

%% Save mesh + labels for VTU export (no rerun needed)
postFile = fullfile(savePath, [ratname '_postForVTU.mat']);

save(postFile, ...
    'ratname', ...
    'V_original','V_UZP_final','V_ED_final', ...
    'E','Fb','Cb', ...
    'LV','RV','S', ...
    'F_LV_pressure','F_RV_pressure','F_base_BC', ...
    'f0','s0','n0', ...
    '-v7.3');

fprintf('Saved VTU export inputs to:\n  %s\n', postFile);



% 


%% Functions

function data = readXDMF(filename)


    % Create a structure to store imported HDF5 data
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
    
    data.Groups(1).Groups(1).Datasets(1).Name = "0";
    data.Groups(1).Groups(1).Datasets(1).Value = h5read(filename, "/Function/FiberDirection/0");
    
    data.Groups(1).Groups(2).Datasets(1).Name = "0";
    data.Groups(1).Groups(2).Datasets(1).Value = h5read(filename, "/Function/Normal/0");
    
    data.Groups(1).Groups(3).Datasets(1).Name = "0";
    data.Groups(1).Groups(3).Datasets(1).Value = h5read(filename, "/Function/SheetDirection/0");
    
    data.Groups(2).Groups(1).Datasets(1).Name = "geometry";
    data.Groups(2).Groups(1).Datasets(1).Value = h5read(filename, "/Mesh/Mesh/geometry");
    
    data.Groups(2).Groups(1).Datasets(2).Name = "topology";
    data.Groups(2).Groups(1).Datasets(2).Value = h5read(filename, "/Mesh/Mesh/topology");
    
    data.Groups(3).Groups(1).Datasets(1).Name = "Values";
    data.Groups(3).Groups(1).Datasets(1).Value = h5read(filename, "/MeshTags/Cell tags/Values");
    
    data.Groups(3).Groups(1).Datasets(2).Name = "topology";
    data.Groups(3).Groups(1).Datasets(2).Value = h5read(filename, "/MeshTags/Cell tags/topology");
    
    data.Groups(3).Groups(2).Datasets(1).Name = "Values";
    data.Groups(3).Groups(2).Datasets(1).Value = h5read(filename, "/MeshTags/Facet tags/Values");
    
    data.Groups(3).Groups(2).Datasets(2).Name = "topology";
    data.Groups(3).Groups(2).Datasets(2).Value = h5read(filename, "/MeshTags/Facet tags/topology");


end
function checkFEBioRunStruct(S)
% Checks common fields runMonitorFEBio expects and prints clear messages.
% Call: checkFEBioRunStruct(febioAnalysis)

req = {'febioPath','run_filename','run_output_path','run_output_name'};
opt = {'fileName_plot','runMode','maxLogCheckTime','disp_on','runFlag'};

fprintf('--- Checking FEBio run struct ---\n');

% 1) Show types/classes for all fields
fns = fieldnames(S);
for i=1:numel(fns)
    v = S.(fns{i});
    c = class(v);
    if isstring(v), c = [c ' (string)']; end
    if iscell(v),   c = [c ' (cell)'];   end
    fprintf('  %-18s : %s\n', fns{i}, c);
end

% 2) Require key fields to exist and be char scalar
for i=1:numel(req)
    k = req{i};
    assert(isfield(S,k), 'Missing field: %s', k);
    v = S.(k);
    if isstring(v), v = char(v); S.(k) = v; end
    if iscell(v),   v = char(v{1}); S.(k) = v; end
    assert(ischar(v) && (isrow(v) || isempty(v)), ...
        'Field %s must be a char row vector (got %s)', k, class(S.(k)));
    assert(~isempty(v), 'Field %s is empty', k);
end

% 3) File/folder existence checks
assert(isfile(S.febioPath),      'febioPath not found: %s', S.febioPath);
assert(isfile(S.run_filename),   'run_filename (.feb) not found: %s', S.run_filename);
if ~exist(S.run_output_path,'dir')
    warning('run_output_path does not exist; creating: %s', S.run_output_path);
    mkdir(S.run_output_path);
end

% 4) Optional plot file: if present, must be char; if absent, suggest one
if isfield(S,'fileName_plot') && ~isempty(S.fileName_plot)
    v = S.fileName_plot;
    if isstring(v), v = char(v); S.fileName_plot = v; end
    if iscell(v),   v = char(v{1}); S.fileName_plot = v; end
    assert(ischar(v) && isrow(v), 'fileName_plot must be char row');
else
    % Suggest a default and report it (does not modify S unless you want to)
    defPlot = fullfile(S.run_output_path, [S.run_output_name '__plot.mat']);
    fprintf('  (info) fileName_plot not set; a typical value is:\n      %s\n', defPlot);
end

% 5) Quick write permission test in output path
testFile = fullfile(S.run_output_path, ['.__write_test_' char(java.util.UUID.randomUUID) '.tmp']);
[fid,msg] = fopen(testFile,'w'); 
assert(fid>0, 'No write permission in run_output_path: %s\n  (%s)', S.run_output_path, msg);
fclose(fid); delete(testFile);

fprintf('✅ Struct looks good.\n\n');
end

function [F_out, n_out, fracNeg] = ensureMinusN(F_in, V_ref, U_last, thresh)
% Check face-avg displacement vs normal sign. If <thresh have U·n<0,
% flip winding so normals point outward (myocardium -> blood).
    n = patchNormal(F_in, V_ref);                 % [nF x 3]
    Uf = (U_last(F_in(:,1),:) + U_last(F_in(:,2),:) + U_last(F_in(:,3),:))/3;
    s = sum(Uf.*n,2);
    fracNeg = mean(s < 0);
    if fracNeg < thresh
        F_out = fliplr(F_in);
        n_out = -n;
    else
        F_out = F_in;
        n_out = n;
    end
end

function vol = closeAndVolume(Fsurf, Vnodes)
% Close an open cavity and return positive volume.
    Edge = patchBoundary(Fsurf);
    [Fc, Vc] = triSurfCloseHoles(double(Fsurf), Vnodes, 0.5, double(Edge));
    Fc = patchNormalFix(Fc);
    vol = patchVolume(Fc, Vc);
    if vol < 0, vol = -vol; end
end

function e = rmse(A,B,dimflag)
% Small RMSE helper that matches your earlier call pattern
    D = A - B; D = D(:);
    e = sqrt(mean(D.^2));
end

function visualizeBC(Fb,V,bcSupportList,F_LV_pressure,F_RV_pressure)
fontSize = 15;
markerSize = 40;
hf=figure;
title('Boundary conditions and loads','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','none',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize/2);
hl(2)=gpatch(F_LV_pressure,V,'r','k',1);
hl(3)=gpatch(F_RV_pressure,V,'g','k',1);
patchNormPlot(F_LV_pressure,V);
patchNormPlot(F_RV_pressure,V);
legend(hl,{'Fixed','LV Pressure','RV Pressure'});


axisGeom(gca,fontSize);
camlight headlight;
drawnow;
end

function Vnod = elemToNodeMean(E, elemVals, nNodes)
% Average element scalar values onto nodes (fast, no loops)
E = double(E);
elemVals = double(elemVals(:));
Vsum   = accumarray(E(:), repmat(elemVals,4,1), [nNodes,1], @sum, 0);
Vcount = accumarray(E(:), 1,                  [nNodes,1], @sum, 0);
Vcount(Vcount==0) = 1;
Vnod = Vsum ./ Vcount;
end

function field_s = lapSmoothScalarField(E, field, nIter, alpha)
%LAPSMOOTHSCALARFIELD  Smooth a nodal scalar field using uniform Laplacian.
%
%   field_s = lapSmoothScalarField(E, field, nIter, alpha)
%
%   E      : [nElem x 4] tetrahedron connectivity
%   field  : [nNodes x 1] scalar at nodes
%   nIter  : smoothing iterations (5–10 recommended)
%   alpha  : smoothing weight (0.3–0.7)

    nNodes = numel(field);
    field_s = field(:);

    % Build adjacency from tet elements
    E = double(E);
    pairs = [E(:,[1 2]); E(:,[1 3]); E(:,[1 4]); ...
             E(:,[2 3]); E(:,[2 4]); E(:,[3 4])];
    A = sparse(pairs(:,1), pairs(:,2), 1, nNodes, nNodes);
    A = A + A.';                        % symmetric adjacency
    A = A - diag(diag(A));

    deg = sum(A,2);
    DinvA = spdiags(1./max(deg,1),0,nNodes,nNodes) * A;

    for k = 1:nIter
        field_s = (1-alpha)*field_s + alpha*(DinvA*field_s);
    end
end