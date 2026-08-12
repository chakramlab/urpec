function [fieldsFileName,T] = urpec_v5( config )
% urpec_v5 Generates a proximity-effect-corrected pattern file for EBL
%
% function [  ] = urpec_v5( config )
% To run urpec and make a run file, see the script run_urpec.
%
% v5 (chakramlab fork, 2026-08-11): supersets v4 with the WEIGHTED
% PROJECTED-GRADIENT solver (config.gradSolver) -- the residual is
% back-propagated through the PSF, so dose is bargained BETWEEN pixels
% under one-sided per-tone constraints (clearing must clear, undercut
% must develop but stay below the PMMA ceiling), optionally weighted
% per feature via config.weightFile. gradSolver=0 (default) behaves
% exactly like urpec_v4, which stays frozen as the pre-gradient engine.
% 
% The corrected file is created by deconvolving a point spread function 
% from an input .dxf or .mat pattern file.
%
% urpec_v4 removes all duplicate vertices in the polygons. 
%
% config is an optional struct with the following optional fields:
%
%   file: pattern file for processing. This can either be a .dxf file or a .mat
%   file. If using a .dxf file, only use 2dpolylines. Mixing 2d and 3d polylines causes problems. 
%   If it is a .mat file, the contents of the file should be a struct
%   called polygons. The polygons struct should have at least these fields:
%       p: an nx2 array of coordinates describing the poylgon.
%       layer: an array of numbers specifying the layer of each polygon
%       according to the convention described above.
%   One unit in the pattern file should be one micron
%   The layer scheme is as follows:
%       The names for all layers should be numbers.
%       Layers 1 and 2 of the input file will
%       both be output to layer 1 of the output file. Layer 1 will not be
%       fractured, and layer 2 will be fractured. Layers 3 and 4 of the input
%       file will be output to layer 2 of the output filed, etc. If the polygons
%       are not fractured, the are written with an average dose.
%
%
%   psfFile: point-spread function file
%
%   dx: spacing in units of microns for the deconvolution. The default is
%   0.01 mcirons. It is also best to have the step size several times
%   larger than the center-center or line-line spacing in npgs. 
%
%   maxIter: maximum number of iterations for the deconvolution. Default is
%   6.
%
%   dvals: doses corresponding to the layers in the output file. Default is
%   1, 1.1, 2.0 in units of the dose to clear.
%   
%   targetPoints: approximate number of points for the simulation. Default
%   is 50e6.
%
%   autoRes: enables auto adjustment of the resolution for ~10min
%   computation time
%
%   fracNum: maximum number of times to divide each polygon during every
%   fracturing iteration.
%
%   fracSize: minimum size for fractured shapes, in units of dx.
%
%   padLen: the size with which to pad the CAD file, in units of microns.
%   The defaul is 5 microns. In general, a larger pad size will improve the
%   accuracy by accounting for long-distance proximity effects, but the
%   computation will take longer.
%
%   outputDir: the directory in which to save the output files.
%
%   savedxf: boolean variable indicating whether or not to save output dxf.
%   Default is false, although autocad is a nice way to view complex files.
%
%   savedc2: boolean variable indicating whether or not to save output dc2.
%   Default is false. If you are using NPGS, set this to true.
%
%   savedose: boolean variable indicating whether or not to save a text file with doses. 
%   Default is false. If you are using NPGS, set this to true.
%
%   npgs: boolean variable indicating whether or not you intend to use
%   NPGS. Default is false. If true, savedc2 and savedose will be set to
%   true.
%
%   overlap: boolean variable indicating how overlaps are handled. If false,
%   urpec will accont for the fact that overlap areas are multiply exposed.
%   If true, urpec will allow over exposure in overlap regions. Default is
%   true.
%
%   triangulate: boolean variable indicating whether or not to triangulate
%   non-convex polygons. Enabling this generates lots of triangles, but all of
%   the polygons will be good, and the fracturing is faster. Default is false.
%
%   convexify: boolean variable indicating whether or not to fracture
%   concave polygons into convex polygons. Default is true.
%
% Call this function without any arguments, or via
% urpec(struct('dx',0.005,'maxIter',6,'dvals',[1:.2:2.4]))
% for example
%
%
% By:
% Adina Ripin 
% Elliot Connors econnors@ur.rochester.edu
% John Nichol jnich10@ur.rochester.edu
%
% Version history
% v2: handles different write fields and writes directly to dc2 format.
% v3: 
%   writes all doses to the same layer but with different colors. 
%   PSF improvements
%   Entirely new fracturing algorithm
% v4: 
%   Code refactoring: new functions shrinkArray and fracturePoly
%   Major speedups in the exposure map creation by using poly2mask
%   Improved fracturing algorithm that will triangulate if needed.

debug=0;

tic

orig_state=warning;
warning('off');

fprintf('urpec is running.\n');

if ~exist('config','var')
    config=struct();
end

config = def(config,'dx',.01);   %Grid spacing in microns. This is can be affected by config.targetPoints and config.autores
config = def(config,'targetPoints',50e6);  %Target number of points for the simulation. 
config = def(config,'autoRes',true);  %auto adjust the resolution
config = def(config,'maxIter',6);  %max number of iterations for the deconvolution
config = def(config,'dvals',linspace(1,2.0,15));  %doses corresponding to different output doses, in units of dose to clear
% layerTargets (chakramlab fork): optional vector of RELATIVE target doses
% per output layer group (group g = input layers 2g-1 and 2g), in units of
% dose to clear. E.g. [1.0 0.4] targets full clearing for group 1 and 40%
% of clearing dose for group 2 (bilayer undercut regions). Empty (default)
% reproduces the original single-target behavior exactly. When using this,
% extend config.dvals downward to cover the lowest target (e.g.
% linspace(0.35,2.0,15)). Overlaps between different-target shapes resolve
% to the HIGHER target.
config = def(config,'layerTargets',[]);
% demandMapWindow (chakramlab fork): [x0 x1 y0 y1] in um -- export the
% CONTINUOUS raw demand map (negatives included, single precision) over
% this window into fields.demandMap/-X/-Y. Empty (default) = no map.
config = def(config,'demandMapWindow',[]);

% nonNeg (chakramlab fork): projected iteration -- clamp doseNew at 0
% after every update, plain residual everywhere. Updates are pointwise
% (Van Cittert), so pixels with unreachable targets (undercut flanks
% blanketed by a neighbor's spillover) simply pin at the floor without
% dragging anything else; every reachable pixel is driven exactly to
% its target, which is also the MINIMAL writable dose there. (A
% band/zero-error variant existed briefly on 2026-08-11 and was
% removed: it froze feasible pixels at their initialization value.)
% relErr=1 additionally scales each pixel's update by 1/target, capped
% at gain 1.8 (fixed-point stability needs gain < 2), so sub-max tones
% descend at comparable FRACTIONAL rates -- the log-loss idea
% translated to this iteration.
config = def(config,'nonNeg',0);
config = def(config,'relErr',0);

% gradSolver (v5): weighted projected-gradient refinement AFTER the
% (warm-start) iteration above converges. Penalties are one-sided per
% tone -- see the gradient loop for the exact forms. bandCeil is the
% PMMA-safety ceiling for undercut/outside pixels; 0.80 is a PLACEHOLDER
% until developer contrast is measured (todo.txt 2026-08-11) -- it sits
% below every clearing threshold of a D0 260-440 sweep at 340 clearing.
% etaGrad: effective per-pixel gain is 2*etaGrad at worst; keep < 1 for
% margin (halve it if the violation printout oscillates).
% weightFile: .mat with polys (cell of Nx2 um) + w (vector); pixels
% outside every weight polygon get w = 1.
config = def(config,'gradSolver',0);
config = def(config,'gradIters',60);
config = def(config,'etaGrad',0.8);
config = def(config,'overBeta',0.25);
config = def(config,'bandCeil',0.80);
config = def(config,'outW',1.0);
config = def(config,'weightFile','');
% bandSlope: gentle constant pull toward the target for undercut pixels
% sitting ABOVE it but inside the band -- the minimal-dose principle.
% Without it the penalty is flat between target and bandCeil and the
% solver leaves free dose in place (the v1 do-nothing-zone lesson,
% 2026-08-11, now applied to the gradient objective).
config = def(config,'bandSlope',0.1);
% clearBand: contrast-derived DEADBAND on the clearing tone -- absorbed
% within clearBand of the target counts as cleared, because finite
% developer contrast clears slightly-under doses too. delta =
% 1 - 10^(-t_scum/gamma): 0.021 at gamma 5.4 (Rooks, interpolated to
% 6.0 C) with 5% tolerable scum. The corridor gets the same bandSlope
% drain (no flat zones -- lesson learned twice). NOTE the corridor and
% develop-T drift spend the same slack budget; the drift reserve is
% the operating-column overdose (todo.txt, deadband entry).
config = def(config,'clearBand',0.021);
% useGPU (v5): run the deconvolution + gradient loops on a CUDA device
% (gpuArray operands make fft2/ifft2/max execute there transparently).
% Needs Parallel Computing Toolbox + a card with ~8 GB free for our
% 50e6-point grids (16 GB comfortable); falls back to CPU with a
% warning if unavailable. Everything is gathered back before
% fracturing/exports, so results are numerically equivalent to CPU up
% to FFT reduction order (harmless at 1%-rung quantization). Per-iter
% diagnostic figures are drawn only on the LAST iteration in GPU mode
% (each draw would round-trip the array over PCIe).
config = def(config,'useGPU',0);
% edgeFree (v5, experiment 2026-08-12): erode the clearing tone's
% penalty mask by edgeFree pixels -- the drawn boundary ring becomes
% DON'T-CARE in the gradient objective. Rationale: holding absorbed >=
% threshold at the last drawn pixel through the PSF forces ~2x written
% edge spikes (the M-profile burrs) whose outward lobes are the width
% overflow. The ring's written dose still serves the interior's
% clearing (no parking pathology: its dose is pulled by neighbors'
% penalties, not its own). Trade: drawn width no longer equals target
% width -- erosion is a designed pre-bias.
config = def(config,'edgeFree',0);
% gradTol (2026-08-12): ADAPTIVE STOPPING for the gradient loop. Stop
% when the total violation changes by less than gradTol (relative)
% over a 10-iteration window; gradIters demotes to a safety cap.
% 0 = disabled (fixed budget, pre-08-12 behavior). Calibration from
% the logged runs: converged plateaus drift ~0.1-0.2%/10 it, active
% descent runs ~1-5%/10 it -- 1e-3 separates them cleanly (v3a shipped
% 40% short of its fixed point under a fixed 60-iteration budget; the
% deadband made convergence detectable by removing the unreachable
% edge pixels from the residual).
config = def(config,'gradTol',0);
% fileSuffix: appended to every output base name so variant runs can
% share one folder without clobbering each other (e.g. '_nn').
config = def(config,'fileSuffix','');
config=def(config,'file',[]);
config=def(config,'psfFile',[]);
config=def(config,'fracNum',4); %Must be >2, otherwise DIVIDEXY has problems.
config=def(config,'fracSize',2);
config=def(config,'padLen',5);
config=def(config,'savedxf',false);
config=def(config,'savedc2',false);
config=def(config,'savedose',false);
config=def(config,'npgs',false);
config=def(config,'overlap',true);
config=def(config,'triangulate',false);
config=def(config,'convexify',true);


if config.npgs
    config.savedose=true;
    config.savedc2=true;
    config.dvals=linspace(config.dvals(1),config.dvals(end),15); %make sure we have 15 colors.
end

%For NPGS: 15 dose values. These jet-like colors are compatible with NPGS
ctab={[0 0 175] [0 0 255] [0 63 255] [0 127 255] [0 191 255] [15 255 239] [79 255 175] [143 255 111] [207 255 047] [255 223 0] [255 159 0] [255 095 0] [255 31 0] [207 0 0] [143 0 0] };

%For other dose configurations, make a jet-like color map.
ncolors=length(config.dvals);
cc=jet(256);
dc=256/ncolors;
if ncolors~=15
    ctab={};
    for ic=1:ncolors
        ctab{ic}=round([cc((ic-1).*dc+1,:)].*255);
    end
end

% ########## Find pattern file ##########

if isempty(config.file)
    %choose and load file
    fprintf('Select your cad file.\n')
    [filename, pathname]=uigetfile({'*.dxf';'*.mat'});
    [pathname,filename,ext] = fileparts(fullfile(pathname,filename));
else
    [pathname,filename,ext] = fileparts(config.file);
    
end

pathname=[pathname filesep];   %filesep, not '\' -- the cluster runs Linux
filename=[filename ext];

config=def(config,'outputDir',pathname);
if config.outputDir(end)~=filesep && config.outputDir(end)~='/'
    config.outputDir=[config.outputDir filesep];
end

% ########## Load pattern from .dxf ##########
if strmatch(ext,'.dxf')
    [lwpolylines,lwpolylayers]=dxf2coord_20(pathname,filename);
    %These are actually the layer names, not the numbers. 
    %If they do not follow the convention, then default to layer 2 and
    %fracturing.
    for i=1:length(lwpolylayers)
        if isempty(str2num(lwpolylayers{i}))
            lwpolylayers{i}='2';
        end
    end
    layerNum=str2num(cell2mat(lwpolylayers));
    
    %splitting each set of points into its own object
    object_num=max(lwpolylines(:,1)); % number of polygons in CAD file
    objects = cell(1,object_num); % cell array of polygons in CAD file
    len = size(lwpolylines,1); % total number of polylines

    % populating 'objects' (cell array of polygons in CAD file)
    count = 1;
    obj_num_ind_start = 1;
    for obj = 1:len % loop over each polyline
        c = lwpolylines(obj,1); % c = object number
        if c ~= count
            objects{count} = lwpolylines(obj_num_ind_start:(obj-1), 1:3); %ID, x, y
            obj_num_ind_start = obj;
            count = count + 1;
            if count == object_num
                objects{count} = lwpolylines(obj_num_ind_start:(len), 1:3);
            end
        end

    end
    if object_num <= 1
        objects{count} = lwpolylines(:, 1:3);
    end
    
    %If something went wrong with the extraction, default to fracturing
    %everything into layer 1.
    if length(layerNum)~=length(objects)
        layerNum=ones(1,length(objects)).*2;
    end
    
    %Make a polygons struct
    polygons=struct();
    for io=1:length(objects)
        polygons(io).layer=layerNum(io);
        
        p=objects{io}(:,2:3);
       
        polygons(io).p=p;
    end
end

% ########## Load pattern from .mat ##########
if strmatch(ext,'.mat')
    d=load([pathname filename]);     
    if isfield(d,'polygons')
        polygons=d.polygons;
        [polygons.dose]=deal(1); %if there is no dose assigned, set it to 1. 
    elseif isfield(d,'fields')
        polygons=d.fields.polygons;
    else
        polygons=struct();
    end
end
fprintf('CAD file imported.\n');

xv=[];
yv=[];

% ########## Analyze polygons ##########

%Find the approximate size of the polygons and try to clean them up. 
%This does not handle duplicate vertices well. After this point, we keep
%track of polygons as structs with fields x and y. 
progressbar('Analyzing polygons');
%TODO: allow for the case that the input struct has x and y fields.
polygons(1).x=[];
polygons(1).y=[];
counter=1;
for ip=length(polygons):-1:1
    progressbar(counter/length(polygons));
    x=polygons(ip).p(:,1);
    y=polygons(ip).p(:,2);
    
    polygons(ip).x=x;
    polygons(ip).y=y;

%     [x1,y1]=fixPoly(x,y);
%     polygons(ip).p=[x1 y1];   
%     polygons(ip).x=x1;
%     polygons(ip).y=y1;

    polys=fixPoly2(polygons(ip));
    polygons=[polygons(1:ip-1) polys polygons(ip+1:end)];
    
    counter=counter+1;
end

for ip=1:length(polygons)
    sx=max(polygons(ip).x)-min(polygons(ip).x);
    sy=max(polygons(ip).y)-min(polygons(ip).y);
    polygons(ip).polysize=min([sx,sy]);
    xv=[xv; polygons(ip).x];
    yv=[yv; polygons(ip).y];
end

%Check for convexity and triangulate non-convex polygons if needed
polygonsTmp=struct();

for fn = fieldnames(polygons)'
   polygonsTmp.(fn{1}) = [];
end

nc=0; %counter for number of non-convex polygons
for ip = 1:length(polygons)
           
    fracture=~mod(polygons(ip).layer,2);

    isConvex = checkConvex(polygons(ip).x,polygons(ip).y);
    if ~isConvex && fracture 
        nc=nc+1;
       
        if config.triangulate || config.convexify
            parent=polygons(ip);
            
            if config.triangulate
                T=triangulatePoly(parent);
            elseif config.convexify
                T=convexify(parent);
            end
            
            T=checkPolys(T,parent);
            T=fixPolys2(T);
            T=checkPolys(T,parent);
    
            for tt=1:length(T)
                polygonsTmp(end+1)=polygons(ip);
                polygonsTmp(end).p=[T{tt}.x(:) T{tt}.y(:)]; %convert to column vectors.
                polygonsTmp(end).x=T{tt}.x(:); 
                polygonsTmp(end).y=T{tt}.y(:);

            end
        else
            polygonsTmp(end+1)=polygons(ip);
        end
        
    else
        polygonsTmp(end+1)=polygons(ip);
    end
            
end
fprintf('Found %d non-convex polygons to be fractured. \n',nc);

polygons=polygonsTmp;
polygons=polygons(2:end); %We initialized it with an empty element.

if debug
    figure(554); clf; hold on;
    for ip=1:length(polygons)
        plot([polygons(ip).x; polygons(ip).x(1)],[polygons(ip).y; polygons(ip).y(1)]);
    end
end

% ########## Creating exposure map ##########

maxX = max(xv);
maxY = max(yv);
minX = min(xv);
minY = min(yv);

dx = config.dx;

fprintf(['Creating 2D binary grid spanning all polygons (spacing = ', num2str(dx), ').\n']);

%Make the simulation area bigger by 5 microns to account for proximity effects
maxXold=maxX;
minXold=minX;
maxYold=maxY;
minYold=minY;
padSize=ceil(5/dx).*dx;
padPoints=padSize/dx;
maxX=maxXold+padSize;
minX=minXold-padSize;
maxY=maxYold+padSize;
minY=minYold-padSize;

xpold = minXold:dx:maxXold;
ypold = minYold:dx:maxYold;

xp = minX:dx:maxX;
yp = minY:dx:maxY;

%Check to make sure we aren't going to use all the memory.
totPoints=length(xp)*length(yp);
fprintf('There are %2.0e points. \n',totPoints);

%Make sure the grid size is appropriate
if config.autoRes && (totPoints<.8*config.targetPoints || totPoints>1.2*config.targetPoints)
    %The following way of changing the resolution will yield a step size
    %that is a power of 2 different that the originally specified step
    %size.
    expand=ceil(log2(sqrt(totPoints/config.targetPoints)));
    dx=dx*2^(expand);

    %This way of changing the resolution will get as close as possible to
    %the target.
%     expand=sqrt(totPoints/config.targetPoints);
%     dx=dx*(expand);
%     dx=round(dx,3);
%     if dx==0
%         dx=0.001;
%     end

    config.dx=dx; %very important for fracturing.
    fprintf('Resetting the resolution to %3.4f.\n',dx);
    padSize=ceil(5/dx).*dx;
    padPoints=padSize/dx;
    maxX=maxXold+padSize;
    minX=minXold-padSize;
    maxY=maxYold+padSize;
    minY=minYold-padSize;
    xp = minX:dx:maxX;
    yp = minY:dx:maxY;
    xpold = minXold:dx:maxXold;
    ypold = minYold:dx:maxYold;
    
    fprintf('There are now %2.0e points.\n',length(xp)*length(yp));
end

xp = minX:dx:maxX;
yp = minY:dx:maxY;

%Make sure the number of points is odd. This is important for deconvolving
%the psf
addX=0;
if ~mod(length(xp),2)
    addX=1;
    xp=[xp xp(end)+dx];
end

addY=0;
if ~mod(length(yp),2)
    addY=1;
    yp=[yp yp(end)+dx];
end

[XP, YP] = meshgrid(xp, yp);

[mp, np] = size(XP);

totgridpts = length(xp)*length(yp);

%polysbin is a 2d array that is zero except inside the polygons
polysbin = zeros(size(XP));
%targetbin (chakramlab fork): per-pixel RELATIVE target dose when
%config.layerTargets is set; overlaps take the higher target.
targetbin = zeros(size(XP));

progressbar(sprintf('Creating the exposure pattern for %d polygons.',length(polygons)))
for ip=1:length(polygons)
    progressbar(ip/length(polygons));

    [xinds,yinds]=shrinkArray(xp,yp,[polygons(ip).x,polygons(ip).y] );

    x=round((polygons(ip).x-xp(xinds(1)))/dx);
    y=round((polygons(ip).y-yp(yinds(1)))/dx);
    subpoly=poly2mask(x',y',length(yinds),length(xinds));
    polysbin(yinds,xinds)=polysbin(yinds,xinds)+subpoly;

    if ~isempty(config.layerTargets)
        grp=min(ceil(polygons(ip).layer/2),length(config.layerTargets));
        tgt=config.layerTargets(max(grp,1));
        targetbin(yinds,xinds)=max(targetbin(yinds,xinds),subpoly.*tgt);
    end

end

[xpts ypts] = size(polysbin);

% ########## Load PSF ##########

if isempty(config.psfFile)
    fprintf('Select point spread function file.\n')
    load(uigetfile('*PSF*'));
else
    load(config.psfFile);
end

minSize=min([max(XP(:))-min(XP(:)),max(YP(:))-min(YP(:))]);
psfRange=round(min([minSize,20]));
%This will break if psfRange is smaller than dx.
npsf=round(psfRange./dx);
[xpsf ypsf]=meshgrid([-npsf:1:npsf],[-npsf:1:npsf]);
xpsf=xpsf.*dx;
ypsf=ypsf.*dx;
rpsf2=xpsf.^2+ypsf.^2;
rpsf=rpsf2.^(1/2);

if ~isfield(psf,'version')
    psf.version=1;
end

switch psf.version
    case 1 %legacy
        eta=psf.eta; %ratio of total back scattered energy to forward-scattered energy.
        alpha=psf.alpha; %alpha, forward scattering range
        beta=psf.beta; %beta, reverse scattering range
        descr=psf.descr;
        
        %compute the forward and backscattered parts separately to enforce
        %the proper ratios.
        psfForward=1/(1+eta).*(1/(pi*alpha^2).*exp(-rpsf2./alpha.^2));
        psfBackscatter=1/(1+eta).*(eta/(pi*beta^2).*exp(-rpsf2./beta.^2));
        
        psf=psfForward./sum(psfForward(:))+eta*psfBackscatter./sum(psfBackscatter(:));
        
    case 2 %version 2
        descr=psf.descr;
        params=psf.params;
        alpha=10^params(1)*1e-3; %alpha, forward scattering range
        beta=10^params(2)*1e-3; %beta, reverse scattering range
        gamma=10^params(3)*1e-3; %gamma, exponential tail range
        
        r=params(4); %ratio of total back scattered energy to forward-scattered energy.
        eta=r-params(5); %eta
        nu=params(5); %gamma
        
        %compute the forward and backscattered parts separately to enforce
        %the proper ratios.
        psfForward= 1/(pi*(1+eta+nu)).*((1/alpha^2).*exp(-(rpsf2)./(alpha^2)));
        psfBackscatter=1/(pi*(1+eta+nu)).*((eta/beta^2).*exp(-(rpsf2)./(beta^2))+(nu/(24*gamma^2)).*exp(-(rpsf./gamma).^(1/2)));
        psf=psfForward./sum(psfForward(:))+(eta+nu)*psfBackscatter./sum(psfBackscatter(:));      
end

%normalize
psf=psf./sum(psf(:));

%Zero pad to at least 10um x 10 um. This is needed in case the sample is
%very small.

%pad along the first dimension
ypad=size(polysbin,1)-size(psf,1);
if ypad>0
    psf=padarray(psf,[ypad/2,0],0,'both');
    padPoints1=padPoints;
    ypdisp=yp;
elseif ypad<0
    polysbin=padarray(polysbin,[-ypad/2,0],0,'both');
    targetbin=padarray(targetbin,[-ypad/2,0],0,'both'); %chakramlab fork
    padPoints1=padPoints-ypad/2;
    ypdisp=[((0:1:-ypad/2-1)+ypad/2).*dx+yp(1) yp ((1:1:-ypad/2).*dx+yp(end))];
end

padPoints1=round(padPoints1);

%pad along the second dimension
xpad=size(polysbin,2)-size(psf,2);
if xpad>0
    psf=padarray(psf,[0,xpad/2],0,'both');
    padPoints2=padPoints;
    xpdisp=xp;
elseif xpad<0
    polysbin=padarray(polysbin,[0,-xpad/2],0,'both');
    targetbin=padarray(targetbin,[0,-xpad/2],0,'both'); %chakramlab fork
    padPoints2=padPoints-xpad/2;
    xpdisp=[((0:1:-xpad/2-1)+xpad/2).*dx+xp(1) xp ((1:1:-xpad/2).*dx+xp(end))];
end

padPoints2=round(padPoints2);

% ########## Deconvolution ##########

dstart=polysbin;
if ~isempty(config.layerTargets)
    %chakramlab fork: per-group target doses. The deconvolution below
    %solves dose*psf = shape, so putting the per-pixel target in shape
    %yields two-tone (or n-tone) correction. Overlap compensation is not
    %target-aware in this mode (overlaps already resolved to max target).
    shape=targetbin;
elseif config.overlap
    %1 inside shapes and 0 everywhere else. Allow overexposure.
    shape=polysbin>0;
else
    %Compensate for overlaps.
    shape=polysbin;
end

dose=dstart;
doseNew=shape; %Initial guess at dose. Just the dose to clear everywhere.
figure(555); clf; imagesc(xpdisp,ypdisp,polysbin);
set(gca,'YDir','norm');
title('CAD pattern');
caxis([0,max(config.dvals)]);
colormap('jet');
colorbar;
drawnow;

%The deconvolution method is to convolve with the psf, add the difference
%between actual dose and desired dose to the programmed dose, and the
%repeat.
meanDose=0;
progressbar('Deconvolving');
%chakramlab fork: fft2(psf) is constant across iterations -- hoist it.
%Bit-identical to computing it inside the loop; ~20% less FFT work.
psfFFT=fft2(psf);
if config.useGPU
    try
        gd=gpuDevice;
        psfFFT=gpuArray(psfFFT);
        doseNew=gpuArray(doseNew);
        shape=gpuArray(shape);
        %loud on purpose: a silent CPU fallback once produced a run that
        %LOOKED like a GPU benchmark (2026-08-11) -- the path taken must
        %be unmissable in every log
        fprintf('*** GPU PATH ACTIVE: %s, %.1f GB free ***\n',...
            gd.Name, gd.AvailableMemory/2^30);
    catch gpuErr
        warning('useGPU requested but unavailable (%s) -- CPU path.',...
            gpuErr.message);
        config.useGPU=0;
    end
end
for iter=1:config.maxIter
    progressbar(iter/config.maxIter);
    %convolve with the point spread function,
    doseActual=ifft2(fft2(doseNew).*psfFFT);

    doseActual=real(fftshift(doseActual)); 
    %The next line is needed becase we are trying to do FFT shift on an
    %array with an odd number of elements. 
    doseActual(2:end,2:end)=doseActual(1:end-1,1:end-1);
    %chakramlab fork: mask with (shape>0), not shape itself -- shape now
    %carries fractional TARGET values when layerTargets is set, and the
    %update below must compare target against the MASKED actual dose.
    %For binary shape this is identical to the original doseActual.*shape.
    doseShape=doseActual.*(shape>0); %actual dose inside shapes only.
    if config.nonNeg
        %chakramlab fork: convergence reported on the max-target
        %(clearing) tone only -- band-tone error is unreachable by
        %design and would poison the metric.
        meanDose=nanmean(doseShape(shape==max(shape(:))))./max(shape(:));
    else
        meanDose=nanmean(doseShape(:))./mean(shape(:));
    end
    
    showFig = ~config.useGPU || iter==config.maxIter;
    if showFig
        figure(556); clf;
        subplot(1,2,2);
        imagesc(xpdisp,ypdisp,doseActual);
        title(sprintf('Actual dose. Iteration %d',iter));
        set(gca,'YDir','norm');
        caxis([0,max(config.dvals)]);
        colormap('jet');
        colorbar;
    end

    if config.nonNeg
        %chakramlab fork: plain residual + projection onto the writable
        %set (dose >= 0). See the nonNeg config comment.
        err=shape-doseShape;
        if config.relErr
            gain=min(1.2./max(shape,1e-6),1.8);
            doseNew=doseNew+gain.*err;
        else
            doseNew=doseNew+1.2*err;
        end
        doseNew=max(doseNew,0);
    else
        doseNew=doseNew+1.2*(shape-doseShape); %Deonvolution: add the difference between the desired dose and the actual dose to doseShape, defined above
    end
    if showFig
        subplot(1,2,1);
        imagesc(xpdisp,ypdisp,doseNew);
        title(sprintf('Programmed dose. Iteration %d',iter));
        set(gca,'YDir','norm');
        caxis([0,max(config.dvals)]);
        colormap('jet');
        colorbar;
    end

    drawnow;
    
end

if meanDose<.98
    warning('Deconvolution not converged. Consider increasing maxIter.');
end

% ########## v5: weighted projected-gradient refinement ##########
if config.gradSolver
    if isempty(config.layerTargets)
        warning('gradSolver needs layerTargets (tone semantics); skipping.');
    else
        tmax=max(shape(:));
        isClear=(shape==tmax)&(shape>0);
        isBand=(shape>0)&(shape<tmax);
        isOut=(shape==0);
        if config.edgeFree>0
            %peel the boundary ring off the clearing penalty mask; the
            %ring joins no other mask (don't-care -- see config note)
            isClear=imerode(isClear,...
                ones(2*config.edgeFree+1));
        end
        %per-pixel weights: 1 everywhere, overridden by weightFile polys
        %(rasterized on THIS padded grid, mirroring the CAD poly2mask).
        %'like' keeps it on the GPU when useGPU is on.
        wmap=ones(size(doseNew),'like',doseNew);
        if ~isempty(config.weightFile)
            WS=load(config.weightFile);
            for kw=1:numel(WS.w)
                pk=WS.polys{kw};
                pxk=round((pk(:,1)-xpdisp(1))./dx);
                pyk=round((pk(:,2)-ypdisp(1))./dx);
                mk=poly2mask(pxk',pyk',size(wmap,1),size(wmap,2));
                wmap(mk)=WS.w(kw);
            end
            fprintf('gradSolver weights: %d polys, w in [%.3g, %.3g]\n',...
                numel(WS.w),min(WS.w),max(WS.w));
        end
        eta=config.etaGrad./max(wmap(:));
        Vhist=zeros(1,config.gradIters);
        Dhist=zeros(1,config.gradIters);
        progressbar('Gradient refinement');
        for it=1:config.gradIters
            progressbar(it/config.gradIters);
            a=ifft2(fft2(doseNew).*psfFFT);
            a=real(fftshift(a));
            a(2:end,2:end)=a(1:end-1,1:end-1);
            %one-sided penalty derivatives (see config comment):
            % clearing: must reach tmax; overshoot costs overBeta
            % undercut: must reach its target; bandCeil caps it
            % outside:  bandCeil caps it (the width-overflow term)
            dphi=zeros(size(a),'like',a);
            cb=(1-config.clearBand)*tmax;
            dphi(isClear)=-2*max(0,cb-a(isClear)) ...
                +2*config.overBeta*max(0,a(isClear)-tmax) ...
                +config.bandSlope*(a(isClear)>cb);
            dphi(isBand)=-2*max(0,shape(isBand)-a(isBand)) ...
                +2*max(0,a(isBand)-config.bandCeil) ...
                +config.bandSlope*(a(isBand)>shape(isBand));
            dphi(isOut)=2*config.outW*max(0,a(isOut)-config.bandCeil);
            g=ifft2(fft2(wmap.*dphi).*psfFFT);
            g=real(fftshift(g));
            g(2:end,2:end)=g(1:end-1,1:end-1);
            doseNew=max(doseNew-eta.*g,0);
            vC=gather(sum(max(0,cb-a(isClear)).^2,'all'));
            vB=gather(sum(max(0,shape(isBand)-a(isBand)).^2,'all'));
            vX=gather(sum(max(0,a(isBand)-config.bandCeil).^2,'all') ...
                +sum(max(0,a(isOut)-config.bandCeil).^2,'all'));
            Vhist(it)=vC+vB+vX;
            %drain mass: written dose sitting in over-ceiling band/out
            %pixels (the reducible reservoir). Its FLATNESS joins the
            %stop test -- the bandSlope drain is a preference, not a
            %violation, so V alone can flatline mid-drain (seen on the
            %u>=0.28 grid rows, 2026-08-12).
            Dhist(it)=gather(sum(doseNew((a>config.bandCeil)& ...
                (isBand|isOut)),'all'));
            if mod(it,10)==0 || it==config.gradIters
                fprintf(['grad %3d: clearing %10.1f  ucLow %10.1f  ' ...
                    'ceiling %10.1f  drain %10.1f\n'],it,vC,vB,vX,...
                    Dhist(it));
            end
            if config.gradTol>0 && it>=20
                dV=abs(Vhist(it)-Vhist(it-10));
                dD=abs(Dhist(it)-Dhist(it-10));
                if dV<config.gradTol*max(Vhist(it),1) && ...
                        dD<config.gradTol*max(Dhist(it),1)
                    fprintf(['grad: CONVERGED at iter %d ' ...
                        '(dV/V = %.2e, dDrain/D = %.2e over 10 ' ...
                        'iters)\n'],it,dV/max(Vhist(it),1),...
                        dD/max(Dhist(it),1));
                    break;
                end
            end
        end
    end
end

%v5 GPU path ends here: everything downstream (fracturing, exports)
%stays on the CPU. gather() is a no-op cost next to the loops above.
%psfFFT must come home too: the windowed-export convolutions reuse it,
%and a gpuArray operand would silently turn doseAbsMap/doseAbsIdealMap
%into gpuArray handles in the saved fields.mat (caught by Eddie's farm
%run, 2026-08-11 -- scipy could not read the maps).
if config.useGPU
    doseNew=gather(doseNew);
    shape=gather(shape);
    psfFFT=gather(psfFFT);
end

dd=doseNew;
ss=shape;

%Unpad arrays
doseNew=dd(padPoints1+1:end-padPoints1-1*addY,padPoints2+1:end-padPoints2-1*addX);
shape=ss(padPoints1+1:end-padPoints1-1,padPoints2+1:end-padPoints2-1);
mp=size(doseNew,1);
np=size(doseNew,2);

%chakramlab fork: capture the RAW demand (negatives included) BEFORE
%the NaN-ing below discards it. doseNew is exactly 0 outside shapes and
%carries the true (possibly negative) programmed-dose demand inside, so
%the inside-shape mask is shape>0. Exports attached at the save section:
%  demandRawHist*  full-range histogram incl. the negative tail
%  demandMap       optional windowed continuous map (single precision),
%                  config.demandMapWindow = [x0 x1 y0 y1] in um
mmr=min(size(shape,1),size(doseNew,1));
nnr=min(size(shape,2),size(doseNew,2));
rawMaskVals=doseNew(1:mmr,1:nnr);
rawMaskVals=rawMaskVals(shape(1:mmr,1:nnr)>0);
demandRawEdges=(floor(min(rawMaskVals)/0.005)*0.005):0.005:(max(rawMaskVals)+0.01);
demandRawCounts=histcounts(rawMaskVals,demandRawEdges);
clear rawMaskVals;
if ~isempty(config.demandMapWindow)
    doseNewRawS=single(doseNew);   %cropped + attached at the save section
    %targetMap: the PRE-BLUR tone ground truth (targetbin unpadded --
    %exact per-pixel target values, razor-sharp edges straight from the
    %polygons). The fitness checks need it for exact tone masks; the
    %ideal-absorbed map is PSF-blurred at every boundary. With
    %layerTargets unset this is just the binary/overlap shape.
    targetMapS=single(shape);
    %chakramlab fork: written + absorbed maps over the same window.
    %written = per-pixel nearest-rung quantization of the demand, with
    %demand<=0 -> 0 (the zero rung / strip_zero_dose: write nothing) --
    %the idealized chip the writer produces (fracturing approximates it
    %to within one rung). absorbed = written conv psf = the dose the
    %resist actually receives, INCLUDING the spillover outside drawn
    %shapes that the clamped-at-zero correction could not remove.
    writtenPad=interp1(config.dvals,config.dvals,dd,'nearest','extrap');
    writtenPad(dd<=0)=0;
    absPad=ifft2(fft2(writtenPad).*psfFFT);
    absPad=real(fftshift(absPad));
    absPad(2:end,2:end)=absPad(1:end-1,1:end-1);
    writtenMapS=single(writtenPad(padPoints1+1:end-padPoints1-1*addY,padPoints2+1:end-padPoints2-1*addX));
    doseAbsMapS=single(absPad(padPoints1+1:end-padPoints1-1*addY,padPoints2+1:end-padPoints2-1*addX));
    %counterfactual: absorbed dose IF the writer could write the raw
    %demand, negatives included -- the deconvolution's own solution, so
    %this should be flat per tone; actual minus ideal = the clamp cost.
    idealPad=ifft2(fft2(dd).*psfFFT);
    idealPad=real(fftshift(idealPad));
    idealPad(2:end,2:end)=idealPad(1:end-1,1:end-1);
    doseAbsIdealS=single(idealPad(padPoints1+1:end-padPoints1-1*addY,padPoints2+1:end-padPoints2-1*addX));
    clear writtenPad absPad idealPad;
end

%This is needed to not count the dose of places that get zero or NaN dose.
try
    doseNew(doseNew==0)=NaN;
    doseNew(doseNew<0)=NaN;
end
[N,X]=hist(doseNew(:),config.dvals);
figure(557); clf; subplot(2,1,1);
hist(doseNew(:),config.dvals);
xlim([0.8,2.5]);
xlabel('Relative dose');
ylabel('Fractured count');
subplot(2,1,2);
hist(doseNew(:),64);
xlim([0.8,2.5]);
xlabel('Relative dose');
ylabel('Pixel count');

drawnow;
if max(N)==N(end)
    warning('Max dose is the highest dval. Consider changing your dvals.');
end

doseStore=doseNew;

dvals=config.dvals;
nlayers=length(dvals);
dvalsl=dvals-(dvals(2)-dvals(1));
dvalsr=dvals;
layer=[];
figSize=ceil(sqrt(nlayers));
dvalsAct=[];
 
% ########## Fracturing ##########

subField=struct();
[XPold, YPold] = meshgrid(xpold, ypold);
nPolys2Frac=sum(~mod([polygons.layer],2));
nPolysNot2Frac=sum(mod([polygons.layer],2));
progressbar(sprintf('Fracturing/Averaging %d/%d',nPolys2Frac,nPolysNot2Frac));
triCount=0; c1=0; c2=0; f1=0; bc=0;
for ip = 1:length(polygons)
    progressbar(ip/length(polygons));
    
    %keep only the parts of the array we care about.
    [xinds,yinds]=shrinkArray(xpold,ypold,polygons(ip).p);
    
    x=round((polygons(ip).x-xpold(xinds(1)))/dx);
    y=round((polygons(ip).y-ypold(yinds(1)))/dx);
    sp=poly2mask(x',y',length(yinds),length(xinds));
    
    subpoly=0.*XPold;
    subpoly(yinds,xinds)=sp;
    
    shotMap=doseNew.*subpoly;
    shotMap(shotMap==0)=NaN;
    shotMapNew=shotMap(yinds,xinds);
    
    fracture=~mod(polygons(ip).layer,2);
        
    if fracture        
        %First figure out the variation in dose across the polygon
        maxDose=max(shotMapNew(:));
        minDose=min(shotMapNew(:));
        
        %Maximum dose variation inside of one shape. 
        dDose=dvals(2)-dvals(1);
        
        minSize=config.dx*config.fracSize; %Smallest eventual polygon size
        fracNum=config.fracNum; %number of times to fracture each polygon each iteration.
        
        if (maxDose-minDose)<dDose %don't need to fracture
            %the subfield struct holds all of the fractured polygons coming
            %from each of the original polygons.
            subField(ip).poly(1).x=polygons(ip).x;
            subField(ip).poly(1).y=polygons(ip).y;
            [mv,ind]=min(abs(dvals-squeeze(nanmean(shotMapNew(:)))));
            subField(ip).poly(1).dose=ind;
            subField(ip).poly(1).layer=ceil(polygons(ip).layer/2);
        else %Need to fracture. In the future, this should be made into a function if increased complexity is desired.
            [subField(ip).poly,td]=fracturePoly(config,shotMapNew,xinds,yinds,XPold,YPold,polygons(ip),ctab);
            triCount=triCount+td.triCount;
            c1=c1+td.c1;
            c2=c2+td.c2;
            f1=f1+td.f1;
            bc=bc+td.bc;
        end
        
        if debug
            figure(779); clf; hold on;
            for pp=1:length(subField(ip).poly)
                plot(subField(ip).poly(pp).x,subField(ip).poly(pp).y,'Color',ctab{subField(ip).poly(pp).dose}./255)
            end
        end
            
    else %no fracturing        
        subField(ip).poly(1).x=polygons(ip).x;
        subField(ip).poly(1).y=polygons(ip).y;
        [mv,ind]=min(abs(dvals-squeeze(nanmean(shotMapNew(:))))); 
        subField(ip).poly(1).dose=ind;
        subField(ip).poly(1).layer=ceil(polygons(ip).layer/2);       
    end

end

fprintf('triCount: %d, c1: %d, c2: %d, f1: %d, bc: %d \n',triCount,c1,c2,f1,bc);

dvalsAct=dvals;

% ########## Exporting ##########

%save the final files
fields=struct();

outputFileName=[config.outputDir filename(1:end-4) '_' descr config.fileSuffix '.dxf'];

if config.savedxf
    fprintf('Exporting to %s\n',outputFileName);
    try
        FID = dxf_open(outputFileName);
    catch
        fprintf('Either dxf_open isn''t on your path, or you need to close the dxf file, and then type dbcont and press enter. \n')
        keyboard;
        FID = dxf_open(outputFileName);
    end
end
        
polygons=struct();
polygons(1)=[];

%Write in order of objects.
%figure(778); clf; hold on;
%the subfield struct array holds all of the fractured polygons which came
%from the initial polygons.
progressbar('Saving');

for ip=1:length(subField)
    progressbar(ip/length(subField));
    
    for iPoly=1:length(subField(ip).poly)
        if ~isempty(subField(ip).poly(iPoly).x) %Needed because sometimes our fracturing algorithm generates empty polygons.
            i=subField(ip).poly(iPoly).dose;
            X=subField(ip).poly(iPoly).x; 
            Y=subField(ip).poly(iPoly).y;
            Z=X.*0;              
            if config.savedxf
                FID=dxf_set(FID,'Color',ctab{i}./255,'Layer',subField(ip).poly(iPoly).layer);
                dxf_polyline(FID,[X(:); X(1)],[Y(:); Y(1)],[Z(:); Z(1)]); %urpec has no duplicate vertices, but dxf files want the same starting and ending vertex.
            end            
            polygons(end+1).p=[X(:) Y(:)]; %Make sure we have column vectors at the end
            polygons(end).x=X(:);
            polygons(end).y=Y(:);
            polygons(end).color=ctab{i}; 
            polygons(end).layer=subField(ip).poly(iPoly).layer;
            polygons(end).lineType=1;
            polygons(end).dose=subField(ip).poly(iPoly).dose;
        end
        
    end
    
end

if config.savedxf
    dxf_close(FID);
end

if config.savedc2
    fields.cadFile=[filename(1:end-4) '_' descr config.fileSuffix '.dc2']; %used by NPGS
    dc2FileName=[config.outputDir filename(1:end-4) '_' descr config.fileSuffix '.dc2'];
    fprintf('Exporting to %s\n',dc2FileName);
    dc2write(polygons,dc2FileName);
end

%Save doses here
if config.savedose
    fields.doseFile=[filename(1:end-4) '_' descr config.fileSuffix '.txt'];
    doseFileName=[config.outputDir filename(1:end-4) '_' descr config.fileSuffix '.txt'];
    fprintf('Exporting to %s\n',doseFileName);
    fileID = fopen(doseFileName,'w');
    fprintf(fileID,'%3.3f \r\n',dvalsAct);
    fclose(fileID);
end

%Save all of the information in a .mat file for later.
%v5: the full config rides with the results -- a number carries its
%provenance (docs/metrics_and_goodhart.md), and the fitness checks
%need to know WHICH solver produced a fields.mat to grade it fairly.
fields.solverConfig=config;
fields.dvalsAct=dvalsAct;
fields.polygons=polygons;
fields.ctab=ctab;
%chakramlab fork: export the CONTINUOUS pre-quantization demand
%histogram (0.005-wide bins) + the write-nothing pixel count, so ladder
%design (analysis/dvals_autoladder.py) can fit the true distribution
%instead of probe-quantized samples. doseNew here is already NaN for
%pixels outside shapes and pixels wanting dose <= 0; 'shape' marks the
%drawn area, so inside-shape NaNs = the zero-class mass.
demandEdges=0:0.005:(max(doseNew(:),[],'omitnan')+0.01);
fields.demandHistEdges=demandEdges;
fields.demandHistCounts=histcounts(doseNew(:),demandEdges);
%shape/doseNew unpad indices can differ by a pixel -- crop to common
mmz=min(size(shape,1),size(doseNew,1));
nnz2=min(size(shape,2),size(doseNew,2));
fields.demandZeroPx=nnz(shape(1:mmz,1:nnz2)>0 & isnan(doseNew(1:mmz,1:nnz2)));
fields.demandDx=dx;
%raw (negative-preserving) exports captured before the NaN step
fields.demandRawHistEdges=demandRawEdges;
fields.demandRawHistCounts=demandRawCounts;
if ~isempty(config.demandMapWindow)
    w=config.demandMapWindow;
    ixw=find(xpold>=w(1) & xpold<=w(2));
    iyw=find(ypold>=w(3) & ypold<=w(4));
    ixw=ixw(ixw<=size(doseNewRawS,2));
    iyw=iyw(iyw<=size(doseNewRawS,1));
    fields.demandMap=doseNewRawS(iyw,ixw);
    fields.demandMapX=xpold(ixw);
    fields.demandMapY=ypold(iyw);
    fields.writtenMap=writtenMapS(iyw,ixw);
    fields.doseAbsMap=doseAbsMapS(iyw,ixw);
    fields.doseAbsIdealMap=doseAbsIdealS(iyw,ixw);
    %shape's unpad indices can differ from doseNew's by a pixel
    iywT=iyw(iyw<=size(targetMapS,1));
    ixwT=ixw(ixw<=size(targetMapS,2));
    fields.targetMap=targetMapS(iywT,ixwT);
    clear doseNewRawS writtenMapS doseAbsMapS doseAbsIdealS targetMapS;
end

fieldsFileName=[config.outputDir filename(1:end-4) '_' descr config.fileSuffix '_fields.mat'];
fprintf('Exporting to %s\n',fieldsFileName);
try
    save(fieldsFileName,'fields');
catch
    fprintf('Can''t save the fields file. Change your path, and then type dbcont and press enter. \n')
    keyboard;
    save(fieldsFileName,'fields');
end
warning(orig_state);

fprintf('urpec is finished.\n')

T=toc;
fprintf('Elapsed time is %3.3f seconds. \n',T);

end

% Apply a default.
function s=def(s,f,v)
if(~isfield(s,f))
    s=setfield(s,f,v);
end
end

