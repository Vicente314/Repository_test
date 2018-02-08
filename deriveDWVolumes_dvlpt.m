function [ud, fid2] = deriveDWVolumes( fiddir, zfill)
%Modifier la description du fichier
% The function will read in a fid and procpar
% the b matrices are generated from vnmrj2.2c sems.c
% the raw data and computed DTI parameter maps are written out.
%
%
% inputs
%       fiddir -
%                string containing the locations of the fid files.
%       zfill  -
%                1x2 matrix to zero fill the data
%                defaults to no zero filling.
%                set to [] to use no zero filling
%                and use the SNR weighting value.
%       method -
%                0 - Justin's unweighted (same result as Richard's)
%                1 - Justin's first weighted
%                2 - Justin's second weighted (Best)
%
% outputs
%       fdf files.  Writes out fdf files for both the
%       raw DWI images, and numerous DTI parameter maps.
%
% Routine was designed to be generic enough to solve
% for any number of directions and b-values.
%
% --------------------------------------------
% Matt Budde - BMRL - WUSTL
% 06_0619 - first version
% 06_0623 - corrected units - allow weighting of calculations by SNR.
% Justin Haldar
% 06_0705 - made lots of changes - gave explicit choice for algorithm, gave
% computational speed option, cleaned up various inefficient spots in the
% code.
% 07_0220 - changed to input data from vnmrj (float format) and to use the
% calcluated b-matrix from vnmrj 2.2c, if the bvalss variable exists in the
% procpar

% Notes:
% subfunctions are self-contained
% Tested with 25 direction datasets with multiple b-values.
% not yet tested with single b-value 25 direction data.
%
% ToDo:
% read and parse array parameter to allow full flexibility
%
% setup files directories and read in parameters

%------------------MODIF LUIS------------------
%Tout ce qui est pr�vu pour le 3D est supprim� (ns c'est imagerie 2D)
%Tout ce qui d�pend de write_nifti est supprim� (tjs fix� � 0 de tte fa�on)
%Tout ce qui est relatif � l'existence d'un masque est supprim�
%Tout ce qui est sur "animal", une variable bizarre
%Tout ce qui est sur le spatial filter (variable im_filter)
%la fonction gaussian2D inutile est supprimee

if (nargin < 1 || isempty(fiddir))
    fiddir = pwd;
end

fidfile = strcat(fiddir,filesep,'fid');

if ~exist(fidfile,'file')
    error('%s could not be found, check paths');
end
%% Lecture du fichier procpar et initialisation de parametres
procparfile = strcat(fiddir,filesep,'procpar');
[pathtobrains, fidname] = fileparts(fiddir);

ud = readprocpar(procparfile);
ud.G = [ud.dpe;ud.dro;ud.dsl]';

ud.diffdirs = size(ud.G,1);

if nargin < 2 || isempty(zfill)
    ud.fn = ud.np/2;
    ud.fn1 = ud.nv;
else %Perhaps not good to use zfill as introduces ripples into image; better to upsample later (in ImageJ)
    ud.fn = zfill(1);
    ud.fn1= zfill(2);
    try ud.fn2 = zfill(3);catch; end %if 3D data
end

% sort bmatrices by b-value, data gets sorted later
[ud.sliceorder, ud.orderIndex] = sort(ud.pss);


%% Get k-space data
% Data structure varies depending on the sequence, the strucure has been
% identified for two sequences. If sequence is different, check structure by debugging
% Let's determine sequence
textfile=strcat(fiddir,filesep,'text');
SEQ=textread(textfile,'%s','delimiter','\n');

if strcmp(SEQ,'Diffusion Weighted (dual) Spin-echo Multi-slice Imaging sequence')
    load_fidfile=fastloadj(fidfile);
    reshape_fidfile=reshape(load_fidfile,ud.np/2,ud.ns,ud.nv,ud.diffdirs);
    fid2 = permute(reshape_fidfile,[1 3 2 4]);
    
elseif strcmp(SEQ,'Spin-echo Multi-slice Imaging sequence')
    load_fidfile=fastloadj(fidfile);
    perm_fidfile = permute(load_fidfile,[1 3 2]); %192x(128x31)x72 complex double; donc l'acquisition de la machine, avec les 31 directions
    reshape_fidfile=reshape(perm_fidfile,ud.np/2,ud.diffdirs,ud.nv,ud.ns);
    %reshape_fidfile=reshape(load_fidfile,ud.np/2,ud.diffdirs,ud.nv,ud.ns);
    fid2 = permute(reshape_fidfile,[1 3 4 2]);
else
    error('Sequence is not recognized. Please complete script to take this new sequence in charge');
end

fid2 = fid2(:,:,ud.orderIndex,:); %fid2 est un "4D complex double"

%zfill
%img =   zeros(ud.np/2,ud.nv,size(fid2,4),size(fid2,3));
img = zeros(ud.fn,ud.fn1,size(fid2,3),size(fid2,4));
%% Derive DW images from raw data (LA)

for sliceno = 1:ud.ns
    
    fid = fid2(:,:,sliceno,:);

    % perform fft and save raw data as fdfs
    dtiImages = fftshift(fftshift(fft2(fid,ud.fn,ud.fn1),1),2); %2D data
    
    %dtiImages = abs(dtiImages); %Negative values introduced for no good reason
    %dtiImages = permute(dtiImages,[2 1 3]);
    %img(:,:,sliceno,:)= dtiImages;
    
    img(:,:,sliceno,:)= abs(dtiImages);
end

dirOut_dw = strrep(pathtobrains,'Brains','DWVolumes');
dirOut = [dirOut_dw,filesep,'no_corr',filesep,fidname];

%ATTENTION CHEKCS THE WHOLE SEARCH PATH, WE ADD './'
if ~exist(['./' dirOut],'dir')  
    mkdir(dirOut);
end

dirOut_tensor = [dirOut_dw,filesep,'tensor'];
if ~exist(['./' dirOut_tensor],'dir')  
    mkdir(dirOut_tensor);
end

%% Repositionnement a l'origine (LA)
%Phase-Encode direction

rpe = ud.lpe / ud.fn1 ; %Resolution dans la direction phase-encode
decalage_axsag = ud.ppe/rpe; %Decalage selon l'axe sagital

if decalage_axsag>=0
    decalage_axsag=round(decalage_axsag);
else
    decalage_axsag= round(ud.fn1 + decalage_axsag);
end
rearrang_sag = [decalage_axsag+1:ud.fn1,1:decalage_axsag];


%Read-out direction

%         rro = ud.lro / (ud.np/2) ; %Resolution dans la direction phase-encode
%         decalage_axcor = ud.pro/rro; %Decalage dans la transversale
%         if decalage_axcor>=0
%             decalage_axcor=round(decalage_axcor);
%         else
%             decalage_axcor= round((ud.np/2) + decalage_axcor);
%         end
%         rearrang_cor = [decalage_axcor+1:(ud.np/2),1:decalage_axcor];

%% Sauvegarde
% Pour output un volume 3D pour chaque direction de gradient
tmp=img;

for dirno=1:size(tmp,4)
    
    img = single(tmp(:,rearrang_sag,:,dirno));
    
    tensor_vol(:,:,:,dirno)=img;
    
    vxl_sz1=ud.lro/size(img,1)*10;
    vxl_sz2=ud.lpe/size(img,2)*10;
    vxl_sz3=ud.thk;
    
    
    % DWI volumes
    out_nii=make_nii(img,[vxl_sz1 vxl_sz2 vxl_sz3]);
    
    fileout = [dirOut,filesep,fidname,'_DWI',num2str(dirno-1),'.nii'];
    save_nii(out_nii,fileout);
    
    
end
    % DTI volume
    out_nii =make_nii(tensor_vol, [vxl_sz1 vxl_sz2 vxl_sz3]);
    
    out_nii.hdr.dime.dim(1)=5; %nbre de dimensions
    out_nii.hdr.dime.dim(5)=1; %pas une serie temporelle, correction de make_nii
    out_nii.hdr.dime.dim(6)=31;
    
    fileout = [dirOut_tensor,filesep,fidname,'_DTI.nii'];

    save_nii(out_nii,fileout);


return
end

