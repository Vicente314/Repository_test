function deriveTensorandScalarMaps( fiddir,method,ud,fid2,mask_path)
%deriveTensorandScalarMaps
%   Estimates the tensor and produces different scalar maps : FA, RAD, L1,
%   V1
%% Initialisation de variables et de paths
debug=0;
%Boolean to output figures (in tif format now)
if (nargin < 1 || isempty(fiddir))
    fiddir = pwd;
end
if nargin < 3 %Choice of estimation method
    method = 0; %0;
end

fidfile = strcat(fiddir,filesep,'fid');

if ~exist(fidfile,'file')
    error('%s could not be found, check paths');
end

% procparfile = strcat(fiddir,filesep,'procpar');
[pathtobrains, fidname] = fileparts(fiddir);

if (method==0)
    datdirname = [fidname '_vnmrj' '.map'];
else
    datdirname = [fidname '_vnmrj' '_algorithm' num2str(method) '.map'];
end
datdir = fullfile(pathtobrains,datdirname);

%% load du masque adequat
maskOk=1;
try
    mask = load_nii(mask_path);
catch
    maskOk=0;
    fprintf('\t\tNo mask available\n')
end
%% Lecture du workspace avec les donnees du fid et initialisation de parametres

if (size(ud.bvalss,2)>1)
    fprintf('\t\tUsing VnmrJ derived b-matrix...\n');
    [ud.sortedb,ud.bsortorder] = sort(ud.bvalue);
    ud.bmat = zeros(3,3,size(ud.bvalue,2));
    ud.bmat(1,1,:) = ud.bvalpp(ud.bsortorder);
    ud.bmat(1,2,:) = ud.bvalrp(ud.bsortorder);
    ud.bmat(1,3,:) = ud.bvalsp(ud.bsortorder);
    ud.bmat(2,1,:) = ud.bmat(1,2,:);
    ud.bmat(2,2,:) = ud.bvalrr(ud.bsortorder);
    ud.bmat(2,3,:) = ud.bvalrs(ud.bsortorder);
    ud.bmat(3,1,:) = ud.bmat(1,3,:);
    ud.bmat(3,2,:) = ud.bmat(2,3,:);
    ud.bmat(3,3,:) = ud.bvalss(ud.bsortorder);
    ud.diffdirs = size(ud.bvalue,2);
    ud.bvalue = ud.bvalue./1000;
    ud.bmat = ud.bmat./1000;
else
    fprintf('\t\tUsing b-matrix Computed from dsl,dro,dpe & gdiff...\n');
    [ud.sortedb,ud.bsortorder] = sort(ud.G(:,1).^2 + ud.G(:,2).^2 + ud.G(:,3).^2);
    % Units here are the gradients in G/cm, the gyromagnetic ratio in
    % rad*MHz/T, and the timings in seconds.
    %[ud.bmat] = dtiGradientsToBMatrices((diag(ud.gdiff)*ud.G(ud.bsortorder,:)).',2*pi*42.58,0,ud.tdelta,ud.tDELTA)/10;
    [ud.bmat] = dtiGradientsToBMatrices((diag(ud.gdiff)*ud.G(ud.bsortorder,:)).',2*pi*42.58,1.25e-4,ud.tdelta,ud.tDELTA)/10;
    
end


%% read information from fid header
fp=fopen(fidfile,'r','ieee-be');
mainHdr.nblocks = fread( fp, 1, 'int32');
mainHdr.ntraces = fread( fp, 1, 'int32');
mainHdr.np = fread( fp, 1, 'int32');
mainHdr.ebytes = fread( fp, 1, 'int32');
mainHdr.tbytes = fread( fp, 1, 'int32');
mainHdr.bbytes = fread( fp, 1, 'int32');
mainHdr.transf = fread( fp, 1, 'int16');
mainHdr.status = fread( fp, 1, 'int16');
statusbits = dec2bin(mainHdr.status,8);
% check if data is stored as int32 (vnmr) or float (vnmrj)
if str2num(statusbits(5))==1
    precision = 'float';
else
    precision = 'int32';
end

mainHdr.spare1 = fread( fp, 1, 'int32');
mainHdr.file_head = 32;
mainHdr.block_head = 28;
mainHdr.blocksize = mainHdr.ntraces*mainHdr.tbytes + mainHdr.block_head;
ud.mainHdrSize = 32;

fclose(fp); % close fid file

%% create directories to save data
if ~debug
    warning('off')
    if (~mkdir(pathtobrains,datdirname))
        error('Could not create data directory, check permissions');
    end
    %create the raw dwi directories
    for jj = 1:ud.ns
        dirname = num2str(100+jj);
        status = mkdir(datdir,strcat('slc_',dirname(2:3)));
    end
    % create the calculated tensor directories
    status=mkdir(datdir,'fdfs_tensor');
    warning('on')
    fdfdir = fullfile(datdir,'fdfs_tensor');
    
    outlist = {'eigen1','eigen2','eigen3','rad','ra','tr','v1x','v1y','v1z','v2x','v2y','v2z','v3x','v3y','v3z','validity','residual'};
    nii_outlist = {'B0','rad','fa','md','eigen1','residual'};
    
    warning('off')
    
    for jj = 1:size(outlist,2)
        status = mkdir(fdfdir,strcat(outlist{jj},'.dat'));
    end
    warning('on')
end
%% Initialise imgs

img_FA = zeros(ud.fn1,ud.fn,size(fid2,3));
img_eigen1 = zeros(ud.fn1,ud.fn,size(fid2,3));
img_RAD =  zeros(ud.fn1,ud.fn,size(fid2,3));
img_MD =  zeros(ud.fn1,ud.fn,size(fid2,3));
img_B0 =   zeros(ud.fn1,ud.fn,size(fid2,3));
img_residual = zeros(ud.fn1,ud.fn,size(fid2,3));
img_RGB = zeros(ud.fn1,ud.fn,size(fid2,3),3);

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

%%
fprintf('\t\tCalculating tensor...\n')

for sliceno = 1:ud.ns
    fprintf(['\t\t\tSlice ' num2str(sliceno) '...\n'])
    
    fid = fid2(:,:,sliceno,:);
    
    dtiImages = fftshift(fftshift(fft2(fid,ud.fn,ud.fn1),1),2);
    dtiImages = dtiImages(:,rearrang_sag,:);
    
    for k=1:ud.diffdirs
        if ~debug
            
            dirnum = find(ud.sliceorder==ud.pss(sliceno));
            dirname = num2str(100+dirnum);
            
            tmp22 = abs(dtiImages(:,:,k));
            tmp22 = rot90(tmp22,2)';
            
            filenum = num2str(100+k);
            filename = strcat(datdir,filesep,'slc_',dirname(2:3),filesep,'image',filenum(2:3),'.fdf');
            f=fopen(filename,'w','ieee-be');
            fhead = makeFDFheader(ud,dirnum,k);
            fwrite(f,fhead,'uint8');
            fwrite(f,tmp22,'float');
            fclose(f);
            
        end
    end
    
    % sort the data by b-value
    dtiImages = dtiImages(:,:,ud.bsortorder);
    %Find b = 0 directions and how many there are
    %
    %     size( niiDWI.img)
    % Calculate the tensor and write the FDF files for each parameter.
    %fprintf('  -Calculating Tensors\n');
    %% Calcul des tenseurs
    
    %%
    switch method
        case 0
            b0Image = mean(dtiImages(:,:,1:ud.nbzero),3);
            dtiImages = dtiImages(:,:,(ud.nbzero+1):end);
            b0Matrix  = ud.bmat(:,:,1);
            b2Matrices = ud.bmat(:,:,(ud.nbzero+1):end);
            for index1 = 1:ud.diffdirs-ud.nbzero
                b2Matrices(:,:,index1) = b2Matrices(:,:,index1)-b0Matrix;
            end
            [lambdas,eigenVectors,residualImage] = dtiLeastSquares2(b0Image, permute(dtiImages(:,:,:),[3 1 2]), b2Matrices);
        case 1
            b0Image = mean(dtiImages(:,:,1:ud.nbzero),3);
            dtiImages = dtiImages(:,:,(ud.nbzero+1):end);
            b0Matrix  = ud.bmat(:,:,1);
            b2Matrices = ud.bmat(:,:,(ud.nbzero+1):end);
            for index1 = 1:ud.diffdirs-ud.nbzero
                b2Matrices(:,:,index1) = b2Matrices(:,:,index1)-b0Matrix;
            end
            [lambdas, eigenVectors,residualImage] = dtiLeastSquaresW(b0Image, permute(dtiImages(:,:,:),[3 1 2]), b2Matrices);
        case 2 %not updated for multiple nbzero values
            [lambdas, eigenVectors,residualImage] = dtiLeastSquaresW2(permute(dtiImages(:,:,:),[3 1 2]), ud.bmat);
        case 3 %added possibility of using a mask
            if maskOk
                [lambdas, eigenVectors,residualImage] = dtiLeastSquaresW2m(permute(dtiImages(:,:,:),[3 1 2]), ud.bmat,mask.img(:,:,sliceno)');
            else
                [lambdas, eigenVectors,residualImage] = dtiLeastSquaresW2(permute(dtiImages(:,:,:),[3 1 2]), ud.bmat);
            end
    end
    lambdas = lambdas.*(lambdas>=0);
    
    
    %% Writing FDF Files & NifTi files
    if ~debug
        
        %fprintf(['\t\tCreating .map folder...\n']);
        slicenum = sliceno;%find(ud.sliceorder==ud.pss(sliceno)); PQ?? Je sais pas...
        filenum = num2str(100+slicenum);
        
        %FDF Files
        for jj=1:size(outlist,2)
            filename = strcat(fdfdir,filesep,outlist{jj},'.dat',filesep,outlist{jj},'_',filenum(2:3),'.fdf');
            switch outlist{jj}
                case 'eigen1'
                    tmp22 = squeeze(lambdas(1,:,:));
                case 'eigen2'
                    tmp22 = squeeze(lambdas(2,:,:));
                case 'eigen3'
                    tmp22 = squeeze(lambdas(3,:,:));
                case 'rad'
                    tmp22 = squeeze(lambdas(2,:,:) + lambdas(3,:,:))./2;
                case 'ra'
                    % set RA values to zero where lambda3 < 0 or lambda1, 2, or 3 is Inf
                    meandiff = mean(lambdas,1);
                    tmp22 = sqrt((lambdas(1,:,:)-meandiff).^2 + (lambdas(2,:,:)-meandiff).^2 + (lambdas(3,:,:)-meandiff).^2)./(sqrt(3).*meandiff);
                    tmp22(find(~isfinite(meandiff))) = 0;
                    tmp22(find(lambdas(3,:,:)<0)) = 0;
                    tmp22 = squeeze(tmp22);
                case 'tr'
                    tmp22 = squeeze(sum(lambdas,1));
                case 'v1x'
                    tmp22 = squeeze(eigenVectors(1,1,:,:));
                case 'v1y'
                    tmp22 = squeeze(eigenVectors(2,1,:,:));
                case 'v1z'
                    tmp22 = squeeze(eigenVectors(3,1,:,:));
                case 'v2x'
                    tmp22 = squeeze(eigenVectors(1,2,:,:));
                case 'v2y'
                    tmp22 = squeeze(eigenVectors(2,2,:,:));
                case 'v2z'
                    tmp22 = squeeze(eigenVectors(3,2,:,:));
                case 'v3x'
                    tmp22 = squeeze(eigenVectors(1,3,:,:));
                case 'v3y'
                    tmp22 = squeeze(eigenVectors(2,3,:,:));
                case 'v3z'
                    tmp22 = squeeze(eigenVectors(3,3,:,:));
                case 'validity'
                    tmp22 = squeeze((lambdas(3,:,:)>0).*isfinite(sum(lambdas,1)));
                case 'residual'
                    tmp22 = squeeze(residualImage(:,:));
            end
            
            tmp22 = rot90(tmp22,2)';
            
            %flip horizontally and vertically
            %tmp22 = fliplr(tmp22);
            %tmp22 = flipud(tmp22);
            
            f2=fopen(filename,'w','ieee-be');
            fhead = makeFDFheader(ud,slicenum);
            fwrite(f2,fhead,'uint8');
            fwrite(f2,tmp22,'float');
            fclose(f2);
            
            
        end %parameter loop
        
        %NifTi files
        for kk=1:size(nii_outlist,2)
            switch nii_outlist{kk}
                case 'B0'
                    tmp22 = abs(dtiImages(:,:,1));
                    
                    img_B0(:,:,sliceno) =  tmp22(:,:)';
                case 'rad'
                    tmp22 = squeeze(lambdas(2,:,:) + lambdas(3,:,:))./2;
                    img_RAD(:,:,sliceno) =  tmp22(:,:)';
                case 'fa' %fractional anisotropy & RGB
                    meandiff = mean(lambdas,1);
                    tmp22 = sqrt(3*((lambdas(1,:,:)-meandiff).^2 + (lambdas(2,:,:)-meandiff).^2 + (lambdas(3,:,:)-meandiff).^2))...
                        ./sqrt(2*(lambdas(1,:,:).^2 +lambdas(2,:,:).^2 +lambdas(3,:,:).^2));
                    tmp22(~isfinite(meandiff)) = 0;
                    tmp22(lambdas(3,:,:)<0) = 0;
                    tmp22 = squeeze(tmp22);
                    img_FA(:,:,sliceno) =  tmp22(:,:)';
                    
                    
                    % Color-coded FA (RGB Color map)  LA
                    tmp_R = squeeze(eigenVectors(1,1,:,:));
                    tmp_G = squeeze(eigenVectors(2,1,:,:));
                    tmp_B = squeeze(eigenVectors(3,1,:,:));
                    
                    sum_colors = abs(tmp_R) + abs(tmp_G) + abs(tmp_B);
                    w_R = abs(tmp_R) ./ sum_colors;
                    w_G = abs(tmp_G) ./ sum_colors;
                    w_B = abs(tmp_B) ./ sum_colors;
                    
                    img_RGB(:,:,sliceno,1) =255*(img_FA(:,:,sliceno).*w_R');
                    img_RGB(:,:,sliceno,2) =255*(img_FA(:,:,sliceno).*w_G');
                    img_RGB(:,:,sliceno,3) =255*(img_FA(:,:,sliceno).*w_B');
                    
                    
                case 'md' %mean diffusivity
                    tmp22 = squeeze(sum(lambdas,1))/3;
                    img_MD(:,:,sliceno) =  tmp22(:,:)';
                case 'eigen1'
                    tmp22 = squeeze(lambdas(1,:,:));
                    img_eigen1(:,:,sliceno) =  tmp22(:,:)';
                case 'residual'  %LA
                    tmp22 = squeeze(residualImage(:,:));
                    img_residual(:,:,sliceno) =  tmp22(:,:)';
                    
            end
            
        end
    end %debug
end %slice loop

%% Sauvegarde des cartes scalaires sous format NifTi (.nii) (LA, JT)
% Gestion des paths
pathout = strrep(pathtobrains,'Brains','ScMaps');
pathout = [pathout,filesep,'Native'];

dirFA = [pathout,filesep,'FA',filesep ];
fileoutFA = [dirFA,fidname,'_FA.nii'];
dirMD = [pathout,filesep,'MD',filesep];
fileoutMD = [dirMD,fidname,'_MD.nii'];
dirL1 = [pathout,filesep,'L1',filesep];
fileoutL1 = [dirL1,fidname,'_L1.nii'];
dirRAD = [pathout,filesep,'RAD',filesep];
fileoutRAD = [dirRAD,fidname,'_RAD.nii'];
dirB0 = [pathout,filesep,'B0',filesep];
fileoutB0 = [dirB0,fidname,'_B0.nii'];
dirResidual = [pathout,filesep,'Residual',filesep];
fileoutResidual = [dirResidual,fidname,'_residual.nii'];
dirRGB = [pathout,filesep,'RGB',filesep];
fileoutRGB = [dirRGB,fidname,'_RGB.nii'];

save_scalar_map(dirFA,fileoutFA,img_FA,ud);    %FA
save_scalar_map(dirMD,fileoutMD,img_MD,ud);    %MD
save_scalar_map(dirL1,fileoutL1,img_eigen1,ud);%L1
save_scalar_map(dirRAD,fileoutRAD,img_RAD,ud); %RAD
save_scalar_map(dirB0,fileoutB0,img_B0,ud);    %B0
save_scalar_map(dirResidual,fileoutResidual,img_residual,ud); %residual image
save_scalar_map4D(dirRGB,fileoutRGB,img_RGB,ud); %RGB image

end

function save_scalar_map(dirout,fileout,nii_img,ud)

if ~isdir(dirout);
    mkdir(dirout);
end
fprintf(['\t\tSaving ' fileout '...\n'])
out_nii=make_nii(single(nii_img),[ud.lro/size(nii_img,1)*10,ud.lpe/size(nii_img,2)*10, ud.thk]);
save_nii(out_nii,fileout);

end
function save_scalar_map4D(dirout,fileout,nii_img,ud)

if ~isdir(dirout);
    mkdir(dirout);
end
fprintf(['\t\tSaving ' fileout '...\n'])
out_nii=make_nii(single(nii_img),[ud.lro/size(nii_img,1)*10,ud.lpe/size(nii_img,2)*10, ud.thk]);
save_nii(out_nii,fileout);

end
