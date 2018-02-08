function fieldInhomogeneityCorrection(DWFolder,maskFile)
%fieldInhomogeneityCorrection applies ANTS' N4BiasFieldCorrection

%   Detailed explanation goes here
%% Managing input paths and output paths
[pathToFile, brainName] = fileparts(DWFolder);

dirOut = strrep(pathToFile,'DWVolumes', 'N4CorrectedDWVolumes');
dirOut = [dirOut,filesep,brainName];

if ~exist(['./' dirOut],'dir')
    mkdir(dirOut);
end

cd(DWFolder);

%% Computing Bias Field on B0 volume
im_in = [DWFolder,filesep,brainName,'_DWI0.nii'];
im_out = [dirOut,filesep,brainName,'_N4_DWI0.nii'];
bias_field = [dirOut,filesep, brainName, '_biasField.nii'];

commandN4 = ['$HOME/antsbin/bin/N4BiasFieldCorrection -d 3 '...
    '-i ' im_in ' '...
    '-x ' maskFile ' '...'/home/luis/Documents/Maitrise/TBSS Pipeline/Masks/LPSI09_mask.nii'
    '-o [' im_out ',' bias_field '] ' ...
    '-s 2'];

fid = fopen('tempscript.sh','w');
fprintf(fid,'%s',commandN4);
fclose(fid);
%Make script executable
commandR = 'chmod +x tempscript.sh';
[status,cmdout] = system(commandR,'-echo')
%Execute
commandR = './tempscript.sh';
[status,cmdout] = system(commandR,'-echo')


%% Applying modification to diffusion-weigthed volumes

%L'idee etait d'utiliser le biasfield calcule sur B0 et de l'appliquer sur
%le autres volumes mais il s'agit d'operations plus complexes au'une simple
%division (cf notes). Donc en suspens pour l'instant.

% Usage: ImageMath ImageDimension <OutputImage.ext> [operations and inputs] <Image1.ext> <Image2.ext>
listing=dir('*.nii');


for elt=listing'
    
    im_in = [DWFolder,filesep,elt.name];
    im_out=[dirOut,filesep,strrep(elt.name,'_','_N4_')];
   
    commandA = ['$HOME/antsbin/bin/ImageMath 3 '...
        im_out ' '...
        '/ ' im_in ' ' bias_field ' '...
        ];

    fid = fopen('tempscript.sh','w');
    fprintf(fid,'%s',commandA);
    fclose(fid);

    %Execute
    commandR = './tempscript.sh';
    [status,cmdout] = system(commandR,'-echo')
end

delete('tempscript.sh');

end

