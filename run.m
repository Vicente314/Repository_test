%% PREPARATION DES FICHIERS POUR LANCER UNE ETUDE
%
% DOSSIER s_XXX : Tout les resultats y seront stock�s
%     DOSSIER Brains : avec les dossiers .fid
%     [DOSSIER Masks : Pour l'etape 2]
%         [DOSSIER Native  : avec les masques d'origines]
%     DOSSIER Atlas : avec l'atlas et sa segmentation

%clear all;clc;
%HOMEDIR='/Users/ldealmei/Desktop/Edu/EPM/M.ScA. Biomedical/rDTIDA';
HOMEDIR='/Users/ENGVIC00/Desktop/LUIS/Luis_Akakpo_2016/Package_LASL'; %Modification 1 Vicente

cd(HOMEDIR);
%addpath('/Users/ldealmei/Documents/MATLAB/z_nifti')
addpath('/Users/ENGVIC00/Documents/MATLAB/z_nifti')
addpath('/Users/ENGVIC00/Desktop/LUIS/Luis_Akakpo_2016/DTI/rDTIDA_scripts/ROI_Analysis')
version ='1.0';

%-----------------------------------------------------------------
%----------------------------DEBUG--------------------------------
%for debug purposes
P1=0;
P2=0;
P3=1;
P4=0;

%-----------------------------------------------------------------
%for debug purposes : PI
dw_bool=0;
tensor_bool=1;
%-----------------------------------------------------------------
%for debug purposes : PIII
temp_cons_bool=1;
reg_to_temp_bool=0;
reg_to_atlas_bool=0;
seg_scmap_bool=0;
seg_clean_bool=0;
%-----------------------------------------------------------------
params_p1 =[];
info_p1 =[];
params_p2 =[];
info_p2 =[];
params_p3  =[];
info_p3 =[];
params_p4 =[];
info_p4 =[];

prompt='Study?'; %Vicente: test_data??
study_name=input(prompt,'s');

cd(study_name)
start_global=tic;
FLAG=0;
diary on
%mail_notif=0;
%% PHASE I : SCALAR MAPS CREATION
if P1
   
    fprintf('PHASE I : SCALAR MAPS CREATION\n');
    
    start_P1=tic;
    
    params_p1=struct; %structure array of parameters for the phase I
    
    %----------------------------------------------------------------------
    %----------------------------Parameters--------------------------------
    prompt='Please enter zfill. Just press Enter if you do not want to zero-fill the data.';
    zf=input(prompt);
    zfill=[zf,zf];
    
    params_p1.zfill=zf;
    params_p1.n4_corr=0;
    params_p1.ec_corr=0;
    params_p1.tensor_est_method=3; % method 3 = method 2 + masquage
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    maindir='Brains';
    brains = dir(fullfile(maindir,'*.fid'));
    
    fprintf([num2str(length(brains)) ' folder(s) have been found.\n'])
    
    for i=1:size(brains,1)
        
        fprintf(['\tStarting treatement of ', brains(i).name,'\n'])
        
        fiddir = strcat(maindir,filesep,brains(i).name);
        dwdir=strrep(fiddir,'Brains',['DWVolumes' filesep 'no_corr']);
        dwdir=strrep(dwdir,'.fid','');
        if dw_bool
            [ud, fid2]=deriveDWVolumes_dvlpt(fiddir,zfill);
        end
        %         if params_p1.n4_corr
        %             fieldInhomogeneityCorrection(); %TO COME
        %         end
        %         if params_p1.ec_corr
        %             ECandMotionCorrection(); %TO COME
        %         end
        mask_path=strrep(fiddir,'Brains',['Masks' filesep 'Native']);
        mask_path=strrep(mask_path,'.fid','_mask.nii');
        if tensor_bool
            deriveTensorandScalarMaps(fiddir,params_p1.tensor_est_method,ud,fid2,mask_path);
        end
    end
    dur_P1=toc(start_P1)
    
    info_p1=struct; %structure array contening info relative to the execution of phase I
    info_p1.brains=brains;
    info_p1.dur=dur_P1;
end
%% PHASE I' : MASKS CREATION AND CHOICE OF SCMAP TO PROCESS

%----------------------------------------------------------------------
%----------------------------Parameters--------------------------------

scmap_type='skel';

%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% PHASE II : SCALAR MAPS PROCESSING
if P2
    fprintf('PHASE II : SCALAR MAPS PROCESSING\n');
    start_P2=tic;
    
    params_p2=struct; %structure array of parameters for the phase II
    
    %----------------------------------------------------------------------
    %----------------------------Parameters--------------------------------
    params_p2.iso=0;
    params_p2.mask=0;
    params_p2.rescale=0;
    params_p2.orient=1;
    params_p2.interpol_method='cubic';
    params_p2.res=0.072; %Change selon le zfill!
    params_p2.cur_orient=[4 6 5];% in vivo :[4 3 2];
    params_p2.scmap_type={'RAD','B0','FA','MD','L1','RGB'};
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    for i=1:length(params_p2.scmap_type)
        %Other var
        cur_scmap_type=params_p2.scmap_type{i};
        
        scmaps_folder=['ScMaps' filesep 'Native' filesep cur_scmap_type];
        excluded_files='';
        % Run Steps
        fprintf(['\tTreating  ' cur_scmap_type ' images...\n'])
        
        if params_p2.iso
            fprintf(['\tIsotropising maps to a resolution of ' num2str(params_p2.res) ' ...\n'])
            scmaps_pre01_iso(scmaps_folder,params_p2.res,params_p2.interpol_method)
        end
        if params_p2.mask
            fprintf('\tMasking files...\n')
            
            excluded_files=scmaps_pre02_mask(scmaps_folder,['Masks' filesep 'Native']);
        end
        if params_p2.rescale
            fprintf('\tRescaling by a factor of XXX ...\n')
            scmaps_pre03_rescale(scmaps_folder,[1,1,1])
        end
        if params_p2.orient
            fprintf('\tReorienting scalar maps ...\n')
            scmaps_pre04_orient(scmaps_folder,params_p2.cur_orient)
        end
    end
    dur_P2=toc(start_P2)
    
    info_p2=struct; %structure array contening info relative to the execution of phase I
    info_p2.dur=dur_P2;
    info_p2.ex=excluded_files;
    
end
%% PHASE III : SEGMENTATION
if P3
    fprintf('PHASE III : REGISTRATION AND SEGMENTATION\n');
    
    
    params_p3=struct; %structure array of parameters for the phase II
    
    %%  III-1 : Template Construction (TC)
    
    %-----------------------------------------------------------------
    %-------------------------Parameters------------------------------
    temp_cons_params.temp_server=0;
    %     hostname='64.18.68.105';
    %     usrname='lakakpo';
    %     pswd='';
    
    temp_cons_params.regex='*.nii.gz';
    
    temp_cons_params.init_rig=1;                     
    temp_cons_params.metric='MI';
    temp_cons_params.metric_param=32; %number of bins for MI, radius for CC
    temp_cons_params.update_it=3;                   
    
    temp_cons_params.isdefault=1; %ATTENTION A METTRE A ZERO SI MODIF DES PARAMS CI DESSOUS!! SINON SANS EFFET!!
    temp_cons_params.init_rig_it='500x100x25x10';     %/\
    temp_cons_params.init_rig_conv='1e-8';           %/  \
    
    temp_cons_params.rig_it='500x100x25x10';          %||
    temp_cons_params.rig_conv='1e-8';                 %||
    
    temp_cons_params.aff_it='500x100x25x10';          %||
    temp_cons_params.aff_conv='1e-8';                 %||
    
    temp_cons_params.SyN_it='100x50x25x10';           %||
    temp_cons_params.SyN_conv='1e-9';                 %||
    %-----------------------------------------------------------------
    %-----------------------------------------------------------------
    
    params_p3.tc=temp_cons_params;
    
    start_P3_1=tic;
    if temp_cons_bool
        if temp_cons_params.temp_server
            %A ADAPTER UNE FOIS LE SERVEUR DISPO
            %             fprintf('\tConstructing template on server...\n')
            %             templateConstruction_server(['./ScMaps/Processed/' scmap_type],'','../rDTIDA/Scripts/Segmentation/Shell_Scripts');
        else
            fprintf('\tConstructing template locally...\n')
            try
                templateConstruction_local(scmap_type,temp_cons_params);
            catch ME
               % matlabmail('luis.akakpo@gmail.com',ME.message,'La construction du template a echouee');
                FLAG=1;
            end
        end
    end
    dur_P3_1=toc(start_P3_1)
    
    %% III-2 Registration to Template (RT)
    
    template=['.' filesep 'segmentation' filesep 'template' filesep 'template0.nii.gz'];
    
    %-----------------------------------------------------------------
    %-------------------------Parameters------------------------------
    reg_to_temp.regex='LPS*.nii';
    %-----------------------------------------------------------------
    %-----------------------------------------------------------------
    
    reg_to_temp.scmap_type=scmap_type;
    %CES PARAMS DOIVENT ETRE IDENTIQUES SOUS PEINE DE BIAISER LA
    %COMPARAISON
    reg_to_temp.metric=temp_cons_params.metric;
    reg_to_temp.metric_param=temp_cons_params.metric_param; %number of bins for MI, radius for CC
    reg_to_temp.rig_it=temp_cons_params.rig_it;
    reg_to_temp.rig_conv=temp_cons_params.rig_conv;
    reg_to_temp.aff_it=temp_cons_params.aff_it;
    reg_to_temp.aff_conv=temp_cons_params.aff_conv;
    reg_to_temp.SyN_it=temp_cons_params.SyN_it;
    reg_to_temp.SyN_conv=temp_cons_params.SyN_conv;
    
    params_p3.rt=reg_to_temp;
    
    start_P3_2=tic;
    if reg_to_temp_bool
        fprintf(['\tRegistering ' reg_to_temp.regex ' files to template ' template ' ...\n'])
        try
            register_to_template(template,reg_to_temp);
        catch ME
            %matlabmail('luis.akakpo@gmail.com',ME.message,'La registration sur le template a echouee');
            FLAG=1;
        end
    end
    
    dur_P3_2=toc(start_P3_2)
    %% III-3 Registration of Template to Atlas (RTA)
    
    template=['.' filesep 'segmentation' filesep 'template' filesep 'template0.nii.gz'];
    
    %-----------------------------------------------------------------
    %-------------------------Parameters------------------------------
    atlas=['.' filesep 'Atlas' filesep 'p24_average_fa_iso07.nii.gz'];
    
    reg_to_atlas.metric='MI';
    reg_to_atlas.metric_param=32;
    reg_to_atlas.rig_it='500x100x25x10';
    reg_to_atlas.rig_conv='1e-8';
    
    reg_to_atlas.aff_it='500x100x25x10';
    reg_to_atlas.aff_conv='1e-8';
    
    reg_to_atlas.SyN_it='100x50x25x10';
    reg_to_atlas.SyN_conv='1e-9';
    %-----------------------------------------------------------------
    %-----------------------------------------------------------------
    
    params_p3.rta=reg_to_atlas;
    
    start_P3_3=tic;
    if reg_to_atlas_bool
        fprintf(['\tRegistering template' template 'to atlas ' atlas ' ...\n'])
        try
            register_template_to_atlas(atlas,template,reg_to_atlas);
        catch ME
            %matlabmail('luis.akakpo@gmail.com',ME.message,'La registration du template sur l atlas a echouee');
            FLAG=1;
        end
        
    end
    dur_P3_3=toc(start_P3_3)
    
    %% III-4 SEGMENTATION IN NATIVE SPACE (SNS)
    start_P3_4=tic;
    if seg_scmap_bool
        fprintf('\tSegmenting scalar maps...\n')
        try
            segment_scmaps();
        catch ME
            %matlabmail('luis.akakpo@gmail.com',ME.message,'La segmentation des scmaps a echouee');
            FLAG=1;
        end
        
    end
    dur_P3_4=toc(start_P3_4)
    
    if seg_clean_bool
        fprintf('\tCleaning Segmentation...\n')
        clean_seg(scmap_type);
    end
    
    %% info structuration
    info_p3=struct; %structure array contening info relative to the execution of phase III
    %this array is divided in 4 other array structures, one for each
    %substep of phase III
    
    info_p3_tc=struct;
    info_p3_tc.dur=dur_P3_1;
    
    scmaps_dir=['ScMaps' filesep 'Processed' filesep scmap_type ];
    list_cons = dir(fullfile(scmaps_dir,temp_cons_params.regex));
    info_p3_tc.list_cons=list_cons; %list of files used for template construction
    
    info_p3_rt=struct;
    info_p3_rt.dur=dur_P3_2;
    
    list_reg = dir(fullfile(scmaps_dir,reg_to_temp.regex));
    info_p3_rt.list_reg=list_reg; %list of files left to register to template
    
    info_p3_rta=struct;
    info_p3_rta.dur=dur_P3_3;
    info_p3_rta.atlas=atlas;
    %     info_p3_rta.atlas_dim=atlas;
    %     info_p3_rta.atlas_voxdim=atlas;
    %     info_p3_rta.atlas_lbs=atlas;
    info_p3_rta.template=template;
    % info_p3_rta.template_dim=template;
    % info_p3_rta.template_voxdim=template;
    
    info_p3_sns=struct;
    info_p3_sns.dur=dur_P3_4;
    
    
    info_p3.tc=info_p3_tc;
    info_p3.rt=info_p3_rt;
    info_p3.rta=info_p3_rta;
    info_p3.sns=info_p3_sns;
    
end
%% PHASE 4 : MORPHOMETRY
if P4
    fprintf('PHASE IV : MORPHOMETRY\n');
    
    params_p4=struct; %reste vide pr l'instant mais on sait jms
    
    start_P4=tic;
    csvname=['segmentation',filesep,study_name '_morph.csv'];
    try
        output_morph_data(csvname,scmap_type);
    catch ME
        %matlabmail('luis.akakpo@gmail.com',ME.message,'L output des ROI a echoue');
        FLAG=1;
    end
    dur_P4=toc(start_P4)
    info_p4=struct;
    info_p4.dur=dur_P4;
    info_p4.filename=csvname;
end



%% FIN
dur_global=toc(start_global)
%if ~FLAG
  %  if mail_notif
        %matlabmail('luis.akakpo@gmail.com',['Felicitations :) Duree : ' num2str(dur_global) ],'EXECUTION REUSSIE!');
        
 %   end
%end
diary off

%% Write study info file
output_study_info(version,study_name, params_p1, info_p1,params_p2, info_p2, params_p3, info_p3, params_p4, info_p4 )

%MODIFICATION TEST IN GITHUB BY VICENTE
