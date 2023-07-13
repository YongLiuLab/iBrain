function [Error] = iBrain_run(Cfg)
%Generate shell file using Cfg struct 
delete(gcp('nocreate'));
if Cfg.ParallelWorkers==0
    Cfg.ParallelWorkers=1;
end
parpool('local',Cfg.ParallelWorkers);
if exist('gcp.m','file')
    try
        gcp;
    end
elseif parpool('size')==0
    try
        parpool;
    end
end

Error=[];
defaults=spm_get_defaults;
SubjectNum=length(Cfg.SubjectID);
disp(num2str(SubjectNum));
iBrainPath = fileparts(which('iBrain.m'));
%% process structural data
if strcmp(Cfg.StartingDir,'T1Raw') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) && Cfg.IsNeedConvertT1DCM2NII
    if ~exist([Cfg.WorkingDir,filesep,'T1Img'],'dir')
        mkdir([Cfg.WorkingDir,filesep,'T1Img'])
    end
    parfor temp_subject=1:SubjectNum
        temp_Dir=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep];
        output_Dir=[Cfg.WorkingDir,filesep,'T1Img',filesep,Cfg.SubjectID{temp_subject},filesep]; 
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end
        if ~exist([output_Dir,Cfg.SubjectID{temp_subject},'.nii'],'file')
            system([iBrainPath,filesep,'dependenices',filesep,'dcm2nii',filesep,'dcm2niix  -f ',Cfg.SubjectID{temp_subject}, ' -o ',output_Dir,32,temp_Dir]);
        else
            disp(['Already exist nifti image in T1Img path, skip dicom2nii transform for subject: ',Cfg.SubjectID{temp_subject}])
        end
        disp('-------------------------------------------------------------------')
        disp(['Done processing dcm2nii for subject: ',Cfg.SubjectID{temp_subject}])
    end
    Cfg.StartingDir='T1Img';
    disp('-------------------------------------------')
    disp('Done processing dcm2nii for all subjects...')
end


if strcmp(Cfg.StartingDir,'T1Img') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) && Cfg.IsSegment 
    suffix='S';% prefix added to the name of files after segment
    if ~exist([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix],'dir')
        mkdir([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix])
    end
    parfor temp_subject=1:SubjectNum       
        temp_InputStruct=dir([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'*.nii*']);
        if length(temp_InputStruct)>1
            Error=[Error;['More than one image file within subject ', Cfg.SubjectID{temp_subject},' directory when running segmentation']];
        else
            if ~strcmp(temp_InputStruct(1).name,[Cfg.SubjectID{temp_subject},'.nii']) %force subject file name match folder name
                if strcmp(temp_InputStruct(1).name(end-2:end),'.gz')
                    unzipsuccess = gunzip([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,temp_InputStruct(1).name]);
                    if unzipsuccess
                        delete([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,temp_InputStruct(1).name]);
                        movefile([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,temp_InputStruct(1).name(1:end-3)], [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii']);
                    end                    
                else
                    movefile([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,temp_InputStruct(1).name], [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii']);
                end               
            end
        end
        temp_InputFile=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'];
        run_spm_cat(temp_InputFile); 
        output_Dir=[Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix,filesep,Cfg.SubjectID{temp_subject},filesep];          
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end
        if isunix
            system(['mv -f ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'label',32,output_Dir]);
            system(['mv -f ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'report',32,output_Dir]);
            system(['mv -f ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',32,output_Dir]);
            system(['mv -f ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'surf',32,output_Dir]);
        elseif ispc
            system(['move /y ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'label',32,output_Dir,'label']);
            system(['move /y ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'report',32,output_Dir,'report']);
            system(['move /y ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',32,output_Dir,'mri']);
            system(['move /y ',Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'surf',32,output_Dir,'surf']);
        end
        if Cfg.IsReportInIndividualSpace 
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'cat12_MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [output_Dir,'mri',filesep,'mwp1',Cfg.SubjectID{temp_subject},'.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'cat12_MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [output_Dir,'mri',filesep,'mwp2',Cfg.SubjectID{temp_subject},'.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
        end
        disp('------------------------------------------------------------------------')
        disp(['Done processing segmentation for subject: ',Cfg.SubjectID{temp_subject}])
    end
    disp('------------------------------------------------')
    disp('Done processing segmentation for all subjects...')
end
if ~isempty(Error)
    disp(Error);
    return;
end
if strcmp(Cfg.StartingDir,'T1Img') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) && Cfg.IsCoregister
    suffix='C';% prefix added to the name of files after coregister
    if ~exist([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix],'dir')
        mkdir([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix])
    end
    template_Data=Cfg.CoregisterTemplate;
    parfor temp_subject=1:SubjectNum
        temp_InputStruct=dir([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii']);
        if length(temp_InputStruct)>1
            Error=[Error;['More than one image file within subject: ', Cfg.SubjectID{temp_subject},' directory， when running coregister']];
        end
        temp_InputFile=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,temp_InputStruct(1).name];
        output_Dir=[Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix,filesep,Cfg.SubjectID{temp_subject},filesep];
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end       
        system(['python ',iBrainPath,filesep,'Scripts',filesep,'run_reg_single.py --input_file=',temp_InputFile,32,'--template_path=',template_Data,32,'--output_dir=',output_Dir]);
        if Cfg.IsReportInIndividualSpace
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [iBrainPath,filesep,'Atlas',filesep,'hippo_mask.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
        end
        disp('--------------------------------------------------------------------------')
        disp(['Done processing coregistartion for subject: ',Cfg.SubjectID{temp_subject}])
    end
    Cfg.StartingDir='T1ImgC';
    disp('--------------------------------------------------')
    disp('Done processing coregistration for all subjects...')
end
if ~isempty(Error)
    disp(Error);
    return;
end
if strcmp(Cfg.StartingDir,'T1ImgC') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) && Cfg.IsExtractT1AtlasFeature 
    suffix='A';% prefix added to the name of files after coregister
    if ~exist([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix],'dir')
        mkdir([Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix])
    end
    parfor temp_subject=1:SubjectNum
        temp_InputFile=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'DNR_',Cfg.SubjectID{temp_subject},'.nii'];              
        RF = generate_atlas_features(temp_InputFile,Cfg.T1AtlasMaskFile);        
        output_Dir=[Cfg.WorkingDir,filesep,Cfg.StartingDir,suffix,filesep,Cfg.SubjectID{temp_subject},filesep];        
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end          
        parsave([output_Dir,filesep,Cfg.SubjectID{temp_subject},'_AtlasFeatures.mat'],RF);
        disp('------------------------------------------------------------------------------------------')
        disp(['Done processing atlas based feature extraction for subject: ',Cfg.SubjectID{temp_subject}])
    end 
    Cfg.StartingDir=[Cfg.StartingDir,'A'];
    disp('------------------------------------------------------------------')
    disp('Done processing atlas based feature extraction for all subjects...')   
end


if strcmp(Cfg.StartingDir,'T1ImgCA') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir'))  && logical(exist([Cfg.WorkingDir,filesep,'T1ImgS'],'dir')) && Cfg.IsR2SN 
    if ~exist([Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN'],'dir')
        mkdir([Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN'])
    end
    [MaskFilepath,filename]=fileparts(Cfg.T1AtlasMaskFile);
    GM_atlas_img = load_nii([MaskFilepath,filesep,'cat12_',filename,'.nii']);
    hippo_reference_GM_path=[iBrainPath,filesep,'Report_template',filesep,'hippo_BN246_volume_reference.csv'];
    Atlas_reference_GM_path=[iBrainPath,filesep,'Report_template',filesep,'BN246_volume_reference.csv'];
    load([iBrainPath,filesep,'model_data',filesep,'wi_90.mat']); 
    load([iBrainPath,filesep,'model_data',filesep,'train_data.mat']);
    parfor temp_subject=1:SubjectNum
        temp_InputStruct=dir([Cfg.WorkingDir,filesep,'T1ImgS',filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,'mwp1*.nii']);
        if length(temp_InputStruct)>1
            Error=[Error;['More than one image file within subject ', Cfg.SubjectID{temp_subject},' directory when running segmentation']];
        end
        temp_GMFile=[Cfg.WorkingDir,filesep,'T1ImgS',filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,temp_InputStruct(1).name];
        GM_img = load_nii(temp_GMFile);
        GM_volume = get_GM(GM_img,GM_atlas_img);
        temp_load=load([Cfg.WorkingDir,filesep,'T1ImgCA',filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'_AtlasFeatures.mat']);
        test_data = constructR2SN(temp_load,GM_volume,wi);        
        output_Dir=[Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN',filesep,Cfg.SubjectID{temp_subject},filesep];
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end
        img_save_name=[output_Dir,filesep,'RiskBrainMap.nii'];
        %generate subject result variables for python report, include CAT12 segmentation gm wm result path, hippocampus path
        orig_img_path = [Cfg.WorkingDir,filesep,'T1Img',filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'];
        generate_python_T1_report_input(train_data,test_data,GM_img,GM_atlas_img,hippo_reference_GM_path,Atlas_reference_GM_path,Cfg.IsReportInIndividualSpace,orig_img_path,img_save_name);
        %call python to generate PDF report  
        output_Dir=[Cfg.WorkingDir,filesep,'Results',filesep,'Report',filesep,Cfg.SubjectID{temp_subject},filesep];
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end
        
        system(['python ',[iBrainPath,filesep,'Scripts',filesep,'Report_combined.py'],32, '--iBrainPath=',...
            iBrainPath,32,'--WorkingDir=',Cfg.WorkingDir,32,'--SubjectID=',Cfg.SubjectID{temp_subject},32,...
            '--ReportWholeBrain=',num2str(Cfg.IsReportWholeBrain),32,'--ReportInIndividualSpace=',num2str(Cfg.IsReportInIndividualSpace),32,'--MiniOutput=',...
            num2str(Cfg.IsMinimumOutput),32,'--savePath=',output_Dir(1:end-1),32,'--writeMode=','T1',32,'--tempFCsession=',num2str(1)]);

        disp('----------------------------------------------------------------')
        disp(['Done processing R2SN for subject: ',Cfg.SubjectID{temp_subject}])
    end
    disp('----------------------------------------------------')
    disp('Done processing R2SN calculation for all subjects...')
end
if ~isempty(Error)
    disp(Error);
    return;
end
delete(gcp('nocreate'));
end

