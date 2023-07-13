function [Error] = iBrainMiniOutput_run(Cfg)
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
iBrainPath = fileparts(which('iBrain.m'));
if ~exist([Cfg.WorkingDir,filesep,'T1Img'],'dir')
    mkdir([Cfg.WorkingDir,filesep,'T1Img'])
end
if ~exist([Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN'],'dir')
    mkdir([Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN'])
end
%% process structural raw data
if strcmp(Cfg.StartingDir,'T1Raw') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) && Cfg.IsNeedConvertT1DCM2NII
    parfor temp_subject=1:SubjectNum
        disp(['Current processing subject: ', Cfg.SubjectID{temp_subject}])
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
        %kick out nan and extreme values
        
        disp('-------------------------------------------------------------------')
        disp('Done processing dcm2nii...')
    end
    Cfg.StartingDir='T1Img';
end
%% process structural R2SN
load([iBrainPath,filesep,'model_data',filesep,'wi_90.mat']);
load([iBrainPath,filesep,'model_data',filesep,'train_data.mat']);
if strcmp(Cfg.StartingDir,'T1Img') && logical(exist([Cfg.WorkingDir,filesep,Cfg.StartingDir],'dir')) 
    parfor temp_subject=1:SubjectNum %use attribute to judge whether process each step
        disp(['Current processing subject: ',Cfg.SubjectID{temp_subject}])
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
        if Cfg.IsReportInIndividualSpace 
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'cat12_MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,'mwp1',Cfg.SubjectID{temp_subject},'.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'cat12_MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,'mwp2',Cfg.SubjectID{temp_subject},'.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
        end
        disp('-------------------------------------------------------------------')
        disp('Done processing segmentation...')
        output_Dir=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep];
        template_Data=Cfg.CoregisterTemplate;
        system(['python ',iBrainPath,filesep,'Scripts',filesep,'run_reg_single.py --input_file=',temp_InputFile,32,'--template_path=',template_Data,32,'--output_dir=',output_Dir]);
        if Cfg.IsReportInIndividualSpace
            system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'],32,'--input_img=',...
                [iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
                [iBrainPath,filesep,'Atlas',filesep,'hippo_mask.nii'],32,'--savePath=',...
                [Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject}]]); 
        end
        disp('-------------------------------------------------------------------')
        disp('Done processing coregistration...')
        temp_InputFile=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'DNR_',Cfg.SubjectID{temp_subject},'.nii'];              
        RF = generate_atlas_features(temp_InputFile,Cfg.T1AtlasMaskFile); 
        parsave([output_Dir,filesep,Cfg.SubjectID{temp_subject},'_AtlasFeatures.mat'],RF);
        disp('-------------------------------------------------------------------')
        disp('Done processing atlas feature extraction...') 
        [MaskFilepath,filename]=fileparts(Cfg.T1AtlasMaskFile);
        GM_atlas_img = load_nii([MaskFilepath,filesep,'cat12_',filename,'.nii']);
        hippo_reference_GM_path=[iBrainPath,filesep,'Report_template',filesep,'hippo_BN246_volume_reference.csv'];
        Atlas_reference_GM_path=[iBrainPath,filesep,'Report_template',filesep,'BN246_volume_reference.csv'];
        temp_InputStruct=dir([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,'mwp1*.nii']);
        temp_GMFile=[Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,'mri',filesep,temp_InputStruct(1).name];
        GM_img = load_nii(temp_GMFile);
        GM_volume = get_GM(GM_img,GM_atlas_img);
        temp_load=load([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'_AtlasFeatures.mat']);
        test_data = constructR2SN_test(temp_load,GM_volume,wi); 
        output_Dir=[Cfg.WorkingDir,filesep,'Results',filesep,'T1ImgR2SN',filesep,Cfg.SubjectID{temp_subject},filesep];
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end
        img_save_name=[output_Dir,'RiskBrainMap.nii'];
        %generate subject result variables for python report, include CAT12 segmentation gm wm result path, hippocampus path
        orig_img_path = [Cfg.WorkingDir,filesep,'T1Img',filesep,Cfg.SubjectID{temp_subject},filesep,Cfg.SubjectID{temp_subject},'.nii'];
        generate_python_T1_report_input(train_data,test_data,GM_img,GM_atlas_img,hippo_reference_GM_path,Atlas_reference_GM_path,Cfg.IsReportInIndividualSpace,orig_img_path,img_save_name);
        %call python to generate word report
        output_Dir=[Cfg.WorkingDir,filesep,'Results',filesep,'Report',filesep,Cfg.SubjectID{temp_subject},filesep];
        if ~exist(output_Dir,'dir')
            mkdir(output_Dir);
        end

        system(['python ',[iBrainPath,filesep,'Scripts',filesep,'Report_combined.py'],32, '--iBrainPath=',...
            iBrainPath,32,'--WorkingDir=',Cfg.WorkingDir,32,'--SubjectID=',Cfg.SubjectID{temp_subject},32,...
            '--ReportWholeBrain=',num2str(Cfg.IsReportWholeBrain),32,'--ReportInIndividualSpace=',num2str(Cfg.IsReportInIndividualSpace),32,'--MiniOutput=',...
            num2str(Cfg.IsMinimumOutput),32,'--savePath=',output_Dir(1:end-1),32,'--writeMode=','T1',32,'--tempFCsession=',num2str(1)]);
        %remove subject T1Img folder
        rmdir([Cfg.WorkingDir,filesep,'T1Img'], 's');
        disp('----------------------------------------------------------------')
        disp('Done processing R2SN...')
    end
end

delete(gcp('nocreate'));
end

