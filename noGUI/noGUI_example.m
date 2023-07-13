%no GUI batch sample
clear all;clc;
%load configure file
miniOutput=0;%specify if needs to minimise output
%set working directory
iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'noGUI',filesep,'iBrain_Run_T1_from_T1Raw_Cfg.mat']);
%Cfg.WorkingDir=['E:',filesep,'MATLABtoolboxes',filesep,'iBrain_test_data',filesep,'T1Only',filesep];%can replace with your own working directory
Cfg.WorkingDir='/data/DATA_SHARE/OASIS/';
Cfg.StartingDir='T1Img';
Cfg.CoregisterTemplate=[iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii'];
Cfg.T1AtlasMaskFile=[iBrainPath,filesep,'Atlas',filesep,'BN_AtlaS_246_1mm.nii'];
Cfg.IsReportWholeBrain=1;
Cfg.IsReportInIndividualSpace=0;
%update subject list
DirFiles=dir([Cfg.WorkingDir,filesep,Cfg.StartingDir,filesep]);
if ismac
    dropDir=3;
else
    dropDir=2;
end
for temp_sub_ind=1:length(DirFiles)-dropDir
    Cfg.SubjectID{temp_sub_ind}=DirFiles(temp_sub_ind+dropDir).name;
end

%run ibrain batch
if miniOutput
    Cfg.ParallelWorkers=10;
    Cfg.IsMinimumOutput=1;
    iBrainMiniOutput_run(Cfg);
else
    Cfg.ParallelWorkers=10;
    Cfg.IsMinimumOutput=0;   
    iBrain_run(Cfg);
end