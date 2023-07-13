function varargout = iBrain(varargin)
% IBRAIN MATLAB code for iBrain.fig
%      IBRAIN, by itself, creates a new IBRAIN or raises the existing
%      singleton*.
%
%      H = IBRAIN returns the handle to a new IBRAIN or the handle to
%      the existing singleton*.
%
%      IBRAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IBRAIN.M with the given input arguments.
%
%      IBRAIN('Property','Value',...) creates a new IBRAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iBrain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iBrain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iBrain

% Last Modified by GUIDE v2.5 26-May-2023 20:10:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iBrain_OpeningFcn, ...
                   'gui_OutputFcn',  @iBrain_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before iBrain is made visible.
function iBrain_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for iBrain
handles.output = hObject;
iBrainPath=fileparts(which('iBrain.m'));
% add default configurations
handles.Cfg.WorkingDir=pwd;
handles.Cfg.DataProcessDir =handles.Cfg.WorkingDir;
handles.Cfg.SubjectID={};
handles.Cfg.IsOneKeyRunFromT1Raw=0;
handles.Cfg.IsOneKeyRunFromT1Img=0;
handles.Cfg.IsNeedConvertT1DCM2NII=0;
handles.Cfg.IsSegment=0;
handles.Cfg.IsCoregister=0;
handles.Cfg.CoregisterTemplate=[iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii'];
handles.Cfg.IsExtractT1AtlasFeature=0;
handles.Cfg.T1AtlasMaskFile=[iBrainPath,filesep,'Atlas',filesep,'BN_Atlas_246_1mm.nii'];
handles.Cfg.IsR2SN=0;

handles.Cfg.IsReportWholeBrain=0;
handles.Cfg.IsReportInIndividualSpace=0;
handles.Cfg.IsMinimumOutput=0;
handles.Cfg.ParallelWorkers=0;
handles.Cfg.StartingDir='';

axes(handles.bgAxes);
img = imread('icon_in_GUI_with_text.jpg');
imshow(img);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes iBrain wait for user response (see UIRESUME)
% uiwait(handles.figiBrainMain);


% --- Outputs from this function are returned to the command line.
function varargout = iBrain_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in btnSelectWorkingDir.
function btnSelectWorkingDir_Callback(hObject, eventdata, handles)
theDir =handles.Cfg.WorkingDir;
theDir =uigetdir(theDir, 'Please select the Working directory: ');
if ~isequal(theDir, 0)
    SetWorkingDir(hObject,handles, theDir);
end

   
function SetWorkingDir(hObject, handles, ADir)
if 7==exist(ADir,'dir')
    handles.Cfg.WorkingDir =ADir;
    handles.Cfg.DataProcessDir =handles.Cfg.WorkingDir;
end
guidata(hObject, handles);
UpdateDisplay(handles);

% --- Executes during object creation, after setting all properties.
function edtWorkingDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtWorkingDir_Callback(hObject, eventdata, handles)
theDir =get(hObject, 'String');	
SetWorkingDir(hObject,handles, theDir);


% --- Executes during object creation, after setting all properties.
function listSubjectID_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listSubjectID.
function listSubjectID_Callback(hObject, eventdata, handles)
theIndex =get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on key press with focus on listSubjectID and none of its controls.
function listSubjectID_KeyPressFcn(hObject, eventdata, handles)
%Delete the selected item when 'Del' is pressed
key =get(handles.figiBrainMain, 'currentkey');
if seqmatch({key},{'delete', 'backspace'})
    DeleteSelectedSubjectID(hObject, eventdata,handles);
end

function DeleteSelectedSubjectID(hObject, eventdata, handles)	
theIndex =get(handles.listSubjectID, 'Value');
if size(handles.Cfg.SubjectID, 1)==0 || theIndex>size(handles.Cfg.SubjectID, 1)
    return;
end
theSubject =handles.Cfg.SubjectID{theIndex, 1};
tmpMsg=sprintf('Delete the Participant: "%s" ?', theSubject);
if strcmp(questdlg(tmpMsg, 'Delete confirmation'), 'Yes')
    if theIndex>1
        set(handles.listSubjectID, 'Value', theIndex-1);
    end
    handles.Cfg.SubjectID(theIndex, :)=[];
    if size(handles.Cfg.SubjectID, 1)==0
        handles.Cfg.SubjectID={};
    end
    guidata(hObject, handles);
    UpdateDisplay(handles);
end

function DeleteAllSubjects(hObject, eventdata, handles)	
tmpMsg=sprintf('Delete all the participants?');
if strcmp(questdlg(tmpMsg, 'Delete confirmation'), 'Yes')
    handles.Cfg.SubjectID={};
    guidata(hObject, handles);
    UpdateDisplay(handles);
end

function LoadSubIDFromTextFile(hObject, eventdata, handles)
[SubID_Name , SubID_Path]=uigetfile({'*.txt','Subject ID Files (*.txt)';'*.*', 'All Files (*.*)';}, ...
    'Pick the text file for all the subject IDs');
SubID_File=[SubID_Path,SubID_Name];
if ischar(SubID_File)
    if exist(SubID_File,'file')==2
        fid = fopen(SubID_File);
        IDCell = textscan(fid,'%s\n'); %YAN Chao-Gan. For compatiblity of MALLAB 2014b. IDCell = textscan(fid,'%s','\n');
        fclose(fid);
        handles.Cfg.SubjectID=IDCell{1};
        guidata(hObject, handles);
        UpdateDisplay(handles);
    end
end
    
function SaveSubIDToTextFile(hObject, eventdata, handles)
[SubID_Name , SubID_Path]=uiputfile({'*.txt','Subject ID Files (*.txt)';'*.*', 'All Files (*.*)';}, ...
    'Specify a text file to save all the subject IDs');
SubID_File=[SubID_Path,SubID_Name];
if ischar(SubID_File)
    fid = fopen(SubID_File,'w');
    for iSub=1:length(handles.Cfg.SubjectID)
        fprintf(fid,'%s\n',handles.Cfg.SubjectID{iSub});
    end
    fclose(fid);
end


function ReLoadSubjects(hObject, eventdata, handles)	
handles.Cfg.SubjectID={};
guidata(hObject, handles);
UpdateDisplay(handles);
    
% --------------------------------------------------------------------
function onekeyRunT1_ButtonDownFcn(hObject, eventdata, handles)
setT1OnekeyPipelineDisplayStatus(hObject,handles,'on');
setT1StepbystepDisplayStatus(hObject,handles,'off');

% --- Executes on button press in ckboxonekeyFromT1Raw.
function ckboxonekeyFromT1Raw_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    setT1OnekeyPipelineDisplayStatus(hObject,handles,'on');
    handles.Cfg.IsOneKeyRunFromT1Raw = 1;
    handles.Cfg.IsOneKeyRunFromT1Img = 0;
else
    handles.Cfg.IsOneKeyRunFromT1Raw = 0;
    handles.Cfg.IsOneKeyRunFromT1Img = 1;
end
setT1StepbystepDisplayStatus(hObject,handles,'off')
handles.Cfg.IsNeedConvertT1DCM2NII=1;
handles.Cfg.StartingDir='T1Raw';
setT1StepbystepDeafultValue(hObject, handles);

% --- Executes on button press in ckboxonekeyFromT1Img.
function ckboxonekeyFromT1Img_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    setT1OnekeyPipelineDisplayStatus(hObject,handles,'on');
    handles.Cfg.IsOneKeyRunFromT1Img = 1;
    handles.Cfg.IsOneKeyRunFromT1Raw = 0;
else
    handles.Cfg.IsOneKeyRunFromT1Img = 0;
    handles.Cfg.IsOneKeyRunFromT1Raw = 1;
end
setT1StepbystepDisplayStatus(hObject,handles,'off');
handles.Cfg.IsNeedConvertT1DCM2NII=0;
handles.Cfg.StartingDir='T1Img';
setT1StepbystepDeafultValue(hObject, handles);

function setT1StepbystepDisplayStatus(hObject, handles, status)
set(handles.ckboxT1DICOM2NIFTI,'Enable',status);
set(handles.ckboxSegmentation, 'Enable', status);
set(handles.ckboxCoregister, 'Enable', status);
set(handles.ckboxAtlasFeatureExtraction, 'Enable', status);
% set(handles.btnSelectAtlasFile,'Enable', status);
% set(handles.edtAtlasFile ,'Enable', status);
set(handles.ckboxR2SN, 'Enable', status);
guidata(hObject, handles);
UpdateDisplay(handles);


function setT1StepbystepDeafultValue(hObject, handles)
[iBrainPath, fileN, extn] = fileparts(which('iBrain.m'));
handles.Cfg.IsSegment=1;
handles.Cfg.IsCoregister=1;
handles.Cfg.CoregisterTemplate=[iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii.gz'];
handles.Cfg.IsExtractT1AtlasFeature=1;
handles.Cfg.T1AtlasMaskFile=[iBrainPath,filesep,'Atlas',filesep,'BN_Atlas_246_1mm.nii'];
handles.Cfg.IsR2SN=1;
handles.Cfg.ParallelWorkers=1;
guidata(hObject, handles);
UpdateDisplay(handles);


function setT1OnekeyPipelineDisplayStatus(hObject,handles,status)
set(handles.ckboxonekeyFromT1Raw,'Enable',status);
set(handles.ckboxonekeyFromT1Img,'Enable',status);
guidata(hObject,handles);
UpdateDisplay(handles);


% --------------------------------------------------------------------
function stepbystepRunT1_ButtonDownFcn(hObject, eventdata, handles)
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
setT1StepbystepDisplayStatus(hObject,handles,'on')

% --- Executes on button press in ckboxT1DICOM2NIFTI.
function ckboxT1DICOM2NIFTI_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsNeedConvertT1DCM2NII = 1;
else
    handles.Cfg.IsNeedConvertT1DCM2NII = 0;
end
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
guidata(hObject, handles);
UpdateDisplay(handles);

% --- Executes on button press in ckboxSegmentation.
function ckboxSegmentation_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsSegment = 1;
else
    handles.Cfg.IsSegment = 0;
end
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
guidata(hObject, handles);
UpdateDisplay(handles);

% --- Executes on button press in ckboxCoregister.
function ckboxCoregister_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsCoregister = 1;
else
    handles.Cfg.IsCoregister = 0;
end
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
guidata(hObject, handles);
UpdateDisplay(handles);


% --- Executes on button press in ckboxAtlasFeatureExtraction.
function ckboxAtlasFeatureExtraction_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsExtractT1AtlasFeature = 1;    
else
    handles.Cfg.IsExtractT1AtlasFeature = 0;
end
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
guidata(hObject, handles);
UpdateDisplay(handles);

% --- Executes during object creation, after setting all properties.
function edtAtlasFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function edtAtlasFile_Callback(hObject, eventdata, handles)
% theMaskfile =get(hObject, 'String');
% theMaskfile =strtrim(theMaskfile);
% if exist(theMaskfile, 'file')
%     handles.Cfg.T1AtlasMaskFile =theMaskfile;
%     guidata(hObject, handles);
% else
%     errordlg(sprintf('The mask file "%s" does not exist!\n Please re-check it.', theMaskfile));
% end
% guidata(hObject, handles);

% --- Executes on button press in btnSelectAtlasFile.
% function btnSelectAtlasFile_Callback(hObject, eventdata, handles)
% [filename, pathname] = uigetfile({'*.img;*.nii;*.nii.gz','Brain Image Files (*.img;*.nii;*.nii.gz)';'*.*', 'All Files (*.*)';}, ...
%     'Pick a mask');
% if ~(filename==0)
%     handles.Cfg.T1AtlasMaskFile =[pathname filename];
% elseif ~( exist(handles.Cfg.T1AtlasMaskFile, 'file')==2)    
%     set(handles.btnSelectAtlasFile, 'Enable','off');
% end
% guidata(hObject, handles);
% UpdateDisplay(handles);



% --- Executes on button press in ckboxR2SN.
function ckboxR2SN_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsR2SN = 1;
else
    handles.Cfg.IsR2SN = 0;
end
setT1OnekeyPipelineDisplayStatus(hObject,handles,'off');
guidata(hObject, handles);
UpdateDisplay(handles);



% --- Executes on button press in ckboxReportWholeBrain.
function ckboxReportWholeBrain_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsReportWholeBrain = 1;
else
    handles.Cfg.IsReportWholeBrain = 0;
end
guidata(hObject, handles);
UpdateDisplay(handles);


% --- Executes on button press in ckboxReportInIndividualSpace.
function ckboxReportInIndividualSpace_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsReportInIndividualSpace = 1;
else
    handles.Cfg.IsReportInIndividualSpace = 0;
end
guidata(hObject, handles);
UpdateDisplay(handles);


function ckboxMinimumOutput_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsMinimumOutput = 1;
else
    handles.Cfg.IsMinimumOutput = 0;
end
guidata(hObject, handles);
UpdateDisplay(handles);

% --- Executes during object creation, after setting all properties.
function edtParallelWorkers_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtParallelWorkers_Callback(hObject, eventdata, handles)
handles.Cfg.ParallelWorkers=str2double(get(hObject,'String'));
guidata(hObject, handles);
UpdateDisplay(handles);


% --- Executes during object creation, after setting all properties.
function edtStartingDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtStartingDir_Callback(hObject, eventdata, handles)
uiwait(msgbox({'If you do not start with raw DICOM images, you need to specify the Starting Directory Name.';...
        'E.g. "T1ImgCA" means you start with images which have been coregistration and feature extracted based on default or user define atlas.';...
        '';...
        'Abbreviations:';...
        'C - Coregister';...
        'S - Segment';...
        'A - Atlas based feature extracted';...
        '';...
        },'Tips for Starting Directory Name'));

handles.Cfg.StartingDir=get(hObject,'String');
handles=CheckCfgParametersBeforeRun(handles);
guidata(hObject, handles);
UpdateDisplay(handles);


% --- Executes on button press in btnLoadConfig.
function handles=btnLoadConfig_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile({'*.mat'}, 'Load Parameters From');
if ischar(filename)
    load([pathname,filename]);
    SetLoadedData(hObject,handles, Cfg);
end

function SetLoadedData(hObject,handles, Cfg)	
handles.Cfg=Cfg;
guidata(hObject, handles);
UpdateDisplay(handles);
    
% --- Executes on button press in btnSaveConfig.
function btnSaveConfig_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile({'*.mat'}, 'Save Parameters As');
if ischar(filename)
    Cfg=handles.Cfg;
    save(['',pathname,filename,''], 'Cfg');
end

% --- Executes on button press in btnQuitiBrain.
function btnQuitiBrain_Callback(hObject, eventdata, handles)
close(handles.figiBrainMain);

% --- Executes on button press in btnRuniBrain.
function btnRuniBrain_Callback(hObject, eventdata, handles)
[handles, CheckingPass]=CheckCfgParametersBeforeRun(handles);
if CheckingPass==0
    return
end
RawBackgroundColor=get(handles.btnRuniBrain ,'BackgroundColor');
RawForegroundColor=get(handles.btnRuniBrain ,'ForegroundColor');
iBrainPath = fileparts(which('iBrain.m'));
set(handles.btnRuniBrain ,'Enable', 'off','BackgroundColor', 'red','ForegroundColor','green');
handles.Cfg.CoregisterTemplate=[iBrainPath,filesep,'Template',filesep,'MNI152_T1_1mm.nii.gz'];
handles.Cfg.T1AtlasMaskFile=[iBrainPath,filesep,'Atlas',filesep,'BN_Atlas_246_1mm.nii'];
Cfg=handles.Cfg; 
Datetime=fix(clock); 
save([handles.Cfg.DataProcessDir,filesep,'iBrain_AutoSave_',num2str(Datetime(1)),'_',num2str(Datetime(2)),'_',num2str(Datetime(3)),'_',num2str(Datetime(4)),'_',num2str(Datetime(5)),'.mat'], 'Cfg'); %Added by YAN Chao-Gan, 100130.
if handles.Cfg.IsMinimumOutput==0
    [Error]=iBrain_run(handles.Cfg);
else
    [Error]=iBrainMiniOutput_run(handles.Cfg);
end

if ~isempty(Error)
    uiwait(msgbox(Error,'Errors were encountered while processing','error'));
end
set(handles.btnRuniBrain ,'Enable', 'on','BackgroundColor', RawBackgroundColor,'ForegroundColor',RawForegroundColor);
UpdateDisplay(handles);


%% Check if the configuration parameters is correct
function [handles, CheckingPass]=CheckCfgParametersBeforeRun(handles)    
CheckingPass=0;
if isempty (handles.Cfg.WorkingDir)
    uiwait(msgbox('Please set the working directory!','Configuration parameters checking','warn'));
    return
end

if 7==exist([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir],'dir') && ~strcmp(handles.Cfg.StartingDir,'')
    if isempty (handles.Cfg.SubjectID)
        Dir=dir([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir]);
        if strcmpi(Dir(3).name,'.DS_Store')  
            StartIndex=4;
        else
            StartIndex=3;
        end
        for i=StartIndex:length(Dir)
            handles.Cfg.SubjectID=[handles.Cfg.SubjectID;{Dir(i).name}];
        end
    end
    
    if ~(strcmpi(handles.Cfg.StartingDir,'T1Raw') || strcmpi(handles.Cfg.StartingDir,'T1Img')) 
        DirImg=dir([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir,filesep,handles.Cfg.SubjectID{1},filesep,'*.img']);
        if isempty(DirImg)  
            DirImg=dir([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir,filesep,handles.Cfg.SubjectID{1},filesep,'*.nii']);
            if length(DirImg)==0
                DirImg=dir([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir,filesep,handles.Cfg.SubjectID{1},filesep,'*.nii.gz']);% Search .nii.gz and unzip; YAN Chao-Gan, 120806.
                if length(DirImg)==1 
                    gunzip([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir,filesep,handles.Cfg.SubjectID{1},filesep,DirImg(1).name]);
                    Nii  = nifti([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDirName,filesep,handles.Cfg.SubjectID{1},filesep,DirImg(1).name(1:end-3)]);
                    delete([handles.Cfg.WorkingDir,filesep,handles.Cfg.StartingDir,filesep,handles.Cfg.SubjectID{1},filesep,DirImg(1).name(1:end-3)]);
                end
            end
        end
    end
else
    uiwait(msgbox(['Please arrange each subject''s NIFTI images in one directory, and then put them in your defined Starting Directory Name "',handles.Cfg.StartingDir,'" directory under the working directory!'],'Configuration parameters checking','warn'));
    return
end

CheckingPass=1;
UpdateDisplay(handles);

%% Update All the uiControls' display on the GUI
function UpdateDisplay(handles)
set(handles.edtWorkingDir ,'String', handles.Cfg.WorkingDir);

if size(handles.Cfg.SubjectID,1)>0
    theOldIndex =get(handles.listSubjectID, 'Value');
    set(handles.listSubjectID, 'String',  handles.Cfg.SubjectID , 'Value', 1);
    theCount =size(handles.Cfg.SubjectID,1);
    if (theOldIndex>0) && (theOldIndex<= theCount)
        set(handles.listSubjectID, 'Value', theOldIndex);
    end
else
    set(handles.listSubjectID, 'String', '' , 'Value', 0);
end
set(handles.ckboxonekeyFromT1Raw, 'Value', handles.Cfg.IsOneKeyRunFromT1Raw);
set(handles.ckboxonekeyFromT1Img, 'Value', handles.Cfg.IsOneKeyRunFromT1Img);
set(handles.ckboxT1DICOM2NIFTI, 'Value', handles.Cfg.IsNeedConvertT1DCM2NII);
set(handles.ckboxSegmentation, 'Value', handles.Cfg.IsSegment);
set(handles.ckboxCoregister, 'Value', handles.Cfg.IsCoregister);
set(handles.ckboxAtlasFeatureExtraction, 'Value', handles.Cfg.IsExtractT1AtlasFeature);
% set(handles.edtAtlasFile ,'String', handles.Cfg.T1AtlasMaskFile);
set(handles.ckboxR2SN, 'Value', handles.Cfg.IsR2SN);


set(handles.ckboxReportWholeBrain, 'Value', handles.Cfg.IsReportWholeBrain);
set(handles.ckboxReportInIndividualSpace, 'Value', handles.Cfg.IsReportInIndividualSpace);
set(handles.ckboxMinimumOutput, 'Value', handles.Cfg.IsMinimumOutput);
set(handles.edtParallelWorkers ,'String', num2str(handles.Cfg.ParallelWorkers));
set(handles.edtStartingDir ,'String', handles.Cfg.StartingDir);

drawnow;
