%install all required python and matlab toolbox
iBrainPath=fileparts(which('iBrain.m'));
addpath(genpath(iBrainPath));
if isunix
    system('python -m pip install antspy');
else
    %obtain python version, and install related whl file based on python,
    %only support python version>=3.7
    [~,cmdout]=system('python --version');
    pattern='\d*\.\d*';
    pythonversion=regexp(cmdout,pattern,'match');
    pythonversion=pythonversion{1};
    pythonversion(find(pythonversion=='.'))=[];
    system(['python -m pip install antspyx-0.3.8-cp',pythonversion,'-win_amd64.whl'])%can replace to download from website
end
system('python -m pip install -r requirements.txt');
