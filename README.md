# iBrain introduction
iBrain toolbox is designed for multi-scale structural feature extraction and report generation for neurodegenerative diseases. The toolbox provide simplified GUI and batch script which support subject level brain report from original DICOM or NIFTI T1 image. iBrain is written by MATLAB and Python scripts which can be successfully on Window and Linux operation system. 
<div align=center>
<img src=https://github.com/YongLiuLab/iBrain/assets/20011474/465e0401-19ce-4b5b-8e04-b1a90de9bc2b>
</div>
<div align=center>
<img src=https://github.com/YongLiuLab/iBrain/assets/20011474/e330cebb-50d2-42d4-b246-086a5999c407>
</div>
# Pre-requests
• SPM12 and CAT12: SPM is freely available to the (neuro) imaging community andrepresents the implementation of the theoretical concepts of Statistical Parametric Mapping in a complete analysis package, our toolbox use SPM for basic data I/O options. CAT12 is a comprehensive MATLAB toolbox specifically designed for the analysis and processing of structural brain MRI data, our toolbox use CAT12 to generate basic structural features. It provides a range of advanced tools and functions to facilitate voxel-based morphometry (VBM) and other morphometric analyses. SPM12 toolbox needs to be copy into dependenices folder, and CAT12 needs to be put in dependenices/spm12/toolbox folder. 
• Brant: Brant is a MATLAB toolbox designed for the analysis and processing of neuroimaging data. It offers a comprehensive set of tools and functions to facilitate various types of analyses and investigations related to brain connectivity and network properties, our toolbox use it for additional data I/O options. Brant needs to be but into dependenices folder. 

Note: Apart from the pre-request toolboxes, you need to install Python 3.8 or 3.9 on your computer and run set_up.m to install required python packages to run iBrain successfully on your own computer. To install antspy in windows, you need to download compete version of .whl file from following website:  https://github.com/ANTsX/ANTsPy/actions/runs/4524955037 or antspy for python version 3.9 from link https://www.dropbox.com/scl/fi/swjvcz1z3xwk3gahq9a7r/antspyx-0.3.8-cp39-cp39-win_amd64.whl?rlkey=ij4ce471x32lpss85sylsv28j&dl=0. Another large file train_data.mat is required and needs to be downloaded and put into folder model_data from link: https://www.dropbox.com/scl/fi/kev2l332xtc31g22qrwnp/train_data.mat?rlkey=l33m7kmu6y81w13o17zofohjd&dl=0

You can also download full version of iBrain from link: https://www.dropbox.com/scl/fo/7xsh5ewq5jz94lbgvx7dj/h?rlkey=0tyzrj9arg921j4fh6a93y0kg&dl=0
