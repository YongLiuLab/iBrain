#!/bin/bash
nohup matlab -nodisplay -nosplash -nojvm -nodesktop -r "addpath(genpath('/data/yhzhang/matlab_toolbox/NIfTI_20140122'));batch_extract_subject_GMWM_from_server;" 1>process_log.txt 2>process_err.txt &
