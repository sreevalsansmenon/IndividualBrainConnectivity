# IndividualBrainConnectivity
The codes shared in this folder can be used to replicate the results explained in “A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Sex Using Intrinsic Individual Brain Connectivity”. 

Initially preprocessed data has to be downloaded from HCP website. For 3T, resting state file named “rfMRI_REST_1/2_LR/RL_hp2000_clean.nii.gz”  has to be downloaded for sessions 1_LR, 1_RL, 2_LR and 2_RL and task file named “tfMRI_task_LR/RL.nii.gz” for tasks EMOTION, GAMBLING, LANGUAGE, MOTOR, RELATIONAL, SOCIAL and WM LR and RL configurations . For 7T, resting state file named “rfMRI_REST1/2_7T_PA/AP_hp2000_clean.nii.gz”  has to be downloaded for sessions 1 PA, 1 AP, 2 PA and 2 AP. 

The next step is to download the stanford functional atlas which can be downloaded from the link https://findlab.stanford.edu/functional_ROIs.html. To extract time series for the parcelled ROIs we can use following FSl command:
fslmeants -i <input_path> -m <mask_path> -o <output_path> 
Input_path can be pointed to downloaded fMRI data,, mask_path to the ROI mask downloaded and output_path to the location where timeseries text file has to saved.

Finally you have to save the timeseries as a matlab file of size 90 x timepoints (use matlab, for our program 343RfmriTimeseries consists of 1372 variables named subject_session.mat(90*1200) and a subject variable with 1372 cell to identify session and subject, similarly for other resting state and task data). Now we can use the matlab files provided here to replicate the results. The main program can be called to output the results presented.

** MAKE USE OF PARFOR FOR COMPUTATIONAL EFFICIANCY **
