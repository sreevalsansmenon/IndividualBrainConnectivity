# IndividualBrainConnectivity
The codes shared in this folder can be used to replicate the results explained in “Menon, S. S. &amp; Krishnamurthy, K. A comparison of static and dynamic functional connectivities for identifying subjects and biological sex using intrinsic individual brain connectivity. Sci. Reports 9:5729, https://doi.org/10.1038/s41598-019-42090-4 (2019)”.

Initially preprocessed data has to be downloaded from the HCP website. For 3T, resting state file named “rfMRI_REST_1/2_LR/RL_hp2000_clean.nii.gz” has to be downloaded for sessions 1_LR, 1_RL, 2_LR and 2_RL and task file named
“tfMRI_task_LR/RL.nii.gz” for tasks EMOTION, GAMBLING, LANGUAGE, MOTOR, RELATIONAL, SOCIAL and WM LR and RL configurations. For 7T, resting state file named “rfMRI_REST1/2_7T_PA/AP_hp2000_clean.nii.gz” has to be downloaded for
sessions 1 PA, 1 AP, 2 PA and 2 AP.

The next step is to download the Stanford functional atlas which can be downloaded from the link https://findlab.stanford.edu/functional_ROIs.html. To extract time series for the parcelled ROIs, we can use the following FSl command:

fslmeants -i &lt;input_path&gt; -m &lt;mask_path&gt; -o &lt;output_path&gt;

Input_path can be pointed to the downloaded fMRI data, mask_path to the ROI mask downloaded and output_path to the location where the time series text file has to be saved.

Finally, you have to save the time series as a MATLAB file of size 90 x time points (use MATLAB, for our program 343RfmriTimeseries consists of 1372 variables named subject_session.mat (90*1200) and a subject variable with 1372 cells to identify session and subject, similarly for other resting state and task data). Now we can use the MATLAB files provided here to replicate the results. The main program can be called to output the results presented.

** MAKE USE OF PARFOR FOR COMPUTATIONAL EFFICIENCY **

** Codes shared here are in a redundant coding fashion for better clarity **
