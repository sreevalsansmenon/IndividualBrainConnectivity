% This program can be used to run the subprograms and output the results
% Menon, S. S., & Krishnamurthy, K. (2019). A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Biological Sex Using Intrinsic Individual Brain Connectivity. Scientific Reports, 9(1), 5729. https://doi.org/10.1038/s41598-019-42090-4

clc                         % Main Program
clear                       % Clear workspace
addpath('./Data')           % Add path to data
mkdir('Data_Save')          % Create folder to save data
addpath('./Data_Save')      % Add path to save data
resting_state_fingerprint   % Outputs results explained in Subject identification using 3T and 7T rfMRI data
edge_analysis               % Outputs results explained in Edge analysis in subject identification
rest_task_fingerprint       % Outputs results explained in Subject identification using 3T tfMRI data
sex_identification          % Outputs results explained in Sex identification using 3T rfMRI data
sex_differential            % Outputs results explained in Differential power of edges in Sex identification
