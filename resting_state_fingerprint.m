% Subject identification using 3T and 7T rfMRI data
% Menon, S. S., & Krishnamurthy, K. (2019). A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Biological Sex Using Intrinsic Individual Brain Connectivity. Scientific Reports, 9(1), 5729. https://doi.org/10.1038/s41598-019-42090-4
% This program finds and save the static and dynamic fc and performs connectome fingerprinting (Results section 1).  Note the data input to the program is 4 resting state data available for subjects
% Here we have subject variable as cell with names LR1/LR2/RL1/RL2_subject row wise
% thus subject variable is 1372*1 cell (4*343) and timeseries was saved as
% LR1/LR2/RL1/RL2_subject.mat(1200*90) (1372 variables)
% For 7T replace with 7T timeseries data 
%% Pearson sFC 3T
% This section of program loads timeseries and calculates Pearson sFC
clear                                                       % Clear workspace
load 343RfmriTimeseries.mat                                 % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                       % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                    % Temporal resolution
fnq=1/(2*TR);                                               % Nyquist frequency
flp = 1/(60*1);                                             % lowpass frequency of filter (Hz)
fhi = 0.15;                                                 % highpass
Wn=[flp/fnq fhi/fnq];                                       % butterworth bandpass non-dimensional frequency
k=2;                                                        % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                 % construct the filter
clear fnq flp fhi Wn k

for s=1:size(Subjects,1)                                    % Loop for finding correlation     
    Subject=char(Subjects(s));                              % Loading subject
    eval(['signaldata = ' Subject ';']);                    % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                         % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);            % Perform filtering
    pearsonfc= Correlation(signaldata,1,'corr',1200);       % Finding pearson correlation
    subfc=pearsonfc(uppertriangle);                         % Computing upper triangle elements
    eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);            % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject % Clear unused variables

s=textread('HCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                        % variable initialization
for i=1:size(Subjects,1)/4                                  % Looping over number of subjects
            a=[4*i-3:4*i];                                  % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                       % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));               % Reads session name                                         
                eval(['ints = ' Subject ';']);              % Loads session FC 
                su(:,j)=double(ints);                       % Saving the FC to a variable su
                eval(['clear ' Subject])                    % Clear subject variable
            end                                             % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);  % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);  % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                          % Increment k variable
end                                                         % End of loop
Subjects=[];                                                % Clear subject variable
Subjects=Subj';                                             % Saving subject ids 
clear s su Subj Subject a i ints j k                        % Clear unused variable
save('./Data_Save/Pearson_Rest_StaticFC.mat')               % Save pearson FC

%% Partial sFC 3T
% This section of program loads timeseries and calculates Partial sFC
clear                                                       % Clear workspace
load 343RfmriTimeseries.mat                                 % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                       % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                    % Temporal resolution
fnq=1/(2*TR);                                               % Nyquist frequency
flp = 1/(60*1);                                             % lowpass frequency of filter (Hz)
fhi = 0.15;                                                 % highpass
Wn=[flp/fnq fhi/fnq];                                       % butterworth bandpass non-dimensional frequency
k=2;                                                        % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                 % construct the filter
clear fnq flp fhi Wn k                                      % Clear variables

for s=1:size(Subjects,1)                                    % Loop for saving the primary data for comparison     
    Subject=char(Subjects(s));                              % Loading subject
    eval(['signaldata = ' Subject ';']);                    % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                         % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);            % Perform filtering
    partialfc= Correlation(signaldata,1,'icov',1200);       % Finding partial correlation
    subfc=partialfc(uppertriangle);                         % Computing upper triangle elements
    eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);            % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn partialfc s Subject % Clear unused variables

s=textread('HCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                        % variable initialization
for i=1:size(Subjects,1)/4                                  % Looping over number of subjects
            a=[4*i-3:4*i];                                  % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                       % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));               % Reads session name                                         
                eval(['ints = ' Subject ';']);              % Loads session FC 
                su(:,j)=double(ints);                       % Saving the FC to a variable su
                eval(['clear ' Subject])                    % Clear subject variable
            end                                             % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);  % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);  % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                          % Increment k variable
end                                                         % End of loop
Subjects=[];                                                % Clear subject variable
Subjects=Subj';                                             % Saving subject ids 
clear s su Subj Subject a i ints j k                        % Clear unused variable
save('./Data_Save/Partial_Rest_StaticFC.mat')               % Save partial FC
clear                                                       % Clear workspace

%% Dynamic FC 3T
% This section of program loads timeseries and calculates Pearson dFC
load 343RfmriTimeseries.mat                                     % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
windowsize=85;                                                  % Window size
nstates=4;                                                      % Number of states
% Bandpass filter settings
TR=0.72;                                                        % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(windowsize*TR);                                        % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variable
for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    t=1;i=1;                                                    % Initialiaztion of Variables
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    while t+windowsize<=1201                                    % Loop to perform window analysis
        dfc=Correlation(zscore(signaldata(t:t+windowsize-1,:)')',1,'corr',windowsize);  % Finding Dynamic FC matrix over window length w
        wdfc(s,i,:)=dfc(uppertriangle);                         % Saving the upper triangular elements of DFC
        i=i+1;t=t+5;                                            % Incrementing Variables
    end                                                 
    X=squeeze(wdfc(s,:,:));                                     % Combing the DFCs for the subject
    [IDX, C, SUMD, D]=kmeans(X,nstates,'Distance','cosine','Replicates',20,'Display','off'); % Performin kmeans to find repeating nstates 
    
    for k=1:nstates                                             % Loop to find the max repeated states
        a(k)=sum(IDX==k);                                       % Saving number of times each states being repeated
    end                                                     
    
    [b,a]=sort(a);                                              % Finding the order of states being repeated maximum in ascending order                                 
    
    Co=C(a(1),:);                                               % Saving the first highest repeated state
    for k=2:nstates                                             % Loop for saving the repeated states for the subject
        Co=horzcat(Co,C(a(k),:));                               % Saving the centroid values for states in ascending order
    end                                                	
    subfc=Co';                                                  % Saving the dynamic FC of the subject 
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject C Co D IDX task wdfc X % Clear variable
s=textread('HCP_Subjects_Rest.txt');                            % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Pearson_Rest_DynamicFC.mat')                  % Save partial FC
clear                                                           % Clear workspace
%% Pearson sFC 7T
% This section of program loads timeseries and calculates Pearson sFC for
% 7T data
clear                                                           % Clear workspace
load 7TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=1;                                                           % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(60*1);                                                 % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variable

for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    pearsonfc= Correlation(signaldata,1,'corr',900);            % Finding pearson correlation
    subfc=pearsonfc(uppertriangle);                             % Computing upper triangle elements
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);                % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject % Clear unused variables

s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Pearson_Rest7T_StaticFC.mat')                 % Save pearson FC

%% Partial sFC 3T
% This section of program loads timeseries and calculates Partial sFC for
% 7T data
clear                                                           % Clear workspace
load 7TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=1;                                                           % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(60*1);                                                 % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variable

for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    partialfc= Correlation(signaldata,1,'icov',900);            % Finding partial correlation
    subfc=partialfc(uppertriangle);                             % Computing upper triangle elements
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);                % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn partialfc s Subject % Clear unused variables

s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Partial_Rest7T_StaticFC.mat')                 % Save partial FC
clear                                                           % Clear workspace

%% Dynamic FC 7T
% This section of program loads timeseries and calculates Pearson dFC for
% 7T data
load 7TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
windowsize=85;                                                  % Window size
nstates=4;                                                      % Number of states
% Bandpass filter settings
TR=1;                                                           % Temporal resolution                                                       
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(windowsize*TR);                                        % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k
for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    t=1;i=1;                                                    % Initialiaztion of Variables
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    while t+windowsize<=900                                     % Loop to perform window analysis
        dfc=Correlation(zscore(signaldata(t:t+windowsize-1,:)')',1,'corr',windowsize);  % Finding Dynamic FC matrix over window length w
        wdfc(s,i,:)=dfc(uppertriangle);                         % Saving the upper triangular elements of DFC
        i=i+1;t=t+5;                                            % Incrementing Variables
    end                                                 
    X=squeeze(wdfc(s,:,:));                                     % Combing the DFCs for the subject
    [IDX, C, SUMD, D]=kmeans(X,nstates,'Distance','cosine','Replicates',20,'Display','off'); % Performin kmeans to find repeating nstates 
    
    for k=1:nstates                                             % Loop to find the max repeated states
        a(k)=sum(IDX==k);                                       % Saving number of times each states being repeated
    end                                                     
    
    [b,a]=sort(a);                                              % Finding the order of states being repeated maximum in ascending order                                 
    
    Co=C(a(1),:);                                               % Saving the first highest repeated state
    for k=2:nstates                                             % Loop for saving the repeated states for the subject
        Co=horzcat(Co,C(a(k),:));                               % Saving the centroid values for states in ascending order
    end                                                	
    subfc=Co';                                                  % Saving the dynamic FC of the subject 
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject C Co D IDX task wdfc X % Clear variable
s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Pearson_Rest7T_DynamicFC.mat')                % Save partial FC
clear
%% Pearson sFC 37T
% This section of program loads timeseries and calculates Pearson sFC for
% 164 3T subjects
clear                                                           % Clear workspace
load 3TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                        % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(60*1);                                                 % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variable

for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    pearsonfc= Correlation(signaldata,1,'corr',1200);           % Finding pearson correlation
    subfc=pearsonfc(uppertriangle);                             % Computing upper triangle elements
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);                % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject % Clear unused variables

s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Pearson_Rest3T_StaticFC.mat')                 % Save pearson FC

%% Partial sFC 3T
% This section of program loads timeseries and calculates Partial sFC for
% 164 3T subjects
clear                                                           % Clear workspace
load 3TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                        % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(60*1);                                                 % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variables

for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata       
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    partialfc= Correlation(signaldata,1,'icov',1200);           % Finding partial correlation
    subfc=partialfc(uppertriangle);                             % Computing upper triangle elements
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);                % Saving the result
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn partialfc s Subject % Clear unused variables

s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Partial_Rest3T_StaticFC.mat')                 % Save partial FC
clear                                                           % Clear workspace

%% Dynamic FC 3T
load 3TRfmriTimeSeries_164.mat                                  % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                           % Finding the Upper Diagonal Elements Index 
windowsize=85;                                                  % Window size
nstates=4;                                                      % Number of states
% Bandpass filter settings
TR=0.72;                                                        % Temporal resolution
fnq=1/(2*TR);                                                   % Nyquist frequency
flp = 1/(windowsize*TR);                                        % lowpass frequency of filter (Hz)
fhi = 0.15;                                                     % highpass
Wn=[flp/fnq fhi/fnq];                                           % butterworth bandpass non-dimensional frequency
k=2;                                                            % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                     % construct the filter
clear fnq flp fhi Wn k                                          % Clear variables

for s=1:size(Subjects,1)                                        % Loop for saving the primary data for comparison     
    t=1;i=1;                                                    % Initialiaztion of Variables
    Subject=char(Subjects(s));                                  % Loading subject
    eval(['signaldata = ' Subject ';']);                        % Extract timeseries data to signaldata
    signaldata=detrend(signaldata);                             % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);                % Perform filtering
    while t+windowsize<=1201                                    % Loop to perform window analysis
        dfc=Correlation(zscore(signaldata(t:t+windowsize-1,:)')',1,'corr',windowsize);  % Finding Dynamic FC matrix over window length w
        wdfc(s,i,:)=dfc(uppertriangle);                         % Saving the upper triangular elements of DFC
        i=i+1;t=t+5;                                            % Incrementing Variables
    end                                                 
    X=squeeze(wdfc(s,:,:));                                     % Combing the DFCs for the subject
    [IDX, C, SUMD, D]=kmeans(X,nstates,'Distance','cosine','Replicates',20,'Display','off'); % Performin kmeans to find repeating nstates 
    
    for k=1:nstates                                             % Loop to find the max repeated states
        a(k)=sum(IDX==k);                                       % Saving number of times each states being repeated
    end                                                     
    
    [b,a]=sort(a);                                              % Finding the order of states being repeated maximum in ascending order                                 
    
    Co=C(a(1),:);                                               % Saving the first highest repeated state
    for k=2:nstates                                             % Loop for saving the repeated states for the subject
        Co=horzcat(Co,C(a(k),:));                               % Saving the centroid values for states in ascending order
    end                                                	
    subfc=Co';                                                  % Saving the dynamic FC of the subject 
    eval(['clear ' Subject ';']);                               % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject C Co D IDX task wdfc X
s=textread('7T3THCP_Subjects_Rest.txt');                        % Loads text file with subject ids
k=0;                                                            % variable initialization
for i=1:size(Subjects,1)/4                                      % Looping over number of subjects
            a=[4*i-3:4*i];                                      % Finding the session row for a subject(Note it depends how you save timeseries)
            for j=1:4                                           % Looping over four sessions of a subject
                Subject=char(Subjects(a(j)));                   % Reads session name                                         
                eval(['ints = ' Subject ';']);                  % Loads session FC 
                su(:,j)=double(ints);                           % Saving the FC to a variable su
                eval(['clear ' Subject])                        % Clear subject variable
            end                                                 % Loop terminates
            Subj(k+1)=cellstr([ 'Rest1_' num2str(s(i)) ]);      % Variable name for day 1 FC
            Subj(k+2)=cellstr([ 'Rest2_' num2str(s(i)) ]);      % Variable name for day 2 FC
            eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,3))./2;' ]);% Finding day 1 FC
            eval([cell2mat(Subj(k+2)) '= (su(:,2)+su(:,4))./2;' ]);% Finding day 2 FC
            k=k+2;                                              % Increment k variable
end                                                             % End of loop
Subjects=[];                                                    % Clear subject variable
Subjects=Subj';                                                 % Saving subject ids 
clear s su Subj Subject a i ints j k                            % Clear unused variable
save('./Data_Save/Pearson_Rest3T_DynamicFC.mat')                % Save partial FC
clear
%% Subject identification using pearson sFC 3T3T
% This section predicts subjects using pearson sFC for 3T data
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest_StaticFC.mat                              % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects 
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    % Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    % Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    % Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  3T-3T Pearson sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction


%% Subject identification using partial sFC 3T3T
% This section predicts subjects using partial sFC for 3T data
clear                                                           % Clear variables
for pr=1:1000                                                   % Perform 1000 runs
    load Partial_Rest_StaticFC.mat                              % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %% Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    %% Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  3T-3T Partial sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
clear                                                           % Clear workspace

%% Subject identification using dFC 3T3T
% This section predicts subjects using pearson dFC for 3T data
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest_DynamicFC.mat                             % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects 
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the dynamic connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %% Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    %% Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for 3T-3T Pearson dFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
clear                                                           % Clear workspace
%% Subject identification using pearson sFC 7T7T
% This section predicts subjects using pearson sFC for 7T data
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest7T_StaticFC.mat                            % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2 
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    % Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    % Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    % Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  7T-7T Pearson sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction

%% Subject identification using partial sFC 7T7T
 % This section predicts subjects using partial sFC for 7T data
 clear                                                          % Clear workspace
for pr=1:1000                                                   % Perform 1000 runs
    load Partial_Rest7T_StaticFC.mat                            % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %% Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    %% Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  7T-7T Partial sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
clear                                                           % Clear workspace

%% Subject identification using dFC 7T7T
% This section predicts subjects using pearson dFC for 7T data
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest7T_DynamicFC.mat                           % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2 
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the dynamic connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %% Prediction for unknown data

    for s=1:size(refpredsessions,2)                             % Loop for predicting data for comparison  
        Subject=char(Subjects(refpredsessions(s)));             % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    %% Comparing the predicted results with actual subject 

    actual=ceil([1:1*ntrials]./1);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for 7T-7T Pearson dFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
             
%% Subject identification using pearson sFC 3T7T
% This section predicts subjects using pearson sFC for 7T-3T data
clear                                                           % Clear workspace
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest7T_StaticFC.mat                            % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2 
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    % Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    % Prediction for unknown data
        load Pearson_Rest3T_StaticFC.mat                        % Loads timeseries data
    for s=1:size(sessiondata,2)                                 % Loop for predicting data for comparison  
        Subject=char(Subjects(sessiondata(s)));                 % Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    % Comparing the predicted results with actual subject 

    actual=ceil([1:2*ntrials]./2);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  7T-3T Pearson sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction

%% Subject identification using partial sFC 3T7T
 % This section predicts subjects using partial sFC for 7T-3T data
 clear                                                          % Clear workspace
for pr=1:1000                                                   % Perform 1000 runs
    load Partial_Rest7T_StaticFC.mat                            % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    % Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %Prediction for unknown data
    load Partial_Rest3T_StaticFC.mat                            % Loads timeseries data

    for s=1:size(sessiondata,2)                             	% Loop for predicting data for comparison  
        Subject=char(Subjects(sessiondata(s)));             	% Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    % Comparing the predicted results with actual subject 

    actual=ceil([1:2*ntrials]./2);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for  7T-3T Partial sFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
clear                                                           % Clear workspace

%% Subject identification using dFC 3T7T
% This section predicts subjects using pearson dFC for 7T-3T data
for pr=1:1000                                                   % Perform 1000 runs
    load Pearson_Rest7T_DynamicFC.mat                           % Loads timeseries data
    ntrials=size(Subjects,1)/2;                                 % ntrials is the number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];% Finding random subjects from size(Subjects,1)/2 
    for i=1:ntrials                                             % Loop  to find correspoding column in subject data
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
    end                                                         

    basesessions=[0:2:2*ntrials-2]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject
    predsessions=[1:2*ntrials];predsessions(basesessions)=[];   % Finding session that are not in the baseline
    refpredsessions=sessiondata(predsessions);                  % Getting reference column of session to be predicted

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['fcdata = ' Subject ';']);                        % Extract timeseries data to signaldata       
        comparefc(:,s)=fcdata;                                  % Saving the dynamic connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                         

    %% Prediction for unknown data
    load Pearson_Rest3T_DynamicFC.mat                           % Loads timeseries data

    for s=1:size(sessiondata,2)                             	% Loop for predicting data for comparison  
        Subject=char(Subjects(sessiondata(s)));             	% Loading subject                 
        eval(['fcdata = ' Subject ';']);                        % Extract subjectFC   
        pcorr=1-pdist2(comparefc',fcdata','cosine');            % Finding cosine similarity
        [corvalue(s),pred(s)]=max(pcorr);                       % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end  

    %% Comparing the predicted results with actual subject 

    actual=ceil([1:2*ntrials]./2);                              % Actual Subject Data                                        
    accu=sum(actual==pred);                                     % Comparing prediction with actual data
    r(pr)=100*accu/s;                                           % Percentage accuracy
    clear i s signaldata                                        % Clearing the data which is no longer needed                       

end

fprintf('The mean accuracy for 7T-3T Pearson dFC is  %d percent\n', mean(r)); % Output the accuracy in prediction
clear                                                           % Clear workspace                            
