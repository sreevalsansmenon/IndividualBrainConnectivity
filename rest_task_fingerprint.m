%  Subject identification using 3T tfMRI data 
% A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Sex Using Intrinsic Individual Brain Connectivity
% Sreevalsan S Menon, K Krishnamurthy
% This program  performs Subject identification using 3T tfMRI data
% (Results section 3). save task.mat files as in resting state(642/task * 7files)

%% Rest Pearson sFC
%  This part of program finds pearson correlation for 321 subjects 3T rfMRI  
clear                                                       % Clear workspace
load 321RfmriTimeseries.mat                                 % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                       % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                    % Temporal resolution
fnq=1/(2*TR);                                               % Nyquist frequency
flp = 1/(60*1);                                             % lowpass frequency of filter (Hz)
fhi = 0.15;                                                 % highpass
Wn=[flp/fnq fhi/fnq];                                       % butterworth bandpass non-dimensional frequency
k=2;                                                        % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                 % construct the filter
clear fnq flp fhi Wn k                                      % Clear workspace

for s=1:size(Subjects,1)                                    % Loop for saving the primary data for comparison     
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

s=textread('HCP_Subjects_Task.txt');                        % Loads text file with subject ids
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
save('./Data_Save/Pearson_Rest_Task_StaticFC.mat')          % Save pearson FC

%% Rest Partial sFC
%  This part of program finds partial correlation for 321 subjects 3T rfMRI  
clear                                                       % Clear workspace
load 321RfmriTimeseries.mat                                 % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                       % Finding the Upper Diagonal Elements Index 
% Bandpass filter settings
TR=0.72;                                                    % Temporal resolution
fnq=1/(2*TR);                                               % Nyquist frequency
flp = 1/(60*1);                                             % lowpass frequency of filter (Hz)
fhi = 0.15;                                                 % highpass
Wn=[flp/fnq fhi/fnq];                                       % butterworth bandpass non-dimensional frequency
k=2;                                                        % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                 % construct the filter
clear fnq flp fhi Wn k                                      % Clear workspace

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

s=textread('HCP_Subjects_Task.txt');                        % Loads text file with subject ids
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
save('./Data_Save/Partial_Rest_Task_StaticFC.mat')                           % Save partial FC
clear                                                       % Clear workspace

%% Rest Dynamic FC
%  This part of program finds dynamic pearson correlation for 321 subjects 3T rfMRI  
load 321RfmriTimeseries.mat                                 % Loads timeseries data; the Subjects cell shows the subject name and session name
uppertriangle=find(triu(ones(90),1));                       % Finding the Upper Diagonal Elements Index 
windowsize=85;
nstates=4;
% Bandpass filter settings
TR=0.72;                                                    % Temporal resolution
fnq=1/(2*TR);                                               % Nyquist frequency
flp = 1/(windowsize*TR);                                    % lowpass frequency of filter (Hz)
fhi = 0.15;                                                 % highpass
Wn=[flp/fnq fhi/fnq];                                       % butterworth bandpass non-dimensional frequency
k=2;                                                        % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);                                 % construct the filter
clear fnq flp fhi Wn k                                      % Clear workspace
for s=1:size(Subjects,1)                                    % Loop for saving the primary data for comparison     
    t=1;i=1;                                                % Initialiaztion of Variables
    Subject=char(Subjects(s));                              % Loading subject
    eval(['signaldata = ' Subject ';']);                    % Extract timeseries data to signaldata
    signaldata=detrend(signaldata);                         % Normalize the data
    signaldata=filtfilt(bfilt,afilt,signaldata);            % Filter the signaldata
    while t+windowsize<=1201                                % Loop to perform window analysis
        dfc=real(Correlation(zscore(signaldata(t:t+windowsize-1,:)')',1,'corr',windowsize));  % Finding Dynamic FC matrix over window length w
        wdfc(s,i,:)=dfc(uppertriangle);                     % Saving the upper triangular elements of DFC
        i=i+1;t=t+5;                                        % Incrementing Variables
    end                                                 
    X=squeeze(wdfc(s,:,:));                                 % Combing the DFCs for the subject
    [IDX, C, SUMD, D]=kmeans(X,nstates,'Distance','cosine','Replicates',20,'Display','off'); % Performin kmeans to find repeating nstates 
    
    for k=1:nstates                                         % Loop to find the max repeated states
        a(k)=sum(IDX==k);                                   % Saving number of times each states being repeated
    end                                                     
    
    [b,a]=sort(a);                                          % Finding the order of states being repeated maximum in ascending order                                 
    
    Co=C(a(1),:);                                           % Saving the first highest repeated state
    for k=2:nstates                                         % Loop for saving the repeated states for the subject
        Co=horzcat(Co,C(a(k),:));                           % Saving the centroid values for states in ascending order
    end                                                	
    subfc=Co';                                              % Saving the dynamic FC of the subject 
    eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    eval([cell2mat(Subjects(s)) '= subfc' ';']);
end 
clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject C Co D IDX task wdfc X
s=textread('HCP_Subjects_Task.txt');                        % Loads text file with subject ids
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
save('./Data_Save/Pearson_Rest_Task_DynamicFC.mat')         % Save partial FC
 
%% Task Pearson sFC
%  This part of program finds pearson correlation for 321 subjects 3T tfMRI
clear                                                       % Clear workspace
timepoints=[176 253 316 284 232 274 405];                   % Timepoints in each task
tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'};% Task labels
for ta=1:size(tasks,2)                                      % Loop over task
    load([cell2mat(tasks(ta)) '.mat'])                      % Loads timeseries data; the Subjects cell shows the subject name and session name
    uppertriangle=find(triu(ones(90),1));                   % Finding the Upper Diagonal Elements Index 
    for s=1:size(Subjects,1)                                % Loop for saving the primary data for comparison     
        Subject=char(Subjects(s));                          % Loading subject
        eval(['signaldata = ' Subject ';']);                % Extract timeseries data to signaldata       
        signaldata=detrend(signaldata);                     % Normalize the data
        pearsonfc= Correlation(signaldata,1,'corr',timepoints(ta));% Finding Pearson correlation
        subfc=real(pearsonfc(uppertriangle));               % Finding upper triangular part
        eval(['clear ' Subject ';']);                       % Clearing the data which is no longer needed
        eval([cell2mat(Subjects(s)) '= subfc' ';']);        % Saving subject FC
    end 
    clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject  % Clear workspace
    s=textread('HCP_Subjects_Task.txt');                                        % Read task subjects
    k=0;                                                                        % Initialize k
    for i=1:size(Subjects,1)/2                              % Loop over subjects
                a=[2*i-1:2*i];                              % Finding session numbers
                for j=1:2                                   % Loop over sessions
                    Subject=char(Subjects(a(j)));           % Load subject name
                    eval(['ints = ' Subject ';']);          % Save subject FC to ints variable
                    su(:,j)=double(ints);                   % Save subject FC 
                    eval(['clear ' Subject])                % Clear subject variable
                end                                             
                Subj(k+1)=cellstr([ 'Task_' num2str(s(i)) ]);% Adding subject name
                eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,2))./2;' ]);% Find and save task FC 
                k=k+1;                                      % Increment k variable
    end                                                     
    Subjects=[];                                            % Clear subject variable
    Subjects=Subj';                                         % Save new subject names
    clear s su Subj Subject a i ints j k task               % Clear variables
    save(['./Data_Save/' cell2mat(tasks(ta)) '_Pearson_sFC.mat'])% Save pearson sFC
end


%% Task Partial sFC
%  This part of program finds partial correlation for 321 subjects 3T tfMRI
clear                                                       % Clear workspace
timepoints=[176 253 316 284 232 274 405];                   % Timepoints in each task
tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'};% Task labels
for ta=1:size(tasks,2)                                      % Loop over task
    load([cell2mat(tasks(ta)) '.mat'])                      % Loads timeseries data; the Subjects cell shows the subject name and session name
    uppertriangle=find(triu(ones(90),1));                   % Finding the Upper Diagonal Elements Index 
    for s=1:size(Subjects,1)                                % Loop for saving the primary data for comparison     
        Subject=char(Subjects(s));                          % Loading subject
        eval(['signaldata = ' Subject ';']);                % Extract timeseries data to signaldata       
        signaldata=detrend(signaldata);                     % Normalize the data
        pearsonfc= Correlation(signaldata,1,'icov',timepoints(ta));% Finding Partial correlation
        subfc=real(pearsonfc(uppertriangle));               % Finding upper triangular part
        eval(['clear ' Subject ';']);                       % Clearing the data which is no longer needed
        eval([cell2mat(Subjects(s)) '= subfc' ';']);        % Saving subject FC
    end 
    clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject  % Clear workspace
    s=textread('HCP_Subjects_Task.txt');                                        % Read task subjects
    k=0;                                                                        % Initialize k
    for i=1:size(Subjects,1)/2                              % Loop over subjects
                a=[2*i-1:2*i];                              % Finding session numbers
                for j=1:2                                   % Loop over sessions
                    Subject=char(Subjects(a(j)));           % Load subject name
                    eval(['ints = ' Subject ';']);          % Save subject FC to ints variable
                    su(:,j)=double(ints);                   % Save subject FC 
                    eval(['clear ' Subject])                % Clear subject variable
                end                                             
                Subj(k+1)=cellstr([ 'Task_' num2str(s(i)) ]);% Adding subject name
                eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,2))./2;' ]);% Find and save task FC 
                k=k+1;                                      % Increment k variable
    end                                                     
    Subjects=[];                                            % Clear subject variable
    Subjects=Subj';                                         % Save new subject names
    clear s su Subj Subject a i ints j k task               % Clear variables
    save(['./Data_Save/' cell2mat(tasks(ta)) '_Partial_sFC.mat'])% Save partial sFC
end

%% Task Pearson dFC
%  This part of program finds partial correlation for 321 subjects 3T tfMRI
clear                                                       % Clear workspace
timepoints=[176 253 316 284 232 274 405];                   % Timepoints in each task
tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'};% Task labels
for ta=1:size(tasks,2)                                      % Loop over tasks
    load([cell2mat(tasks(ta)) '.mat'])                      % Loads timeseries data; the Subjects cell shows the subject name and session name
    uppertriangle=find(triu(ones(90),1));                   % Finding the Upper Diagonal Elements Index 
    windowsize=85;                                          % Window size
    nstates=4;                                              % Number of states to extract in kmeans
    for s=1:size(Subjects,1)                                % Loop for saving the primary data for comparison     
        t=1;i=1;                                            % Initialiaztion of Variables
        Subject=char(Subjects(s));                          % Loading subject
        eval(['signaldata = ' Subject ';']);                % Extract timeseries data to signaldata
        signaldata=detrend(signaldata);                     % Normalize the data
        while t+windowsize<=timepoints(ta)                  % Loop to perform window analysis
            dfc=Correlation(zscore(signaldata(t:t+windowsize-1,:)')'/sqrt(windowsize),1,'corr',windowsize);  % Finding Dynamic FC matrix over window length w
            wdfc(s,i,:)=dfc(uppertriangle);                 % Saving the upper triangular elements of DFC
            i=i+1;t=t+5;                                    % Incrementing Variables
        end                                                 
        X=squeeze(wdfc(s,:,:));                             % Combing the DFCs for the subject
        [IDX, C, SUMD, D]=kmeans(X,nstates,'Distance','cosine','Replicates',20,'Display','off'); % Performin kmeans to find repeating nstates 

        for k=1:nstates                                     % Loop to find the max repeated states
            a(k)=sum(IDX==k);                               % Saving number of times each states being repeated
        end                                                     

        [b,a]=sort(a);                                      % Finding the order of states being repeated maximum in ascending order                                 

        Co=C(a(1),:);                                       % Saving the first highest repeated state
        for k=2:nstates                                     % Loop for saving the repeated states for the subject
            Co=horzcat(Co,C(a(k),:));                       % Saving the centroid values for states in ascending order
        end                                                	
        subfc=Co';                                          % Saving the dynamic FC of the subject 
        eval(['clear ' Subject ';']);                       % Clearing the data which is no longer needed
        eval([cell2mat(Subjects(s)) '= subfc' ';']);        % Save subject FC
    end 
    clear signaldata subfc uppertriangle TR bfilt afilt Wn pearsonfc s Subject % Clear variables
    s=textread('HCP_Subjects_Task.txt');                    % Read task subjects
    k=0;                                                    % Initialize k=0
    for i=1:size(Subjects,1)/2                              % Loop over subjects
                a=[2*i-1:2*i];                              % Finding session numbers
                for j=1:2                                   % Loop over sessions
                    Subject=char(Subjects(a(j)));           % Load subject name
                    eval(['ints = ' Subject ';']);          % Save subject FC to ints variable
                    su(:,j)=double(ints);                   % Save subject FC 
                    eval(['clear ' Subject])                % Clear subject variable
                end                                             
                Subj(k+1)=cellstr([ 'Task_' num2str(s(i)) ]);% Adding subject name
                eval([cell2mat(Subj(k+1)) '= (su(:,1)+su(:,2))./2;' ]);% Find and save task FC 
                k=k+1;                                      % Increment k variable
    end                                                     
    Subjects=[];                                            % Clear subject variable
    Subjects=Subj';                                         % Save new subject names
    clear s su Subj Subject a i ints j k task               % Clear variables
    save(['./Data_Save/' cell2mat(tasks(ta)) '_Pearson_dFC.mat'])% Save pearson dFC
end
%% Subject identification using pearson sFC
% This part of program performs subject identification using pearson sFC
% from task data
clear                                                       % Clear workspace
for tri=1:1000                                              % Perform 1000 iterations
    load Pearson_Rest_Task_StaticFC.mat                     % Loads timeseries data; the Subjects cell shows the subject name and session name 
    ntrials=size(Subjects,1)/2;                             % Finding number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];refpredsessions=[];% Finding random subjects 
    for i=1:ntrials                                         % Loop  to find correspoding column in subject  and reference prediction
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
        refpredsessions=horzcat(refpredsessions,randsubj(i));
    end                                                         

    basesessions=[0:2:2*ntrials-1]+randi([1 2],1,ntrials);  % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                  % Getting reference column of subject

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['subfc = ' Subject ';']);                         % Extract timeseries data to signaldata       
        comparefc(:,s)=subfc;                                   % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                        
    clear Subjects signaldata                                   % Clear variables

    %% Prediction Task 
    tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'}; % Save task labels
    for ta=1:size(tasks,2)                                                      % Loop over tasks
    load([cell2mat(tasks(ta)) '_Pearson_sFC.mat'])              % Loading task Data
    for i=1:size(refpredsessions,2)                             % Loop to predict task subject
        Subject=char(Subjects(refpredsessions(i)));             % Loading subject                 
        eval(['taskfc = ' Subject ';']);                        % Extract timeseries data to signaldata   
        pcorr=1-pdist2(comparefc',taskfc','cosine');            % finding cosine similarity
        [ecorvalue(i),epred(i)]=max(pcorr);                     % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end
    clear Subjects signaldata                                   % Clearing some variables
    actual=ceil([1:1*ntrials]./1);                              % Actual subject data                                               
    accu=sum(actual==epred);                                    % Comparing prediction with actual data
    eval([cell2mat(tasks(ta)) '(tri)  =  100*accu/i ;']);       % saving accuracy value
    end
end
fprintf('The Emotion identification accuracy for pearson sFC is  %d percent\n', mean(Emotion));     % Output the accuracy in prediction
fprintf('The Gambling identification accuracy for pearson sFC is  %d percent\n', mean(Gambling));   % Output the accuracy in prediction
fprintf('The Language identification accuracy for pearson sFC is  %d percent\n', mean(Language));   % Output the accuracy in prediction
fprintf('The Motor identification accuracy for pearson sFC is  %d percent\n', mean(Motor));         % Output the accuracy in prediction
fprintf('The Relational identification accuracy for pearson sFC is  %d percent\n', mean(Relational)); % Output the accuracy in prediction
fprintf('The Social identification accuracy for pearson sFC is  %d percent\n', mean(Social));       % Output the accuracy in prediction
fprintf('The WM identification accuracy for pearson sFC is  %d percent\n', mean(WM));               % Output the accuracy in prediction
%% Subject identification using partial sFC
% This part of program performs subject identification using partial sFC
% from task data
clear                                                       % Clear workspace
for tri=1:1000                                              % performing 1000 trials
    load Partial_Rest_Task_StaticFC.mat                     % Loads timeseries data; the Subjects cell shows the subject name and session name 
    ntrials=size(Subjects,1)/2;                             % Finding number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];refpredsessions=[];% Finding random subjects 
    for i=1:ntrials                                             % Loop  to find correspoding column in subject  and reference prediction
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
        refpredsessions=horzcat(refpredsessions,randsubj(i));
    end                                                         

    basesessions=[0:2:2*ntrials-1]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['subfc = ' Subject ';']);                         % Extract timeseries data to signaldata       
        comparefc(:,s)=subfc;                                   % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                        
    clear Subjects signaldata                                   % Clear variables

    %% Prediction Task 
    tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'}; % Task labels
    for ta=1:size(tasks,2)                                      % Loop over tasks
    load([cell2mat(tasks(ta)) '_Partial_sFC.mat'])              % Loading task Data
    for i=1:size(refpredsessions,2)                             % Loop to predict task subject
        Subject=char(Subjects(refpredsessions(i)));             % Loading subject                 
        eval(['taskfc = ' Subject ';']);                        % Extract timeseries data to signaldata   
        pcorr=1-pdist2(comparefc',taskfc','cosine');            % Find cosine similarity
        [ecorvalue(i),epred(i)]=max(pcorr);                     % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end
    clear Subjects signaldata                                   % Clearing some variables
    actual=ceil([1:1*ntrials]./1);                              % Actual subject data                                               

    accu=sum(actual==epred);                                    % Comparing prediction with actual data
    eval([cell2mat(tasks(ta)) '(tri)  =  100*accu/i ;']);       % Finding percentage accuracy
    end
end
fprintf('The Emotion identification accuracy for partial sFC is  %d percent\n', mean(Emotion)); % Output the accuracy in prediction
fprintf('The Gambling identification accuracy partial sFC is  %d percent\n', mean(Gambling));   % Output the accuracy in prediction
fprintf('The Language identification accuracy partial sFC is  %d percent\n', mean(Language));   % Output the accuracy in prediction
fprintf('The Motor identification accuracy partial sFC is  %d percent\n', mean(Motor));         % Output the accuracy in prediction
fprintf('The Relational identification accuracy partial sFC is  %d percent\n', mean(Relational)); % Output the accuracy in prediction
fprintf('The Social identification accuracy partial sFC is  %d percent\n', mean(Social));       % Output the accuracy in prediction
fprintf('The WM identification accuracy partial sFC is  %d percent\n', mean(WM));               % Output the accuracy in prediction

%% Subject identification using pearson dFC
% This part of program performs subject identification using pearson dFC
% from task data
clear                                                       % Clear workspace
for tri=1:1000                                              % Perform 1000 trials
    load Pearson_Rest_Task_DynamicFC.mat                    % Loads timeseries data; the Subjects cell shows the subject name and session name 
    ntrials=size(Subjects,1)/2;                             % Find number of subjects
    randsubj=randperm(size(Subjects,1)/2,ntrials);sessiondata=[];refpredsessions=[];% Finding random subjects
    for i=1:ntrials                                         % Loop  to find correspoding column in subject  and reference prediction
        sessiondata=horzcat(sessiondata,[2*randsubj(i)-1:2*randsubj(i)]);
        refpredsessions=horzcat(refpredsessions,randsubj(i));
    end                                                         

    basesessions=[0:2:2*ntrials-1]+randi([1 2],1,ntrials);      % Randomizing subject session baseline for comparison
    refsessions=sessiondata(basesessions);                      % Getting reference column of subject

    %% Saving initial connectivity values of subjects which can be used for prediction 

    for s=1:size(refsessions,2)                                 % Loop for saving the primary data for comparison     
        Subject=char(Subjects(refsessions(s)));                 % Loading subject
        eval(['subfc = ' Subject ';']);                         % Extract timeseries data to signaldata       
        comparefc(:,s)=subfc;                                   % Saving the static connectivity of the subject 
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end                                                        
    clear Subjects signaldata                                   % Clear variables

    %% Prediction  Task 
    tasks={'Emotion','Gambling','Language','Motor','Relational','Social','WM'}; % Task labels
    for ta=1:size(tasks,2)                                      % Looping over tasks
    load([cell2mat(tasks(ta)) '_Pearson_dFC.mat'])              % Loading Task Data
    for i=1:size(refpredsessions,2)                             % Loop to predict task subject
        Subject=char(Subjects(refpredsessions(i)));             % Loading subject                 
        eval(['taskfc = ' Subject ';']);                        % Extract timeseries data to signaldata   
        pcorr=1-pdist2(comparefc',taskfc','cosine');            % Find cosine similarity
        [ecorvalue(i),epred(i)]=max(pcorr);                     % Predicting the subject
        eval(['clear ' Subject ';']);                           % Clearing the data which is no longer needed
    end
    clear Subjects signaldata                                   % Clearing some variables
    actual=ceil([1:1*ntrials]./1);                              % Actual subject data                                               

    accu=sum(actual==epred);                                    % Comparing prediction with actual data
    eval([cell2mat(tasks(ta)) '(tri)  =  100*accu/i ;']);       % Find task accuracy percentage
    end
end
fprintf('The Emotion identification accuracy for pearson dFC is  %d percent\n', mean(Emotion)); % Output the accuracy in prediction
fprintf('The Gambling identification accuracy pearson dFC is  %d percent\n', mean(Gambling));   % Output the accuracy in prediction
fprintf('The Language identification accuracy pearson dFC is  %d percent\n', mean(Language));   % Output the accuracy in prediction
fprintf('The Motor identification accuracy pearson dFC is  %d percent\n', mean(Motor));         % Output the accuracy in prediction
fprintf('The Relational identification accuracy pearson dFC is  %d percent\n', mean(Relational)); % Output the accuracy in prediction
fprintf('The Social identification accuracy pearson dFC is  %d percent\n', mean(Social));       % Output the accuracy in prediction
fprintf('The WM identification accuracy pearson dFC is  %d percent\n', mean(WM));               % Output the accuracy in prediction
