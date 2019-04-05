% sex Identification
% Menon, S. S., & Krishnamurthy, K. (2019). A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Biological Sex Using Intrinsic Individual Brain Connectivity. Scientific Reports, 9(1), 5729. https://doi.org/10.1038/s41598-019-42090-4
% This program can be used to reproduce results explained in sex
% identification ( Results section 4). A leave one out strategy has been
% implemened in thi program
%% Pearson sFC
% This section of program identifies sex of a subject using Pearson sFC
% network measure
clear                                                       % Clear workspace
for trial=1:1000                                            % Perform 1000 trials
ntrials=343;                                                % Number of subjects
subjectsex=textread('MF343Subjects.txt');                   % Reading subject information whether male (1) or female (2)
load Pearson_Rest_StaticFC.mat                              % Loads timeseries data; the Subjects cell shows the subject name and session name 
malesubjects=find(subjectsex(:,2)==1)';                     % Male Subjects
femalesubjects=find(subjectsex(:,2)==2)';                   % Female Subjects
randmales=randperm(138,138);                                % Finding random males
randfemales=randperm(205,205);                              % Finding random females
randsubject=[malesubjects(randmales) femalesubjects(randfemales)];  % Finding random subjects 
randsubsession=2*randsubject-randi([0 1],1,343);            % Finding randon sessions as base


% Saving FC by leaving a subject (if 50 leave 1st and use 49 to find representative male/female fc)

for trials=1:ntrials                                        % Loop to leave a subject and find FC
    maletimeseries=0;femaletimeseries=0;m=0;f=0;            % Initializing timeseries for male and female
    for s=1:size(randsubsession,2)                          % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(randsubsession(s)));          % Loading session name                    
        eval(['fcdata = ' Subject ';']);                    % Loading session data           
        if randsubsession(s)~=2*ntrials-1 && randsubsession(s)~=2*ntrials % Leaving out a subject
        if subjectsex(randsubject(s),2)==1                  % Checking whether subject is male or not
            maletimeseries=maletimeseries+fcdata;           % Concat male timeseries
            m=m+1;                                          % Incrementing male variable
        else
            femaletimeseries=femaletimeseries+fcdata;       % Concat female timeseries
            f=f+1;                                          % Incrementing female varible
        end  
        end
        clear fcdata  Subject                               % Clearing variable
    end    
    standardfc(trials,:,1)= maletimeseries/m;               % Saving the male fc to comparison matrix
    standardfc(trials,:,2)=femaletimeseries/f;              % Saving the female fc to comparison matrix  
end

% Finding the sex for left out subject 

accu(1:2:686)=subjectsex(:,2)';accu(2:2:686)=subjectsex(:,2)';  % Actual sex for the subjects
for trials=1:686                                                % Loop to find sex of left out of subject
    Subject=char(Subjects(trials));                             % Loading subject name and session for left out                                   
    eval(['fcdata = ' Subject ';']);                            % Loading the timeseries                             
    pcorr=1-pdist2(squeeze(standardfc(ceil(trials/2),:,:))',fcdata','cosine');% Correlationg left out FC with male/female representative FC                           
    [corvalue(trials),pred(trials)]=max(pcorr);                 % Predicting the subject                                                 
    clear signaldata                                            % Clearing variable
end

% Prediction comparison

k(1,:)=accu;k(2,:)=pred;                                        % Finding number of male female predictions
m=sum(k(1,:)==1&k(2,:)==1);                                     % males correctly classified
f=sum(k(1,:)==2&k(2,:)==2);                                     % female correctly classified
number=sum(accu==pred);                                         % number of accurate predictions
b(trial)=(m+f)/686;                                             % Finding percentage accuracy
end
fprintf('The mean accuracy for sex identification by  Pearson sFC is %d \n', 100*mean(b));% The number of males females predicted
%% Partial sFC
% This section of program identifies sex of a subject using Partial sFC
% network measure
clear                                                       % Clear workspace
for trial=1:10
ntrials=343;                                                % Number of subjects
subjectsex=textread('MF343Subjects.txt');                   % Reading subject information whether male (1) or female (2)
load Partial_Rest_StaticFC.mat                              % Loads timeseries data; the Subjects cell shows the subject name and session name 
malesubjects=find(subjectsex(:,2)==1)';                     % Male Subjects
femalesubjects=find(subjectsex(:,2)==2)';                   % Female Subjects
randmales=randperm(138,138);                                % Finding random males
randfemales=randperm(205,205);                              % Finding random females
randsubject=[malesubjects(randmales) femalesubjects(randfemales)];  % Finding random subjects from 200
randsubsession=2*randsubject-randi([0 1],1,343);            % Finding randon sessions as base


% Saving FC by leaving a subject (if 50 leave 1st and use 49 to find representative male/female fc)

for trials=1:ntrials                                            % Loop to leave a subject and find FC
    maletimeseries=0;femaletimeseries=0;m=0;f=0;                % Initializing timeseries for male and female
    for s=1:size(randsubsession,2)                              % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(randsubsession(s)));              % Loading session name                    
        eval(['fcdata = ' Subject ';']);                        % Loading session data           
        if randsubsession(s)~=2*ntrials-1 && randsubsession(s)~=2*ntrials % Leave out subject
        if subjectsex(randsubject(s),2)==1                      % Checking whether subject is male or not
            maletimeseries=maletimeseries+fcdata;               % Concat male timeseries
            m=m+1;                                              % increment male variable
        else 
            femaletimeseries=femaletimeseries+fcdata;           % Concat female timeseries
            f=f+1;                                              % Increment female variable
        end  
        end
        clear fcdata  Subject                                   % Clearing variable
    end    
    standardfc(trials,:,1)= maletimeseries/m;                   % Saving the male fc to comparison matrix
    standardfc(trials,:,2)=femaletimeseries/f;                  % Saving the female fc to comparison matrix  
end

% Finding the sex for left out subject 
accu(1:2:686)=subjectsex(:,2)';accu(2:2:686)=subjectsex(:,2)';  % Actual sex for the subjects
for trials=1:686                                                % Loop to find sex of left out of subject
    Subject=char(Subjects(trials));                             % Loading subject name and session for left out                                   
    eval(['fcdata = ' Subject ';']);                            % Loading the timeseries                             
    pcorr=1-pdist2(squeeze(standardfc(ceil(trials/2),:,:))',fcdata','cosine');               % Correlationg left out FC with male/female representative FC                           
    [corvalue(trials),pred(trials)]=max(pcorr);                 % Predicting the subject                                                 
    clear signaldata                                            % Clearing variable
end

% Prediction comparison
k(1,:)=accu;k(2,:)=pred;                                        % Finding number of male female predictions
m=sum(k(1,:)==1&k(2,:)==1);                                     % males correctly classified
f=sum(k(1,:)==2&k(2,:)==2);                                     % female correctly classified
number=sum(accu==pred);                                         % number of accurate predictions
b(trial)=(m+f)/686;                                             % Finding percentage accuracy
end
fprintf('The mean accuracy for sex identification by  Partial sFC is %d \n', 100*mean(b));% The number of males females predicted

%% Dynamic FC
% This section of program identifies sex of a subject using Pearson dFC
% network measure
clear                                                       % Clear workspace
for trial =1:10
ntrials=343;                                                % number of subjects
subjectsex=textread('MF343Subjects.txt');                   % Reading subject information whether male (1) or female (2)
load Pearson_Rest_DynamicFC.mat                             % Loads timeseries data; the Subjects cell shows the subject name and session name 
malesubjects=find(subjectsex(:,2)==1)';                     % Male Subjects
femalesubjects=find(subjectsex(:,2)==2)';                   % Female Subjects
randmales=randperm(138,138);                                % Finding random males
randfemales=randperm(205,205);                              % Finding random females
randsubject=[malesubjects(randmales) femalesubjects(randfemales)];    % Finding random subjects from 200
randsubsession=2*randsubject-randi([0 1],1,343);            % Finding randon sessions as base


% Saving FC by leaving a subject (if 50 leave 1st and use 49 to find representative male/female fc)

for trials=1:ntrials                                        % Loop to leave a subject and find FC
    maletimeseries=0;femaletimeseries=0;m=0;f=0;            % Initializing timeseries for male and female
    for s=1:size(randsubsession,2)                          % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(randsubsession(s)));          % Loading session name                    
        eval(['fcdata = ' Subject ';']);                    % Loading session data           
        if randsubsession(s)~=2*ntrials-1 && randsubsession(s)~=2*ntrials % Leave out subject
        if subjectsex(randsubject(s),2)==1                  % Checking whether subject is male or not
            maletimeseries=maletimeseries+fcdata;           % Concat male timeseries
            m=m+1;                                          % Increment male variable
        else
            femaletimeseries=femaletimeseries+fcdata;       % Concat female timeseries
            f=f+1;                                          % Incremant female variable    
        end  
        end
        clear fcdata  Subject                               % Clearing variable
    end    
    standardfc(trials,:,1)= maletimeseries/m;               % Saving the male fc to comparison matrix
    standardfc(trials,:,2)=femaletimeseries/f;              % Saving the male fc to comparison matrix  
end

% Finding the sex for left out subject 
accu(1:2:686)=subjectsex(:,2)';accu(2:2:686)=subjectsex(:,2)';% Actual sex for the subjects

for trials=1:686                                            % Loop to find sex of left out of subject
    Subject=char(Subjects(trials));                         % Loading subject name and session for left out                                   
    eval(['fcdata = ' Subject ';']);                        % Loading the timeseries                             
    pcorr=1-pdist2(squeeze(standardfc(ceil(trials/2),:,:))',fcdata','cosine');               % Correlationg left out FC with male/female representative FC                           
    [corvalue(trials),pred(trials)]=max(pcorr);             % Predicting the subject                                                 
    clear signaldata                                        % Clearing variable
end

% Prediction comparison
k(1,:)=accu;k(2,:)=pred;                                    % Finding number of male female predictions
m=sum(k(1,:)==1&k(2,:)==1);                                 % males correctly classified
f=sum(k(1,:)==2&k(2,:)==2);                                 % female correctly classified
number=sum(accu==pred);                                     % number of accurate predictions
be(trial)=(m+f)/686;                                        % Finding percentage accuracy
end
fprintf('The mean accuracy for sex identification by Pearson dFC is %d \n', 100*mean(be));% The number of males females predicted
