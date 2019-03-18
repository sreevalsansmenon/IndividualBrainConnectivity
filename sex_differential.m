% sex Identification differential power of edges
% This program can be used to reproduce results explained in differential power of edges in sex
% identification ( Results section 5).
% A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Biological Sex Using Intrinsic Individual Brain Connectivity
% Sreevalsan S Menon, K Krishnamurthy

%% Pearson sFC
% This section of program finds edges differential power in sex
% identification using pearson sFC
clear                                               % Clear workspace
subjectsex=textread('MF343Subjects.txt');           % Reading subject information whether male (1) or female (2)
load Pearson_Rest_StaticFC.mat                      % Loads timeseries data; the Subjects cell shows the subject name and session name 
intm=[];intf=[];m=0;f=0;                            % Initializing timeseries for male and female
for s=1:686                                         % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            intm=[intm fcdata];                     % Concat male timeseries
            m=m+1;                                  % Increment male variable
        else
            intf=[intf fcdata];                     % Concat female timeseries
            f=f+1;                                  % Increment male variable           
        end     
        clear fcdata  Subject                       % Clearing variable
 end    
femalefc=mean(intf,2);                              % Mean female FC 
malefc=mean(intm,2);                                % Mean male FC

% Probability calculation
    
    mama=[];mafe=[];fema=[];fefe=[];m=0;f=0;mprob=zeros(4005,1);fprob=zeros(4005,1);% Initializing variables
    for s=1:686                                     % Loop to find phi for male/female                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            mama=fcdata.*malefc;                    % male male product value
            mafe=fcdata.*femalefc;                  % Male female product value
            mprob=mprob+double(mama>mafe);          % Probability male edge> male female
            m=m+1;                                  % Incrementing male variable
        else
            fefe=fcdata.*femalefc;                  % female female product value
            fema=fcdata.*malefc;                    % Male female product value
            fprob=fprob+double(fefe>fema);          % Probability female edge> male female
            f=f+1;                                  % Incremant female variable                
        end     
        clear fcdata  Subject                       % Clearing variable
    end    
fprob=fprob/f;mprob=mprob/m;                        % Finding male and female probabilities

% Female differential power
temp=fprob;                                         % Adding differential power to a temporary variable
y=prctile(temp,95);                                 % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
I=find(triu(ones(90),1));                           % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Pearson sFC - Female')        % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Male differntal power
temp=mprob;                                         % Adding differential power to a temporary variable
y=prctile(mprob,95);                                % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Pearson sFC - Male')          % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

%% Partial sFC
% This section of program finds edges differential power in sex
% identification using partial sFC
clear                                               % Clear workspace
subjectsex=textread('MF343Subjects.txt');           % Reading subject information whether male (1) or female (2)
load Partial_Rest_StaticFC.mat                      % Loads timeseries data; the Subjects cell shows the subject name and session name 
intm=[];intf=[];m=0;f=0;                            % Initializing timeseries for male and female
for s=1:686                                         % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            intm=[intm fcdata];                     % Concat male timeseries
            m=m+1;                                  % Increment male variable
        else
            intf=[intf fcdata];                     % Concat female timeseries
            f=f+1;                                  % Increment male variable           
        end     
        clear fcdata  Subject                       % Clearing variable
 end    
femalefc=mean(intf,2);                              % Mean female FC 
malefc=mean(intm,2);                                % Mean male FC

% Probability calculation
    
mama=[];mafe=[];fema=[];fefe=[];m=0;f=0;mprob=zeros(4005,1);fprob=zeros(4005,1);% Initializing variables
    for s=1:686                                     % Loop to find phi for male/female                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            mama=fcdata.*malefc;                    % male male product value
            mafe=fcdata.*femalefc;                  % Male female product value
            mprob=mprob+double(mama>mafe);          % Probability male edge> male female
            m=m+1;                                  % Incrementing male variable
        else
            fefe=fcdata.*femalefc;                  % female female product value
            fema=fcdata.*malefc;                    % Male female product value
            fprob=fprob+double(fefe>fema);          % Probability female edge> male female
            f=f+1;                                  % Incremant female variable                
        end     
        clear fcdata  Subject                       % Clearing variable
    end    
fprob=fprob/f;mprob=mprob/m;                        % Finding male and female probabilities

% Female differential power
temp=fprob;                                         % Adding differential power to a temporary variable
y=prctile(temp,95);                                 % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
I=find(triu(ones(90),1));                           % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Partial sFC - Female')        % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Male differntal power
temp=mprob;                                         % Adding differential power to a temporary variable
y=prctile(mprob,95);                                % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Partial sFC - Male')          % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

%% Dynamic FC
% This section of program finds edges differential power in sex
% identification using pearson dFC
clear                                               % Clear workspace
subjectsex=textread('MF343Subjects.txt');           % Reading subject information whether male (1) or female (2)
load Pearson_Rest_DynamicFC.mat                     % Loads timeseries data; the Subjects cell shows the subject name and session name 
intm=[];intf=[];m=0;f=0;                            % Initializing timeseries for male and female
for s=1:686                                         % Loop to find FC for ntrial-1 subjects                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            intm=[intm fcdata];                     % Concat male timeseries
            m=m+1;                                  % Increment male variable
        else
            intf=[intf fcdata];                     % Concat female timeseries
            f=f+1;                                  % Increment male variable           
        end     
        clear fcdata  Subject                       % Clearing variable
 end    
femalefc=mean(intf,2);                              % Mean female FC 
malefc=mean(intm,2);                                % Mean male FC

% Probability calculation
    
mama=[];mafe=[];fema=[];fefe=[];m=0;f=0;mprob=zeros(4005*4,1);fprob=zeros(4005*4,1);% Initializing variables
    for s=1:686                                     % Loop to find phi for male/female                              
        Subject=char(Subjects(s));                  % Loading session name                    
        eval(['fcdata = ' Subject ';']);            % Loading session data           
        if subjectsex(ceil(s/2),2)==1               % Checking whether subject is male or not
            mama=fcdata.*malefc;                    % male male product value
            mafe=fcdata.*femalefc;                  % Male female product value
            mprob=mprob+double(mama>mafe);          % Probability male edge> male female
            m=m+1;                                  % Incrementing male variable
        else
            fefe=fcdata.*femalefc;                  % female female product value
            fema=fcdata.*malefc;                    % Male female product value
            fprob=fprob+double(fefe>fema);          % Probability female edge> male female
            f=f+1;                                  % Incremant female variable                
        end     
        clear fcdata  Subject                       % Clearing variable
    end    
fprob=fprob/f;mprob=mprob/m;                        % Finding male female edge probabilities
temp=fprob;fprob=[];fprob=temp(1:4005)+temp(4006:8010)+temp(8011:12015)+temp(12016:16020);fprob=fprob/4; % Finding probability for edges female
temp=mprob;mprob=[];mprob=temp(1:4005)+temp(4006:8010)+temp(8011:12015)+temp(12016:16020);mprob=mprob/4; % Finding probability for edges male

% Female differential power
temp=fprob;                                         % Adding differential power to a temporary variable
y=prctile(temp,95);                                 % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
I=find(triu(ones(90),1));                           % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Pearson dFC - Female')        % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Male differntal power
temp=mprob;                                         % Adding differential power to a temporary variable
y=prctile(mprob,95);                                % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                     % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                         % Initializing network variable
nnzb=nnz(temp)*2;                                   % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);         % Finding percentage edges in network 1
network(2,2)=sum(sum(df1(4:8,4:8)))/(nnzb);      network(2,1)=sum(sum(df1(1:3,4:8)))/nnzb;% Finding percentage edges in network 2
network(3,3)=sum(sum(df1(9:17,9:17)))/(nnzb);    network(3,1)=sum(sum(df1(1:3,9:17)))/nnzb;     network(3,2)=sum(sum(df1(4:8,9:17)))/nnzb;% Finding percentage edges in network 3
network(4,4)=sum(sum(df1(18:19,18:19)))/(nnzb);  network(4,1)=sum(sum(df1(1:3,18:19)))/nnzb;    network(4,2)=sum(sum(df1(4:8,18:19)))/nnzb;    network(4,3)=sum(sum(df1(9:17,18:19)))/nnzb;% Finding percentage edges in network 4
network(5,5)=sum(sum(df1(20:26,20:26)))/(nnzb);  network(5,1)=sum(sum(df1(1:3,20:26)))/nnzb;    network(5,2)=sum(sum(df1(4:8,20:26)))/nnzb;    network(5,3)=sum(sum(df1(9:17,20:26)))/nnzb;   network(5,4)=sum(sum(df1(18:19,20:26)))/nnzb;% Finding percentage edges in network 5
network(6,6)=sum(sum(df1(27:32,27:32)))/(nnzb);  network(6,1)=sum(sum(df1(1:3,27:32)))/nnzb;    network(6,2)=sum(sum(df1(4:8,27:32)))/nnzb;    network(6,3)=sum(sum(df1(9:17,27:32)))/nnzb;   network(6,4)=sum(sum(df1(18:19,27:32)))/nnzb;  network(6,5)=sum(sum(df1(20:26,27:32)))/nnzb;% Finding percentage edges in network 6
network(7,7)=sum(sum(df1(33:38,33:38)))/(nnzb);  network(7,1)=sum(sum(df1(1:3,33:38)))/nnzb;    network(7,2)=sum(sum(df1(4:8,33:38)))/nnzb;    network(7,3)=sum(sum(df1(9:17,33:38)))/nnzb;   network(7,4)=sum(sum(df1(18:19,33:38)))/nnzb;  network(7,5)=sum(sum(df1(20:26,33:38)))/nnzb;  network(7,6)=sum(sum(df1(27:32,33:38)))/nnzb;% Finding percentage edges in network 7
network(8,8)=sum(sum(df1(39:50,39:50)))/(nnzb);  network(8,1)=sum(sum(df1(1:3,39:50)))/nnzb;    network(8,2)=sum(sum(df1(4:8,39:50)))/nnzb;    network(8,3)=sum(sum(df1(9:17,39:50)))/nnzb;   network(8,4)=sum(sum(df1(18:19,39:50)))/nnzb;  network(8,5)=sum(sum(df1(20:26,39:50)))/nnzb;  network(8,6)=sum(sum(df1(27:32,39:50)))/nnzb;  network(8,7)=sum(sum(df1(33:38,39:50)))/nnzb;% Finding percentage edges in network 8
network(9,9)=sum(sum(df1(51:54,51:54)))/(nnzb);  network(9,1)=sum(sum(df1(1:3,51:54)))/nnzb;    network(9,2)=sum(sum(df1(4:8,51:54)))/nnzb;    network(9,3)=sum(sum(df1(9:17,51:54)))/nnzb;   network(9,4)=sum(sum(df1(18:19,51:54)))/nnzb;  network(9,5)=sum(sum(df1(20:26,51:54)))/nnzb;  network(9,6)=sum(sum(df1(27:32,51:54)))/nnzb;  network(9,7)=sum(sum(df1(33:38,51:54)))/nnzb;  network(9,8)=sum(sum(df1(39:50,51:54)))/nnzb;% Finding percentage edges in network 9
network(10,10)=sum(sum(df1(55:56,55:56)))/(nnzb); network(10,1)=sum(sum(df1(1:3,55:56)))/nnzb;   network(10,2)=sum(sum(df1(4:8,55:56)))/nnzb;   network(10,3)=sum(sum(df1(9:17,55:56)))/nnzb;  network(10,4)=sum(sum(df1(18:19,55:56)))/nnzb; network(10,5)=sum(sum(df1(20:26,55:56)))/nnzb; network(10,6)=sum(sum(df1(27:32,55:56)))/nnzb; network(10,7)=sum(sum(df1(33:38,55:56)))/nnzb; network(10,8)=sum(sum(df1(39:50,55:56)))/nnzb; network(10,9)=sum(sum(df1(51:54,55:56)))/nnzb;% Finding percentage edges in network 10
network(11,11)=sum(sum(df1(57:62,57:62)))/(nnzb); network(11,1)=sum(sum(df1(1:3,57:62)))/nnzb;   network(11,2)=sum(sum(df1(4:8,57:62)))/nnzb;   network(11,3)=sum(sum(df1(9:17,57:62)))/nnzb;  network(11,4)=sum(sum(df1(18:19,57:62)))/nnzb; network(11,5)=sum(sum(df1(20:26,57:62)))/nnzb; network(11,6)=sum(sum(df1(27:32,57:62)))/nnzb; network(11,7)=sum(sum(df1(33:38,57:62)))/nnzb; network(11,8)=sum(sum(df1(39:50,57:62)))/nnzb; network(11,9)=sum(sum(df1(51:54,57:62)))/nnzb;    network(11,10)=sum(sum(df1(55:56,57:62)))/nnzb;% Finding percentage edges in network 11
network(12,12)=sum(sum(df1(63:69,63:69)))/(nnzb); network(12,1)=sum(sum(df1(1:3,63:69)))/nnzb;   network(12,2)=sum(sum(df1(4:8,63:69)))/nnzb;   network(12,3)=sum(sum(df1(9:17,63:69)))/nnzb;  network(12,4)=sum(sum(df1(18:19,63:69)))/nnzb; network(12,5)=sum(sum(df1(20:26,63:69)))/nnzb; network(12,6)=sum(sum(df1(27:32,63:69)))/nnzb; network(12,7)=sum(sum(df1(33:38,63:69)))/nnzb; network(12,8)=sum(sum(df1(39:50,63:69)))/nnzb; network(12,9)=sum(sum(df1(51:54,63:69)))/nnzb;    network(12,10)=sum(sum(df1(55:56,63:69)))/nnzb;    network(12,11)=sum(sum(df1(57:62,63:69)))/nnzb;% Finding percentage edges in network 12
network(13,13)=sum(sum(df1(70:79,70:79)))/(nnzb); network(13,1)=sum(sum(df1(1:3,70:79)))/nnzb;   network(13,2)=sum(sum(df1(4:8,70:79)))/nnzb;   network(13,3)=sum(sum(df1(9:17,70:79)))/nnzb;  network(13,4)=sum(sum(df1(18:19,70:79)))/nnzb; network(13,5)=sum(sum(df1(20:26,70:79)))/nnzb; network(13,6)=sum(sum(df1(27:32,70:79)))/nnzb; network(13,7)=sum(sum(df1(33:38,70:79)))/nnzb; network(13,8)=sum(sum(df1(39:50,70:79)))/nnzb; network(13,9)=sum(sum(df1(51:54,70:79)))/nnzb;    network(13,10)=sum(sum(df1(55:56,70:79)))/nnzb;    network(13,11)=sum(sum(df1(57:62,70:79)))/nnzb;    network(13,12)=sum(sum(df1(63:69,70:79)))/nnzb;% Finding percentage edges in network 13
network(14,14)=sum(sum(df1(80:90,80:90)))/(nnzb); network(14,1)=sum(sum(df1(1:3,80:90)))/nnzb;   network(14,2)=sum(sum(df1(4:8,80:90)))/nnzb;   network(14,3)=sum(sum(df1(9:17,80:90)))/nnzb;  network(14,4)=sum(sum(df1(18:19,80:90)))/nnzb; network(14,5)=sum(sum(df1(20:26,80:90)))/nnzb; network(14,6)=sum(sum(df1(27:32,80:90)))/nnzb; network(14,7)=sum(sum(df1(33:38,80:90)))/nnzb; network(14,8)=sum(sum(df1(39:50,80:90)))/nnzb; network(14,9)=sum(sum(df1(51:54,80:90)))/nnzb;    network(14,10)=sum(sum(df1(55:56,80:90)))/nnzb;    network(14,11)=sum(sum(df1(57:62,80:90)))/nnzb;    network(14,12)=sum(sum(df1(63:69,80:90)))/nnzb; network(14,13)=sum(sum(df1(70:79,80:90)))/nnzb;% Finding percentage edges in network 14
df1=tril(network)+tril(network,-1)';df1=df1.*100;% Transferring lower part to upper triangular part and multiplying by 100 for percentage

% Plot differential power percentages within and between netowrks 
figure;imagesc(df1);colormap(jet);axis image;colorbar;  % Plot edge contribution percentages
title('Differential power Pearson dFC - Male')          % Add title to the figure
ax = gca;ax.XTick =1:14; ax.YTick =1:14;                % Add tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};         % Create x tick labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};   % Create y tick labels
set(gca, 'FontSize', 15);xtickangle(90);                % Set label font and orientation
textStrings = num2str(df1(:), '%0.2f');                 % Create text strings 
textStrings = strtrim(cellstr(textStrings));            % Remove white spaces
[x, y] = meshgrid(1:14);                                % Create mesh grid for RSNs
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','Fontsize',15);% Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end
                                                                
