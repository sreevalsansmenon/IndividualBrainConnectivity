%  Edge analysis in subject identification
% This program  performs edge analysis (Results section 2).  Note this has
% to be ran after resting_state_fingerprint function
% Menon, S. S., & Krishnamurthy, K. (2019). A Comparison of Static and Dynamic Functional Connectivities for Identifying Subjects and Biological Sex Using Intrinsic Individual Brain Connectivity. Scientific Reports, 9(1), 5729. https://doi.org/10.1038/s41598-019-42090-4

%% Pearson sFC edge consistency and variation
% This section of program calculates edge consistency and variation for
% Pearson sFC
clear                                           % Clear workspace
load Pearson_Rest_StaticFC.mat                  % Load pearson sFC 
for i=1:686                                     % Loop over available session data
        Subject=char(Subjects(i));              % Finding subject                        
        eval(['ints = ' Subject ';']);          % Storing subject FC value
        su(:,i)=ints;                           % Concatinating subject FCs
end                                             % Loop terminates                          
st=std(su,0,2);                                 % Finding standard deviation of edges

% Finding Consistency 
temp=st;                                        % Storing variation to a temperorary variable 
y=prctile(st,5);                                % Finding 5 percentile for consistency
temp(temp>y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
dfcx1=zeros(90,90); dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Initializing and finding the active edges in consistency
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the consistency 
figure;                                         % Opens new figure to plot consistency
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting consistency 
title('Pearson sFC edge consistency')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};     % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Finding Variability
temp=st;                                        % storing variation to a temperorary variable 
y=prctile(st,95);                               % Finding 95 percentile for variability
temp(temp<y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
dfcx1=zeros(90,90); dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';              % Initializing and finding the active edges in variability
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the variability 
figure;                                         % Opens new figure to plot variability
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting variability
title('Pearson sFC edge variability')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};     % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

%% Partial sFC edge consistency and variation
% This section of program calculates edge consistency and variation for
% Partial sFC
clear                                           % Clear workspace
load Partial_Rest_StaticFC.mat                  % Load pearson sFC 
for i=1:686                                     % Loop over available session data
        Subject=char(Subjects(i));              % Finding subject                        
        eval(['ints = ' Subject ';']);          % Storing subject FC value
        su(:,i)=ints;                           % Concatinating subject FCs
end                                             % Loop terminates                          
st=std(su,0,2);                                 % Finding standard deviation of edges

% Finding Consistency 
temp=st;                                        % Storing variation to a temperorary variable 
y=prctile(st,5);                                % Finding 5 percentile for consistency
temp(temp>y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
dfcx1=zeros(90,90); dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Initializing and finding the active edges in consistency
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the consistency 
figure;                                         % Opens new figure to plot consistency
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting consistency 
title('Partial sFC edge consistency')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};      % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Finding Variability
temp=st;                                        % storing variation to a temperorary variable 
y=prctile(st,95);                               % Finding 95 percentile for variability
temp(temp<y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
dfcx1=zeros(90,90); dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';              % Initializing and finding the active edges in variability
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the variability 
figure;                                         % Opens new figure to plot variability
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting variability
title('Partial sFC edge variability')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};     % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

%% Dynamic dFC edge consistency and variation
% This section of program calculates edge consistency and variation for
% Pearson dFC
clear                                           % Clear workspace
load Pearson_Rest_DynamicFC.mat                 % Load pearson sFC 
for i=1:686                                     % Loop over available session data
        Subject=char(Subjects(i));              % Finding subject                        
        eval(['ints = ' Subject ';']);          % Storing subject FC value
        su(:,i)=ints;                           % Concatinating subject FCs
end                                             % Loop terminates                          
st=std(su,0,2);                                 % Finding standard deviation of edges

% Finding Consistency 
temp=st;                                        % Storing variation to a temperorary variable 
y=prctile(st,5);                                % Finding 5 percentile for consistency
temp(temp>y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
afc=temp(1:4005)+temp(4006:8010)+temp(8011:12015)+temp(12016:16020);afc(afc~=0)=1;dfcx1=zeros(90,90);dfcx1(I)=afc;df1=triu(dfcx1)+triu(dfcx1,1)'; % Initializing and finding the active edges in consistency
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the consistency 
figure;                                         % Opens new figure to plot consistency
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting consistency 
title('Pearson dFC edge consistency')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};     % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

% Finding Variability
temp=st;                                        % storing variation to a temperorary variable 
y=prctile(st,95);                               % Finding 95 percentile for variability
temp(temp<y)=0;temp(temp~=0)=1;                 % Thresholding for percentile value
I=find(triu(ones(90),1));                       % Finding upper triangular elements
afc=temp(1:4005)+temp(4006:8010)+temp(8011:12015)+temp(12016:16020);afc(afc~=0)=1;dfcx1=zeros(90,90);dfcx1(I)=afc;df1=triu(dfcx1)+triu(dfcx1,1)';% Initializing and finding the active edges in variability
% Following codes find percentage number of edges in and between networks diagonal will be same so only lower triangular part is calculated here              
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
% Plotting the variability 
figure;                                         % Opens new figure to plot variability
imagesc(df1);colormap(jet);axis image; colorbar;% Plotting variability
title('Pearson dFC edge variability')           % Adding title to the plot
ax = gca;ax.XTick =1:14;ax.YTick =1:14;         % Adding row column tick marks
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRE','V1','RECN','ASAL','VDMN','VISUO'};     % Adding network labels
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding network labels
set(gca, 'FontSize', 15);xtickangle(90);        % Setting font size
textStrings = num2str(df1(:), '%0.2f');         % Create string for percentages 
textStrings = strtrim(cellstr(textStrings));    % Remove white spaces
[x, y] = meshgrid(1:14);                        % Create RSN meshgrid
hStrings = text(x(:), y(:), textStrings(:),'HorizontalAlignment', 'center','FontSize',15); % Add text string to plots
midValue = mean(get(gca, 'CLim'));              % Get positions
textColors = repmat(df1(:) > midValue, 1, 3);set(hStrings,'Color','w'); % Set colors
hold on;                                        % Hold the plot
for i = 1:15                                    % Plotting RSN line grids
   plot([.5,14.5],[i-.5,i-.5],'w-');            % Plot lines
   plot([i-.5,i-.5],[.5,14.5],'w-');            % Plot lines
end

%% Pearson sFC differential power
% This part of program finds differential power for Pearson sFC
clear                                                   % Clear workspace
load Pearson_Rest_StaticFC.mat                          % Load pearson sFC

% Finding Phi and PhiIntra
for i=1:343                                             % Loop over subjects
        Subject=char(Subjects(2*i-1));                  % Find subject day 1 FC                            
        eval(['ints = ' Subject ';']);                  % Save day 1 fc to a variable
        su(:,i)=ints;                                   % Save fc to variable su
        for j=1:343                                     % Loop over subjects to find day 2 session
            Subject=char(Subjects(2*j));                % Load subject day 2 FC                           
            eval(['ints = ' Subject ';']);              % Save fc to a variable
            su(:,j)=ints;                               % Saving value to su
            phi(:,i,j) = su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));        % Finding phi value
            if i==j                                                             % Finding whether subject is same
                phiintra(:,i)=su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));   % Saving phi value for same subject
            end                         
        end 
end
GroupPhi=mean(phiintra,2);                              % Finding mean phi for subjects

% Probability calculation
% Loop below finds probability of phi of a subject greter than between
% subjects
for i=1:343                                             % Looping over subjects 
    b=1:343;                                            % Looping over subjects 
    c=setdiff(b,i);                                     % Remove subject i
    edge(:,i)=zeros(4005,1);                            % Initializing edge variable
    for j=1:size(c,2)                                   % Looping over subjects removing subject i
        k=phiintra(:,i)>squeeze(phi(:,i,c(j)));         % Finding phi(i,i)>phi(i,j)
        l=phiintra(:,i)>squeeze(phi(:,c(j),i));         % Finding phi(i,i)>phi(j,i)
        edge(:,i)=edge(:,i)+k+l;                        % Finding number of times edge was having higher phi
    end
    p(:,i)=edge(:,i)./(2*342);                          % Finding probabilty  
end
avgp=mean(p,2);                                         % Finding differential power
I=find(triu(ones(90),1));                               % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=avgp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding the diferential power in matrix format

% Plot differential power values for edges

figure;imagesc(df1);c=jet;colormap(c);axis image;colorbar;                      % Ploting differential power in new window
title('Differential power Pearson sFC')                                         % Adding title to image
ax = gca;ax.XTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5]; % Adding tick marks for label
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
set(gca,'XTickLabelRotation',90,'fontsize',13)                                  % Seting fontsize and rotation
ax.YTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5];          % Adding tickmarks for label
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
l1=[3.5 3.5;8.5 8.5 ; 17.5  17.5 ;19.5  19.5 ;26.5  26.5 ;32.5  32.5 ;38.5  38.5 ;50.5  50.5 ;54.5  54.5 ;56.5  56.5 ;62.5  62.5 ;69.5  69.5 ;79.5  79.5]; % Creating line centroids
l2=[0.5 90.5;0.5 90.5; 0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5 ];                               % Creating line centroids

hold on                                                 % Hold figure
for i= 1:13                                             % Loop to plot lines
    line(l1(i,:),l2(i,:),'color','w')                   % Plotting line             
    line(l2(i,:),l1(i,:),'color','w')                   % Plotting lines
end 

% Probability Percentile 

temp=avgp;                                              % Adding differential power to a temporary variable
y=prctile(avgp,95);                                     % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                         % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
title('Edge contribution Pearson sFC')                  % Add title to the figure
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

%% Partial sFC differential power
% This part of program finds differential power for Partial sFC
clear                                                   % Clear workspace
load Partial_Rest_StaticFC.mat                          % Load pearson sFC
% Finding Phi and PhiIntra
for i=1:343                                             % Loop over subjects
        Subject=char(Subjects(2*i-1));                  % Find subject day 1 FC                            
        eval(['ints = ' Subject ';']);                  % Save day 1 fc to a variable
        su(:,i)=ints;                                   % Save fc to variable su
        for j=1:343                                     % Loop over subjects to find day 2 session
            Subject=char(Subjects(2*j));                % Load subject day 2 FC                           
            eval(['ints = ' Subject ';']);              % Save fc to a variable
            su(:,j)=ints;                               % Saving value to su
            phi(:,i,j) = su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));        % Finding phi value
            if i==j                                                             % Finding whether subject is same
                phiintra(:,i)=su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));   % Saving phi value for same subject
            end                         
        end 
end
GroupPhi=mean(phiintra,2);                              % Finding mean phi for subjects

% Probability calculation
% Loop below finds probability of phi of a subject greter than between
% subjects
for i=1:343                                             % Looping over subjects 
    b=1:343;                                            % Looping over subjects 
    c=setdiff(b,i);                                     % Remove subject i
    edge(:,i)=zeros(4005,1);                            % Initializing edge variable
    for j=1:size(c,2)                                   % Looping over subjects removing subject i
        k=phiintra(:,i)>squeeze(phi(:,i,c(j)));         % Finding phi(i,i)>phi(i,j)
        l=phiintra(:,i)>squeeze(phi(:,c(j),i));         % Finding phi(i,i)>phi(j,i)
        edge(:,i)=edge(:,i)+k+l;                        % Finding number of times edge was having higher phi
    end
    p(:,i)=edge(:,i)./(2*342);                          % Finding probabilty  
end
avgp=mean(p,2);                                         % Finding differential power
I=find(triu(ones(90),1));                               % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=avgp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding the diferential power in matrix format

% Plot differential power values for edges

figure;imagesc(df1);c=jet;colormap(c);axis image;colorbar;                      % Ploting differential power in new window
title('Differential power Partial sFC')                                         % Adding title to image
ax = gca;ax.XTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5]; % Adding tick marks for label
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
set(gca,'XTickLabelRotation',90,'fontsize',13)                                  % Seting fontsize and rotation
ax.YTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5];          % Adding tickmarks for label
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
l1=[3.5 3.5;8.5 8.5 ; 17.5  17.5 ;19.5  19.5 ;26.5  26.5 ;32.5  32.5 ;38.5  38.5 ;50.5  50.5 ;54.5  54.5 ;56.5  56.5 ;62.5  62.5 ;69.5  69.5 ;79.5  79.5]; % Creating line centroids
l2=[0.5 90.5;0.5 90.5; 0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5 ];                               % Creating line centroids

hold on                                                 % Hold figure
for i= 1:13                                             % Loop to plot lines
    line(l1(i,:),l2(i,:),'color','w')                   % Plotting line             
    line(l2(i,:),l1(i,:),'color','w')                   % Plotting lines
end 

% Probability Percentile 

temp=avgp;                                              % Adding differential power to a temporary variable
y=prctile(avgp,95);                                     % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                         % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
title('Edge contribution Partial sFC')                  % Add title to the figure
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

%% Pearson dFC differential power
% This part of program finds differential power for Pearson dFC
clear                                                   % Clear workspace
load Pearson_Rest_DynamicFC.mat                         % Load pearson sFC
% Finding Phi and PhiIntra
for i=1:343                                             % Loop over subjects
        Subject=char(Subjects(2*i-1));                  % Find subject day 1 FC                            
        eval(['ints = ' Subject ';']);                  % Save day 1 fc to a variable
        su(:,i)=ints;                                   % Save fc to variable su
        for j=1:343                                     % Loop over subjects to find day 2 session
            Subject=char(Subjects(2*j));                % Load subject day 2 FC                           
            eval(['ints = ' Subject ';']);              % Save fc to a variable
            su(:,j)=ints;                               % Saving value to su
            phi(:,i,j) = su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));        % Finding phi value
            if i==j                                                             % Finding whether subject is same
                phiintra(:,i)=su(:,i).*su(:,j)/(norm(su(:,i))*norm(su(:,j)));   % Saving phi value for same subject
            end                         
        end 
end
GroupPhi=mean(phiintra,2);                              % Finding mean phi for subjects

% Probability calculation
% Loop below finds probability of phi of a subject greter than between
% subjects
for i=1:343                                             % Looping over subjects 
    b=1:343;                                            % Looping over subjects 
    c=setdiff(b,i);                                     % Remove subject i
    edge(:,i)=zeros(4005*4,1);                            % Initializing edge variable
    for j=1:size(c,2)                                   % Looping over subjects removing subject i
        k=phiintra(:,i)>squeeze(phi(:,i,c(j)));         % Finding phi(i,i)>phi(i,j)
        l=phiintra(:,i)>squeeze(phi(:,c(j),i));         % Finding phi(i,i)>phi(j,i)
        edge(:,i)=edge(:,i)+k+l;                        % Finding number of times edge was having higher phi
    end
    p(:,i)=edge(:,i)./(2*342);                          % Finding probabilty  
end
davgp=mean(p,2);                                        % Finding differential power for edges
avgp=davgp(1:4005)+davgp(4006:8010)+davgp(8011:12015)+davgp(12016:16020);% Find average differential power fo edges
I=find(triu(ones(90),1));                               % Finding upper triangular positions
dfcx1=zeros(90,90);dfcx1(I)=avgp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding the diferential power in matrix format

% Plot differential power values for edges

figure;imagesc(df1);c=jet;colormap(c);axis image;colorbar;                      % Ploting differential power in new window
title('Differential power Pearson dFC')                                         % Adding title to image
ax = gca;ax.XTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5]; % Adding tick marks for label
ax.XTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
set(gca,'XTickLabelRotation',90,'fontsize',13)                                  % Seting fontsize and rotation
ax.YTick =[2 6 13 18.5 23 29.5 35.5 44.5 52.5 55.5 59.5 66 74.5 85.5];          % Adding tickmarks for label
ax.YTickLabel ={'AUD','BAS','D-DMN','V2','LAN','LECN','SMOTOR','PSAL','PRECUNEUS','V1','RECN','ASAL','VDMN','VISUO'};% Adding labels for tick mark
l1=[3.5 3.5;8.5 8.5 ; 17.5  17.5 ;19.5  19.5 ;26.5  26.5 ;32.5  32.5 ;38.5  38.5 ;50.5  50.5 ;54.5  54.5 ;56.5  56.5 ;62.5  62.5 ;69.5  69.5 ;79.5  79.5]; % Creating line centroids
l2=[0.5 90.5;0.5 90.5; 0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5;0.5 90.5 ];                               % Creating line centroids

hold on                                                 % Hold figure
for i= 1:13                                             % Loop to plot lines
    line(l1(i,:),l2(i,:),'color','w')                   % Plotting line             
    line(l2(i,:),l1(i,:),'color','w')                   % Plotting lines
end 

% Probability Percentile 

temp=avgp;                                              % Adding differential power to a temporary variable
y=prctile(avgp,95);                                     % Finding 95 percentile value
temp(temp<y)=0;temp(temp>=y)=1;                         % Thresholding at 95 percentile
dfcx1=zeros(90,90);dfcx1(I)=temp;df1=triu(dfcx1)+triu(dfcx1,1)';% Finding differential power edges in matrix format 
network=[];                                     % Initializing network variable
nnzb=nnz(temp)*2;                               % Finding number of active edges
network(1,1)=sum(sum(df1(1:3,1:3)))/(nnzb);     % Finding percentage edges in network 1
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
title('Edge contribution Pearson dFC')                  % Add title to the figure
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
