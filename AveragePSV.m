%% Average Data Files
% David Ellsworth
% 2016-10-18
%
% Averages multiple .txt data files into one output.
% 


clc;
clear all;
fprintf('\nRunning\n');
%% Setup

% number of files to average together
% total file number must be a multiple of this
filesToAvg = 5;


files = dir('*.txt');
numFiles = size(files,1);

if numFiles == 0
    error('AveragePSV:noData',['*** Error: No data files found ***' ...
          '\n\tData must be *.txt file.']);
end

% initialize output
testData = importdata(files(1).name,'\t',1);
numHeader = size(testData.colheaders,1);        % number of header rows
numCols = size(testData.data,2);                % number of columns

AvgOutput = cell(length(testData.data)+numHeader,numCols);
AvgOutput(:) = {0};

%% Averaging

for i=1:numFiles/filesToAvg
    
    % reinitialize stuff
    clearvars data
    AvgOutput(:) = {0};
    avgIndex = 1;
    
    % Average in groups of filesToAvg
    fileStart = 1+filesToAvg*(i-1);
    
    for m=fileStart:fileStart+filesToAvg-1
        % Import data file
        data = importdata(files(m).name,'\t',1);
        %disp(['File: ',files(m).name]);
        
        % Average
        for k=1:length(AvgOutput)-numHeader-1
            for n=1:numCols
                AvgOutput{k+numHeader,n} = AvgOutput{k+numHeader,n} + ...
                    (data.data(k,n)-AvgOutput{k+numHeader,n})/avgIndex;
            end
        end
        AvgOutput{length(AvgOutput),numCols} = 999;
        avgIndex = avgIndex+1;
    end
    
    
    % import header rows
    for l=1:numCols
        AvgOutput{1,l} = testData.colheaders{l};
    end
    
    %% Output
    
    % how to format the output
    [nrows, ncols] = size(AvgOutput);
    formatSpec1 = char('');
    formatSpec2 = char('');
    for k = 1:ncols
        formatSpec1 = strcat(formatSpec1,'%s\t');
        formatSpec2 = strcat(formatSpec2,'%6.7g\t');
    end
    formatSpec1 = strcat(formatSpec1,'\n');
    formatSpec2 = strcat(formatSpec2,'\n');
    
    % make a directory to stick files if needed
    directory = strcat(pwd,'\average\');
    if ~exist(directory,'file')
        mkdir(directory);
    end
    
    [~, name] = strtok(files(i*filesToAvg).name,'-');
    fileName = fullfile(directory,strcat('avg',num2str(i,'%03d'),name));
    
    % open file
    fileID1 = fopen(fileName,'w');
    if fileID1 == -1
        fprintf(2,['\n*** Error: Output file is already open.'...
                   ' Close file and press the any key *** \n']);
        pause;
        fileID1 = fopen(fileName,'w');
        if fileID1 == -1
            error(['*** Error: Output file is still open.'...
                   ' Run the program again ***']);
        end
    end
    % write the output
    for row = 1:numHeader
        fprintf(fileID1,formatSpec1,AvgOutput{row,:});
    end
    for row = numHeader+1:nrows
        fprintf(fileID1,formatSpec2,AvgOutput{row,:});
    end
    fclose(fileID1);
    
end

%% The End

fclose('all');
fprintf('\nDone\n');