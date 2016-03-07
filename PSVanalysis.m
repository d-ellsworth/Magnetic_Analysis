%% PSV Data Analysis
% David Ellswoth
% 2015-11-09
%
% Reads data from all .txt files in folder
% Subtracts background voltage and finds average signal over time interval
% Output corrected data to DataOutput.dat
% Output averaged voltages to VoltageOutput.dat

clc;
clear all;
fprintf('\nRunning\n');
%% Setup

tempData = 0;                   % use temperature data?

zeroTime = 45;                  % time to zero average (s)

autoAvg = 1;                    % 0= use avg times below, 1= detect light on from data file

avgStart = 55;                  % when to start data average (s)
avgStop = 95;                   % when to stop data average (s)
cutoff = 300;                   % when to truncate data file (s)

files = dir('*.txt');
numFiles = size(files,1);

if numFiles == 0
    error('PSVanalysis:noData','*** Error: No fitting data found ***\n\tData must be *.txt file.');
end

% initialize outputs
VoltageOutput = cell(numFiles+3,5);
VoltageOutput{1,1} = 'Filter';
VoltageOutput{1,2} = 'Average Voltage';
VoltageOutput{2,2} = 'uV';
VoltageOutput{1,3} = 'StDev';
VoltageOutput{1,4} = 'Percent';
VoltageOutput{1,5} = 'P Error';


testData = importdata(files(1).name,'\t',1);
if autoAvg == 1
    cutoffIndex = find(testData.data(:,5)==999)-1;
else
    [~, cutoffIndex] = min(abs(testData.data(:,1)-cutoff));
end
testData = testData.data(1:cutoffIndex);
DataOutput = cell(length(testData)+3,numFiles*2);

%% Read and modify files

for i=1:numFiles
    %% Import data file
    data = importdata(files(i).name,'\t',1);
    
    % extract part of filename after dash (and remove file extension)
    [~, name] = strtok(files(i).name,'-');
    dataName = name(2:end-4);
    
    %% Find and subtract background
    [~, zeroIndex] = min(abs(data.data(:,1)-zeroTime));
    zeroAvg = mean(data.data(1:zeroIndex,2));
    
    normData(:,1) = data.data(1:cutoffIndex,1);
    normData(:,2) = (data.data(1:cutoffIndex,2)-zeroAvg);
    zeroAvgT = mean(data.data(1:zeroIndex,3));
    normData(:,3) = (data.data(1:cutoffIndex,3)-zeroAvgT);
    normData(:,4) = data.data(1:cutoffIndex,5);
    
    if tempData == 1
        DataOutput{1,3*i-2} = 'Time';
        DataOutput{2,3*i-2} = 's';
        DataOutput{1,3*i-1} = 'Voltage';
        DataOutput{2,3*i-1} = 'uV';
        DataOutput{3,3*i-1} = dataName;
        DataOutput{1,3*i} = 'Delta T';
        DataOutput{2,3*i} = 'C';
    else
        DataOutput{1,2*i-1} = 'Time';
        DataOutput{2,2*i-1} = 's';
        DataOutput{1,2*i} = 'Voltage';
        DataOutput{2,2*i} = 'uV';
        DataOutput{3,2*i} = dataName;
    end
    
    for j=1:length(normData)
        if tempData == 1
            DataOutput{j+3,3*i-2} = normData(j,1);
            DataOutput{j+3,3*i-1} = normData(j,2);
            DataOutput{j+3,3*i} = normData(j,3);
        else
            DataOutput{j+3,2*i-1} = normData(j,1);
            DataOutput{j+3,2*i} = normData(j,2);
        end
    end
    
    %% Find avgerage voltage signal
    if autoAvg == 1
        lightOn = zeros(sum(normData(:,4)~=0));
        for j = 2:cutoffIndex-1
            if normData(j-1,4)==100 && normData(j+1,4)==100
                lightOn(j) = j;
            end
        end
        lightOn = lightOn(lightOn~=0);
        
        avgV = mean(normData(lightOn,2));
        stdV = std(normData(lightOn,2));
    else
        [~, startIndex] = min(abs(normData(:,1)-avgStart));
        [~, stopIndex] = min(abs(normData(:,1)-avgStop));
        
        avgV = mean(normData(startIndex:stopIndex,2));
        stdV = std(normData(startIndex:stopIndex,2));
    end
    
    VoltageOutput{i+3,1} = dataName;
    VoltageOutput{i+3,2} = avgV;
    VoltageOutput{i+3,3} = stdV;
    
end

%% Find percent of max voltage

voltages = abs(cell2mat(VoltageOutput(3:end,2)));
volt_err = abs(cell2mat(VoltageOutput(3:end,3)));
[maxV, max_err] = max(voltages);
for l = 1:length(voltages)
    VoltageOutput{l+3,4} = voltages(l)/maxV *100 ;
    VoltageOutput{l+3,5} = voltages(l)/maxV *100 * sqrt( (volt_err(max_err) / maxV)^2 + (volt_err(l)/voltages(l))^2);
end

%% Output results

[nrows, ncols] = size(DataOutput);
formatSpec1 = char('');
formatSpec2 = char('');
for k = 1:ncols
    formatSpec1 = strcat(formatSpec1,'%s\t');
    formatSpec2 = strcat(formatSpec2,'%6.4g\t');
end
formatSpec1 = strcat(formatSpec1,'\n');
formatSpec2 = strcat(formatSpec2,'\n');

fileID1 = fopen('DataOutput.dat','w');
if fileID1 == -1
    fprintf(2,'\n*** Error: DataOutput.dat already open. Close file and press the any key *** \n');
    pause;
    fileID1 = fopen('DataOutput.dat','w');
    if fileID1 == -1
        error('*** Error: DataOutput.dat is still open. Run the program again ***');
    end
end
fprintf(fileID1,formatSpec1,DataOutput{1,:});
fprintf(fileID1,formatSpec1,DataOutput{2,:});
fprintf(fileID1,formatSpec1,DataOutput{3,:});
for row = 4:nrows
    fprintf(fileID1,formatSpec2,DataOutput{row,:});
end
fclose(fileID1);


fileID2 = fopen('VoltageOutput.dat','w');
if fileID2 == -1
    fprintf(2,'\n*** Error: VoltageOutput.dat already open. Close file and press the any key *** \n');
    pause;
    fileID2 = fopen('VoltageOutput.dat','w');
    if fileID2 == -1
        error('*** Error: VoltageOutput.dat is still open. Run the program again *** ');
    end
end
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{1,:});
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{2,:});
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{3,:});
for row = 4:length(VoltageOutput)
    fprintf(fileID2,'%s\t%6.4f\t%6.4f\t%4.2f\t%4.2f\n',VoltageOutput{row,:});
end
fclose(fileID2);


fclose('all');
fprintf('\nDone\n');
