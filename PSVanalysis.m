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

autoAvg = 1;                    % 0= use avg times below
                                % 1= detect light on from data file

get_mag = 0;                    % use to extract magnetic component of data
                                % Files must be in order:
                                %   x1_+h, x1_-h, x2_+h, x2_-h...
include_nonmag = 0;             % include non-magnetic component in output

avgStart = 55;                  % when to start data average (s)
avgStop = 95;                   % when to stop data average (s)
cutoff = 75;                    % when to truncate data file (s)

files = dir('*.txt');
numFiles = size(files,1);

if numFiles == 0
    error('PSVanalysis:noData',['*** Error: No fitting data found ***'...
          '\n\tData must be *.txt file.']);
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

%% Read files

for i=1:numFiles
    % import data file
    data = importdata(files(i).name,'\t',1);
    
    % extract part of filename after dash (and remove file extension)
    [~, name] = strtok(files(i).name,'-');
    dataName = name(2:end-4);
    
    % setup output headers
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
    
    %% Find and subtract background
    
    [~, zeroIndex] = min(abs(data.data(:,1)-zeroTime));
    zeroAvg = mean(data.data(1:zeroIndex,2));
    zeroAvgT = mean(data.data(1:zeroIndex,3));
    
    fixData(:,1) = data.data(1:cutoffIndex,1);             % time
    fixData(:,2) = data.data(1:cutoffIndex,2)-zeroAvg;     % voltage
    fixData(:,3) = data.data(1:cutoffIndex,3)-zeroAvgT;    % delta T
    fixData(:,4) = data.data(1:cutoffIndex,5);             % light on
    
    %% Prep for voltage averaging
    % period to be averaged is made imaginary if using autoAvg
    
    if autoAvg == 1
        for k = 3:cutoffIndex-2
            % if the light is on for the previous two and next two data
            % points then flag for averaging
            if fixData(k-2,4)==100 && fixData(k+2,4)==100
                fixData(k,2) = fixData(k,2)*1j;
            end
        end
    end
    
    %% Write output
    
    for k=1:length(fixData)
        if tempData == 1
            DataOutput{k+3,3*i-2} = fixData(k,1);
            DataOutput{k+3,3*i-1} = fixData(k,2);
            DataOutput{k+3,3*i} = fixData(k,3);
        else
            DataOutput{k+3,2*i-1} = fixData(k,1);
            DataOutput{k+3,2*i} = fixData(k,2);
        end
    end
    
end

%% Extract magnetic component
% must set get_mag = 1 in setup section
% data must be in order: x1_+h, x1_-h, x2_+h, x2_-h...

if (get_mag==1) && (tempData~=1)
    for i=1:size(DataOutput,2)/4
        % find magnetic and non-magnetic components
        for k = 4:length(DataOutput)
            pos_h_t = DataOutput{k,4*i-3};
            pos_h_v = DataOutput{k,4*i-2};
            neg_h_t = DataOutput{k,4*i-1};
            neg_h_v = DataOutput{k,4*i  };
            
            avgT = (pos_h_t + neg_h_t)/2;
            mag = (pos_h_v - neg_h_v)/2;
            if include_nonmag
                nonMag = (pos_h_v + neg_h_v)/2;
            end
            
            if include_nonmag
                DataOutput{k,4*i-3} = avgT;
                DataOutput{k,4*i-2} = mag;
                DataOutput{k,4*i-1} = avgT;
                DataOutput{k,4*i  } = nonMag;
            else
                DataOutput{k,4*i-3} = avgT;
                DataOutput{k,4*i-2} = mag;
                DataOutput{k,4*i-1} = '';
                DataOutput{k,4*i  } = '';
            end
        end
        % Change output labels
        if include_nonmag
            tmpStr = DataOutput{3,4*i-2};
            DataOutput{3,4*i-2} = strcat(tmpStr(1:end-2),'Mag');
            tmpStr = DataOutput{3,4*i  };
            DataOutput{3,4*i  } = strcat(tmpStr(1:end-2),'Non-Mag');
        else
            tmpStr = DataOutput{3,4*i-2};
            DataOutput{3,4*i-2} = strcat(tmpStr(1:end-2),'Mag');
            DataOutput{k,4*i  } = '';
        end
            
    if ~include_nonmag
        DataOutput{1,4*i-1} = '';
        DataOutput{2,4*i-1} = '';
        DataOutput{3,4*i-1} = '';
        DataOutput{1,4*i  } = '';
        DataOutput{2,4*i  } = '';
        DataOutput{3,4*i  } = '';
    end
    
    end
    % Delete empty columns
    DataOutput(:,find(all(cellfun(@isempty,DataOutput),1)))=[];
end

%% Find average voltage signal
% if using autoAvg: average imaginary numbers, make them real again
% otherwise average over times specified in setup

for i=1:size(DataOutput,2)/2
    if autoAvg == 1
        avgData = zeros(length(DataOutput));
        for k = 4:length(DataOutput)
            if (isreal(DataOutput{k-2,2*i})==0)
                avgData(k) = imag(DataOutput{k,2*i});
                DataOutput{k-2,2*i} = imag(DataOutput{k-2,2*i});
            end
        end
        avgData = avgData(avgData~=0);
        
        avgV = mean(avgData);
        stdV = std(avgData);
    else
        avgData = zeros(length(DataOutput));
        avgTime = zeros(length(DataOutput));
        for k = 4:length(DataOutput)
            avgData(k) = DataOutput{k,i*2};
            avgTime(k) = DataOutput{k,i*2-1};
        end
        
        [~, startIndex] = min(abs(avgTime-avgStart));
        [~, stopIndex] = min(abs(avgTime-avgStop));
        
        avgV = mean(avgData(startIndex:stopIndex));
        stdV = std(avgData(startIndex:stopIndex));
    end
    
    VoltageOutput{i+3,1} = DataOutput{3,2*i};
    VoltageOutput{i+3,2} = avgV;
    VoltageOutput{i+3,3} = stdV;
    
end

%% Find percent of max voltage

voltages = abs(cell2mat(VoltageOutput(3:end,2)));
volt_err = abs(cell2mat(VoltageOutput(3:end,3)));
[maxV, max_err] = max(voltages);
percent_volt = round((voltages./maxV)*10000)/100;
for l = 1:length(voltages)
    VoltageOutput{l+3,4} = voltages(l)/maxV *100 ;
    VoltageOutput{l+3,5} = voltages(l)/maxV *100 * ...
        sqrt((volt_err(max_err)/maxV)^2 + (volt_err(l)/voltages(l))^2);
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
    fprintf(2,['\n*** Error: DataOutput.dat already open.'...
               ' Close file and press the any key *** \n']);
    pause;
    fileID1 = fopen('DataOutput.dat','w');
    if fileID1 == -1
        error(['*** Error: DataOutput.dat is still open.'...
               ' Run the program again ***']);
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
    fprintf(2,['\n*** Error: VoltageOutput.dat already open.'...
               ' Close file and press the any key *** \n']);
    pause;
    fileID2 = fopen('VoltageOutput.dat','w');
    if fileID2 == -1
        error(['*** Error: VoltageOutput.dat is still open.'...
               ' Run the program again *** ']);
    end
end
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{1,:});
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{2,:});
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',VoltageOutput{3,:});
for row = 4:length(VoltageOutput)
    fprintf(fileID2,'%s\t%6.4f\t%6.4f\t%4.2f\t%4.2f\n',...
            VoltageOutput{row,:});
end
fclose(fileID2);

fclose('all');
fprintf('\nDone\n');
