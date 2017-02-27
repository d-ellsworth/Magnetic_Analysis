%% PSV / SSE Analysis
% David Ellsworth
% 2016-11-16

% Do all the analysis to separate PSV and SSE signals.
%
% Dependencies: PSVanaFun.m, AvgFun.m
%
% Takes input .txt data files, averages them, seperates PSV & SSE
% components and outputs analysis.

clc;
clf;
clear all;
fprintf('Running\n');

%% Setup

% avg function defaults
filesToAvgIn = 5;

% analysis function defaults
timeZero = 45;
avgAuto = 1;
magGet = 1;
startAvg = 80;
stopAvg = 95;
stopData = 150;

% what lights are used?
lights = {'LED' 'Hal' '950L' 'IR'};
integratedLight = [21769, 123728, 75654, 2093];
lightScale = 10^5;
normLight = 'LED';
psvGuess = 0.5;
sseGuess = 0.1;

header = {'Time (s)','Voltage (uV)','Correction Factor','null','null'};
lightData = cell(3+length(lights)*2,5);
lightOut = cell(length(lightData),6);

% find where we are and make a path
myPath = mfilename('fullpath');
myPath = myPath(1:end-length(mfilename));
addpath(genpath(myPath));

cd(myPath);

%% Average files

AvgFun(filesToAvgIn)

%% Get magnetic components

cd(strcat(myPath,'average\'));
PSVanaFun(timeZero,avgAuto,magGet,startAvg,stopAvg,stopData);

%% Separate PSV & SSE for each light
rawData = importdata('DataOutput.dat');
for i = 1:length(lights)
    fprintf(['\nSeparating ' ,lights{i},' data\n']);
    % search for appropriate data
    for j = 1:length(rawData.colheaders)
        % find correct data for light
        if ~cellfun('isempty',strfind(rawData.colheaders(j),char(lights(i))))
            if ~cellfun('isempty',strfind(rawData.colheaders(j),'up'))
                %find "up" data
                upData = rawData.data(:,j-1:j);
            elseif ~cellfun('isempty',strfind(rawData.colheaders(j),'down'))
                % find "down" data
                downData = rawData.data(:,j-1:j);
            end
            
        end
    end
    
    % Get intial PSV & SSE components
    timePSV = (upData(:,1) + downData(:,1))/2;
    rawPSV = (upData(:,2) + downData(:,2))/2;
    fixPSV = rawPSV;
    psvFactor = psvGuess;
    
    timeSSE = (upData(:,1) + downData(:,1))/2;
    rawSSE = (upData(:,2) - downData(:,2))/2;
    fixSSE = rawSSE;
    sseFactor = sseGuess;
    
    % make PSV plot
    f1 = figure;
    p1 = plot(timePSV,rawPSV);
    hold on;
    p2 = plot(timePSV,fixPSV);
    set(p1, 'LineWidth', 3.0, 'Color', 'blue');
    xlabel('Time (s)');
    ylabel('Voltage (uV)');
    
    % tune PSV correction factor
    while 1
        fixPSV = rawPSV + (rawSSE * psvFactor);
        
        delete(p2);
        p2 = plot(timePSV,fixPSV);
        set(p2, 'LineWidth', 3.0, 'Color', 'red');
        title(strcat(lights(i),' PSV Current factor: ',num2str(psvFactor)),'FontSize',16);
        legend('Original','Corrected');
        
        newFactor = input('Input PSV correction factor, or "done": ','s');
        if strcmp(newFactor,'done') ||  strcmp(newFactor,'')
            break
        elseif ~isempty(str2double(newFactor))
            psvFactor = str2double(newFactor);
        else
            continue
        end
        figure(f1);
    end
    
    % make SSE plot
    f2 = figure;
    s1 = plot(timeSSE,rawSSE);
    hold on;
    s2 = plot(timeSSE,fixSSE);
    set(s1, 'LineWidth', 3.0, 'Color', 'blue');
    xlabel('Time (s)');
    ylabel('Voltage (uV)');
    
    % tune SSE correction factor
    while 1
        fixSSE = rawSSE + (rawPSV * sseFactor);
        
        delete(s2);
        s2 = plot(timeSSE,fixSSE);
        set(s2, 'LineWidth', 3.0, 'Color', 'red');
        title(strcat(lights(i),' SSE Current factor: ',num2str(sseFactor)),'FontSize',16);
        legend('Original','Corrected');
        
        newFactor = input('Input SSE correction factor, or "done": ','s');
        if strcmp(newFactor,'done') ||  strcmp(newFactor,'')
            break
        elseif (~isempty(str2double(newFactor)))
            sseFactor = str2double(newFactor);
        else
            continue
        end
        figure(f2);
    end
    
    p = get(gcf,'Position');
    set(0,'DefaultFigurePosition',p);
    hold off;
    close all;
    
    %% Output corrected data
    
    fprintf(['Writing ' ,lights{i},' data\n']);
    
    corrDir = strcat(myPath,'corrected\');
    if ~exist(corrDir,'file')
        mkdir(corrDir);
    end
    
    fileNamePSV = strcat(corrDir,strcat('Corrected-',lights{i},'_PSV.txt'));
    fileNameSSE = strcat(corrDir,strcat('Corrected-',lights{i},'_SSE.txt'));
    
    % open file for PSV output
    fileID1 = fopen(fileNamePSV,'w');
    if fileID1 == -1
        fprintf(2,['\n*** Error: Output file is already open.'...
            ' Close file and press the any key *** \n']);
        pause;
        fileID1 = fopen(fileNamePSV,'w');
        if fileID1 == -1
            error(['*** Error: Output file is still open.'...
                ' Run the program again ***']);
        end
    end
    % write PSV output
    fprintf(fileID1,'%s\t%s\t%s\t%s\t%s\n',header{:});
    for row = 1:length(fixPSV)
        fprintf(fileID1,'%6.4g\t%6.4g\t%6.4g\t%6.4g\t%6.4g\n',[timePSV(row),fixPSV(row),psvFactor,0,0]);
    end
    fclose(fileID1);
    
    lightOut{i+3,6} = psvFactor;
    
    % open file for SSE output
    fileID2 = fopen(fileNameSSE,'w');
    if fileID2 == -1
        fprintf(2,['\n*** Error: Output file is already open.'...
            ' Close file and press the any key *** \n']);
        pause;
        fileID2 = fopen(fileNameSSE,'w');
        if fileID2 == -1
            error(['*** Error: Output file is still open.'...
                ' Run the program again ***']);
        end
    end
    % write SSE output
    fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\n',header{:});
    for row = 1:length(fixPSV)
        fprintf(fileID2,'%6.4g\t%6.4g\t%6.4g\t%6.4g\t%6.4g\n',[timeSSE(row),fixSSE(row),sseFactor,0,0]);
    end
    fclose(fileID2);
    
    lightOut{length(lights)+i+3,6} = sseFactor;
    
end


%% Get PSV & SSE component values

fprintf('\nExtracting PSV & SSE values...\n');
cd(corrDir);
PSVanaFun(timeZero,0,0,startAvg,stopAvg,stopData);

% Read data from VoltageOutput.dat
fileID0 = fopen('VoltageOutput.dat');
lightData1 = textscan(fileID0,'%s %s %s %s %s',3,'Delimiter','\t');
lightData2 = textscan(fileID0,'%s %f %f %f %f','Delimiter','\t');

for col=1:5
    for row1=1:length(lightData1{1,1})
        lightData{row1,col}=lightData1{1,col}{row1};
    end
    for row2=1:length(lightData2{1,1})
        lightData{row2+length(lightData1{1,1}),col}=lightData2{1,col}(row2);
    end
end
fclose(fileID0);

%% Divide voltage by area under light curves

lightOut{1,1} = 'Filter';
lightOut{1,2} = 'Voltage';
lightOut{2,2} = 'uV';
lightOut{1,3} = 'Voltage / Area';
lightOut{2,3} = 'a.u.';
lightOut{1,4} = 'Normalized Voltage / Area';
lightOut{2,4} = 'a.u.';
lightOut{1,5} = 'Integrated Light';
lightOut{2,5} = 'a.u.';
lightOut{1,6} = 'Correction Factor';

for i=1:length(lights)
    % find correct data for light
    for j=1:(length(lightData)-3)
        if ~cellfun('isempty',strfind(lightData{j+3},char(lights(i))))
            if ~cellfun('isempty',strfind(lightData{j+3},'PSV'))
                %find "PSV" data
                lightOut{i+3,1} = lightData{j+3}{1};
                lightOut{i+3,2} = lightData{j+3,2};
                lightOut{i+3,3} = lightData{j+3,2} / integratedLight(i) * lightScale;
                lightOut{i+3,5} = integratedLight(i);
                %lightOut{i+3,6} = lightData{j+3,6};
                if ~isempty(strfind(normLight,char(lights(i))))
                    normPSV = lightData{j+3,2} / integratedLight(i) * lightScale;
                end
            elseif ~cellfun('isempty',strfind(lightData{j+3},'SSE'))
                % find "SSE" data
                lightOut{(length(lights)+i)+3,1} = lightData{j+3}{1};
                lightOut{(length(lights)+i)+3,2} = lightData{j+3,2};
                lightOut{(length(lights)+i)+3,3} = lightData{j+3,2} / integratedLight(i) * lightScale;
                lightOut{(length(lights)+i)+3,5} = integratedLight(i);
                %lightOut{(length(lights)+i)+3,6} = lightData{j+3,6};
                if ~isempty(strfind(normLight,char(lights(i))))
                    normSSE = lightData{j+3,2} / integratedLight(i) * lightScale;
                end
            end
        end
    end
end

% Normalize to normLight
for j=1:length(lightOut)-3
    if ~isempty(strfind(lightOut{j+3},'PSV'))
        lightOut{j+3,4} = lightOut{j+3,3} / normPSV;
    elseif ~isempty(strfind(lightOut{j+3},'SSE'))
        lightOut{j+3,4} = lightOut{j+3,3} / normSSE;
    end
end

% open file for normalized output
cd(myPath);

fileID3 = fopen('LightData.dat','w');
if fileID3 == -1
    fprintf(2,['\n*** Error: Output file is already open.'...
        ' Close file and press the any key *** \n']);
    pause;
    fileID3 = fopen('LightData.dat','w');
    if fileID3 == -1
        error(['*** Error: Output file is still open.'...
            ' Run the program again ***']);
    end
end
% write normalized output
fprintf(fileID3,'%s\t%s\t%s\t%s\t%s\t%s\n',lightOut{1,:});
fprintf(fileID3,'%s\t%s\t%s\t%s\t%s\t%s\n',lightOut{2,:});
fprintf(fileID3,'%s\t%s\t%s\t%s\t%s\t%s\n',lightOut{3,:});
for row = 4:length(lightOut)
    fprintf(fileID3,'%s\t%6.6g\t%6.6g\t%6.6g\t%6.6g\t%3.3g\n',lightOut{row,:});
end
fclose(fileID3);

%% The End

fclose('all');
fprintf('\nJob''s done.\n');