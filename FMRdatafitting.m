%% FMR data fitting
% Praveen FMR data fitting 13/Nov/2013
% David Ellsworth update 2014-05-23
%
% Put this .m file in the folder with the .txt files which contain FMR data
%
% Will fit each data set with Lorentz fit (blue trace is final fit)
% Indivitual fitting results output to FittingResults.dat
%
% Also outputs Kittel and Gilbert damping fits for total data set
% to AnalysisResults.dat

clc;
clear all;
fprintf('\nRunning\n');

%% Modify These Things!

% Do you want to plots?
plots = 2;              % 0=no plots, 1=only kittel&damping, 2=all plots
savePlots = 0;          % 0=don't save plots, 1=save plots to .jpg files

% Geometry for Kittel/Gilbert fitting
Geometry = 0;           % 0=in plane, 1=out of plane, 2=no kittel/gilbert fitting (use for polar angle analysis)

% Select data range
SetDataRange = 2;       % 0=use % of data range, 1=use fixed data range, 2=use full data range
                        % set to 0 or 2 if data files have different number of points

DataRange = 150;        % Data points to either side to be considerd (for FixDataRange = 1)
DataFraction = 1/3;     % part of data range to either side to be considered (for FixDataRange = 0)

% Initialize or fix fitting parameters
gammaGuess = 2.8e6;     % initial gamma value for kittel fitting
fixGamma = 0;           % fix gamma value in fitting, use guess value
MsGuess = 1750;         % initial 4piMs value for kittel fitting
fixMs = 0;              % fix 4piMs value for fitting, use guess value

%% Other parameter setup

files=dir('*.txt');
numFiles=size(files,1);

if numFiles == 0
    error('FMRdatafitting:noData','*** Error: No fitting data found ***\n\tData must be *.txt file.');
end

numplot = 0;
FMRfreq = zeros(1,numFiles);
FMRlineWidth = zeros(1,numFiles);
Linewidtherror = zeros(1,numFiles);
FMRField = zeros(1,numFiles);
Fielderror = zeros(1,numFiles);

FittingOutput = cell(numFiles+3,6);
FittingOutput{1,1} = 'Filename';
FittingOutput{1,2} = 'Frequency';
FittingOutput{2,2} = 'GHz';
FittingOutput{1,3} = 'Linewidth';
FittingOutput{2,3} = 'Oe';
FittingOutput{1,4} = 'LineWidthError';
FittingOutput{1,5} = 'FMR Field';
FittingOutput{2,5} = 'Oe';
FittingOutput{1,6} = 'FieldError';


%% read all the files & Lorentz fit em
for i = 1:numFiles
    clc;
    fprintf('\nCurrent Fit - %s\n',files(i).name);
    if plots>1
        numplot = i;i;
    end
    data = importdata(files(i).name, '\t', 1);
        
    if (SetDataRange==0)
        DataRange=round((length(data.data))*DataFraction);
    end
    
    xdataFMR = data.data(1:length(data.data),1);
    ydataFMR = data.data(1:length(data.data),2);
    %% the function to be fit
    % y0+b*x-(16*a/PI)*sqrt(3)*(xr-xl)*(x-(xl+xr)/2)/(4*(x-(xl+xr)/2)^2+3*(xr-xl)^2)^2
    
    %% least squares fit
    [maxval,maxIndx] = max(data.data(1:length(data.data),2));
    [minval,minIndx] = min(data.data(1:length(data.data),2));
    xr = data.data(maxIndx,1); xl = data.data(minIndx,1);
    
    if SetDataRange < 2
        xdataFMR = xdataFMR(round(0.5*(maxIndx+minIndx))-DataRange:round(0.5*(maxIndx+minIndx))+DataRange);
        ydataFMR = ydataFMR(round(0.5*(maxIndx+minIndx))-DataRange:round(0.5*(maxIndx+minIndx))+DataRange);
    end
    
    F = @(x,xdata)x(1)+x(2)*xdata-((16/pi)*x(3))*sqrt(3)*(x(4)-x(5))*(xdata-0.5*(x(4)+x(5)))./((4*(xdata-0.5*(x(4)+x(5))).^2+3*(x(4)-x(5)).^2).^2);
    x0 = [0 0 1 xr xl];
    [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,xdataFMR,ydataFMR);
    
    if plots>1
        figure(i)
        clf;
        plot(xdataFMR,ydataFMR,'r.')
        title(files(i).name,'Interpreter','none')
        % y0 = 0; b =0; a = 1;
        % x(1) = y0;
        % x(2) = b;
        % x(3) = a;
        % x(4) = xr;
        % x(5) = xl;
        hold on
        xrange = data.data(length(data.data),1):0.001:data.data(2,1);
        %plot(xrange,F(x,xrange),'g-');
    end
    
    %% nonlinear fit
    F1 = @(x,xdata)x(1)+x(2)*xdata-((16/pi)*x(3))*sqrt(3)*(x(4)-x(5))*(xdata-0.5*(x(4)+x(5)))./((4*(xdata-0.5*(x(4)+x(5))).^2+3*(x(4)-x(5)).^2).^2);
    x(1,4)=data.data(minIndx,1);
    x(1,5)=data.data(maxIndx,1);
    [beta,r,J,COVB,mse] = nlinfit(xdataFMR,ydataFMR,F1,x);
    
    if plots>1
        hold on
        plot(xrange,F1(beta,xrange),'b-');
        xlabel('Field (Oe)')
        ylabel('Absorption Derivative (A.U.)')
        if savePlots==1
            h=figure(i);
            fig_file = strrep(files(i).name,'.txt','.jpg');
            saveas(h,fig_file,'jpg');
        end
    end
    
    betaci=nlparci(beta,r,J); %Upper and lower bounds for error
    xl_error = max(abs(beta(4)-betaci(4,1)),abs(beta(4)-betaci(4,2)));
    xr_error = max(abs(beta(5)-betaci(5,1)),abs(beta(5)-betaci(5,2)));
    FMRfreq(i) = data.data(3,4);
    FMRlineWidth(i) = abs(beta(4)-beta(5));
    Linewidtherror(i) = sqrt(xl_error^2+xr_error^2);
    FMRField(i) = (beta(4)+beta(5))/2;
    Fielderror(i) = 0.5*sqrt(xl_error^2+xr_error^2);
    
    FittingOutput{i+3,1} = files(i).name;
    FittingOutput{i+3,2} = FMRfreq(i);
    FittingOutput{i+3,3} = FMRlineWidth(i);
    FittingOutput{i+3,4} = Linewidtherror(i);
    FittingOutput{i+3,5} = FMRField(i);
    FittingOutput{i+3,6} = Fielderror(i);
    
end

%% Kittel & Gilbert Fitting

if Geometry~=2
    %% Kittel fitting
    
    if plots>0
        hold off
        numplot = numplot+1;
        figure(numplot)
        clf;
        errorbar(FMRfreq,FMRField,Fielderror,'b.','MarkerSize',20)
        xlabel('Freqency (GHz)')
        ylabel('Field (Oe)')
    end
    
    if (fixGamma==1 && fixMs==1)
        error('Error: need at least 1 free fitting parameter (Gamma line 23 or Ms line 25)');
    elseif (fixGamma==1)
        if (Geometry == 0)
            F2 = @(xk,xdata)gammaGuess*sqrt(xdata.*(xdata+xk(1)));
        elseif (Geometry == 1)
            F2 = @(xk,xdata)gammaGuess*(xdata-xk(1));
        else
            error('Error: Invalid Geometry (line 23)');
        end
        xktrial = MsGuess;
        
        [beta,r,J,~,~] = nlinfit(FMRField,(1e9*FMRfreq),F2,xktrial);
        xrange = min(FMRField):0.001:max(FMRField);
        
        betaci=nlparci(beta,r,J); %Upper and lower bounds for error
        Gamma = gammaGuess;
        GammaError = 0;
        FourPiMs = beta(1);
        FourPiMsError = max(abs(beta(1)-betaci(1,1)),abs(beta(1)-betaci(1,2)));
    elseif (fixMs==1)
        if (Geometry == 0)
            F2 = @(xk,xdata)xk(1)*sqrt(xdata.*(xdata+MsGuess));
        elseif (Geometry == 1)
            F2 = @(xk,xdata)xk(1)*(xdata-MsGuess);
        else
            error('Error: Invalid Geometry (line 23)');
        end
        xktrial = gammaGuess;
        
        [beta,r,J,~,~] = nlinfit(FMRField,(1e9*FMRfreq),F2,xktrial);
        xrange = min(FMRField):0.001:max(FMRField);
        
        betaci=nlparci(beta,r,J); %Upper and lower bounds for error
        Gamma = beta(1);
        GammaError = max(abs(beta(1)-betaci(1,1)),abs(beta(1)-betaci(1,2)));
        FourPiMs = MsGuess;
        FourPiMsError = 0;
    else
        if (Geometry == 0)
            F2 = @(xk,xdata)xk(1)*sqrt(xdata.*(xdata+xk(2)));
        elseif (Geometry == 1)
            F2 = @(xk,xdata)xk(1)*(xdata-xk(2));
        else
            error('Error: Invalid Geometry (line 23)');
        end
        xktrial = [gammaGuess,MsGuess];
        
        [beta,r,J,~,~] = nlinfit(FMRField,(1e9*FMRfreq),F2,xktrial);
        xrange = min(FMRField):0.001:max(FMRField);
        
        betaci=nlparci(beta,r,J); %Upper and lower bounds for error
        Gamma = beta(1);
        GammaError = max(abs(beta(1)-betaci(1,1)),abs(beta(1)-betaci(1,2)));
        FourPiMs = beta(2);
        FourPiMsError = max(abs(beta(2)-betaci(2,1)),abs(beta(2)-betaci(2,2)));
    end
    
    if plots>0
        hold on
        plot(F2(beta,xrange)*10^-9,xrange,'m-','LineWidth',2);
        if savePlots==1
            h=figure(numplot);
            saveas(h,'Kittel.jpg','jpg');
        end
    end
    
    %% Damping fitting
    
    if plots>0
        numplot = numplot+1;
        figure(numplot)
        clf;
        errorbar((FMRfreq),FMRlineWidth,Linewidtherror,'b.','MarkerSize',20)
        xlabel('Freqency (GHz)')
        ylabel('Line Width (Oe)')
    end
    
    F3 = @(xd,xdata)xd(1)*xdata+xd(2);
    xdtrial = [0,0];
    [beta,r,J,COVB,mse] = nlinfit((FMRfreq),FMRlineWidth,F3,xdtrial);
    betaci=nlparci(beta,r,J); %Upper and lower bounds for error
    xrange = min(FMRfreq):0.01:max(FMRfreq);
    
    if plots>0
        hold on
        plot(xrange,F3(beta,xrange),'m-','LineWidth',2);
        if savePlots==1
            h=figure(numplot);
            saveas(h,'Damping.jpg','jpg');
        end
    end
    
    Damping = beta(1)*Gamma*sqrt(3)/(10^9*2);
    DampingError = max(abs(beta(1)-betaci(1,1)),abs(beta(1)-betaci(1,2)))*Gamma*sqrt(3)/(10^9*2);
    ILB = beta(2);
    ILBerror = max(abs(beta(2)-betaci(2,1)),abs(beta(2)-betaci(2,2)));
    
    AnalysisOutput = {'Gamma', 'GammaError', '4piMs', '4piMsError', 'ILB', 'ILB Error', 'Damping', 'DampingError';
        'MHz/Oe', '', 'G', '', 'Oe', '', '-', '';
        Gamma, GammaError, FourPiMs, FourPiMsError, ILB, ILBerror, Damping, DampingError};
    
    clc;
    fprintf('Gamma   = %6.3e ± %-6.3e\n',Gamma,GammaError)
    fprintf('4piMs   = %6.0f ± %-6.0f\n',FourPiMs,FourPiMsError)
    fprintf('ILB     = %6.4g ± %-6.4g\n',ILB,ILBerror)
    fprintf('Damping = %6.4g ± %-6.4g\n',Damping,DampingError)
        
end

%% Outputs

fileID1 = fopen('FittingResults.dat','w');
formatSpec = '%s\t%g\t%g\t%g\t%g\t%g\n';
[nrows,ncols] = size(FittingOutput);
if fileID1 == -1
    fprintf(2,'\n*** Error: FittingResults.dat already open. Close file and press the any key *** \n');
    pause;
    fileID1 = fopen('FittingResults.dat','w');
    if fileID1 == -1
        error('*** Error: FittingResults.dat is still open. Run the program again ***');
    end
end
fprintf(fileID1,'%s\t%s\t%s\t%s\t%s\t%s\t\n',FittingOutput{1,:});
fprintf(fileID1,'%s\t%s\t%s\t%s\t%s\t%s\t\n',FittingOutput{2,:});
for row = 3:nrows
    fprintf(fileID1,formatSpec,FittingOutput{row,:});
end
fclose(fileID1);

if Geometry~=2
    fileID2 = fopen('AnalysisResults.dat','w');
    if fileID2 == -1
        fprintf(2,'\n*** Error: AnalyssisResults.dat already open. Close file and press the any key *** \n');
        pause;
        fileID2 = fopen('AnalysisResults.dat','w');
        if fileID2 == -1
            error('*** Error: AnalyssisResults.dat is still open. Run the program again ***');
        end
    end
    fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',AnalysisOutput{1,:});
    fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',AnalysisOutput{2,:});
    fprintf(fileID2,'%e\t%e\t%f\t%f\t%g\t%g\t%g\t%g\n',AnalysisOutput{3,:});
    fclose(fileID2);
end
fclose('all');
fprintf('\nDone\n');