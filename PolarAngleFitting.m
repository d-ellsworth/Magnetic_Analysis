%% Polar angle linewidth fitting
% Praveen/Houchen/David 14-Mar-2014
%
% Fits data in Data.dat with polar angle dependent gilbert damping eqn
% FittingResults.dat should match output of FMRdatafittingFinal.m (from
% polar angle data fitting) with col(1) (filename) replaced with angle
%
% Reads 4piMs and gamma from AnaysisResults.dat (from freq dependent
% analysis). To enter values manually comment out line 33-35 and uncomment
% line 36-37.
%
% Outputs fitted alpha and Theta vs dH curve to ThetaVsDh.dat
%
% Generates Hfmr & phi for given Theta
% Outputs Theta (range), Hfmr (from calc), and Phi (from calc) to
% ThetaVsHmfr.dat
%
%

clc;
clear all;

%% Input parameters

step = 1;
ThetaC = 0:step:90;

LinewidthData = importdata('FittingResults.dat','\t',1);
Theta = (LinewidthData.data(:,1))';
LineWidth = (LinewidthData.data(:,3))';
f = mean(LinewidthData.data(:,2));

% Comment/uncomment below if using data import or not
SampleData = importdata('AnalysisResults.dat','\t',1);
fourpiMs = SampleData.data(3)/1000;
gamma = SampleData.data(1)*10^-6;
% fourpiMs = 15.5;
% gamma = 2.8;
Ha = 0;             % in kOe

Hini = f/gamma+fourpiMs-Ha;

%% Find Hfmr & Polarangle

% equations to solve
% x(1)is external magnetic field, x(2)is angle phi
funField = @(x,f,gamma,Ha,fourpiMs,theta)[x(1)*sin((theta-x(2))/180*pi)+.5*(fourpiMs-Ha)*sin(2*x(2)/180*pi),...
    (f/gamma)^2-(x(1)*cos((theta-x(2))/180*pi)-fourpiMs*cos(x(2)/180*pi)^2+Ha*cos(x(2)/180*pi)^2)*(x(1)*cos((theta-x(2))/180*pi)+(Ha-fourpiMs)*(cos(x(2)/180*pi)^2-sin(x(2)/180*pi)^2))];

% for data
for i=1:length(Theta)
    thetai=Theta(i);
    x0=[Hini,thetai];
    y=fsolve(@(x) funField(x,f,gamma,Ha,fourpiMs,thetai),x0);
    Hfmr(i)=y(1);
    Phi(i)=y(2);
    clear y;
end

% for calculation
for i=1:size(ThetaC,2)
    thetai=ThetaC(i);
    x0=[Hini,thetai];
    y=fsolve(@(x) funField(x,f,gamma,Ha,fourpiMs,thetai),x0);
    HfmrC(i)=y(1);
    PhiC(i)=y(2);
    clear y;
end


%% Find C
% for data
for i=1:length(Theta)
    Hx = (Hfmr(i)*cos(((Theta(i)-Phi(i))*pi)/180)-fourpiMs*cos((2*Phi(i)*pi)/180));
    Hy = (Hfmr(i)*cos(((Theta(i)-Phi(i))*pi)/180)-fourpiMs*(cos((Phi(i)*pi)/180))^2);
    c(i) = cos(((Theta(i)-Phi(i))*pi)/180)-((Hfmr(i)*(Hx+3*Hy)*(sin(((Theta(i)-Phi(i))*pi)/180)^2))/(Hx*(Hx+Hy)));
    OneOverC(i) = 1/c(i);
end

% for calculation
for i=1:length(ThetaC)
    Hx = (HfmrC(i)*cos(((ThetaC(i)-PhiC(i))*pi)/180)-fourpiMs*cos((2*PhiC(i)*pi)/180));
    Hy = (HfmrC(i)*cos(((ThetaC(i)-PhiC(i))*pi)/180)-fourpiMs*(cos((PhiC(i)*pi)/180))^2);
    C(i) = cos(((ThetaC(i)-PhiC(i))*pi)/180)-((HfmrC(i)*(Hx+3*Hy)*(sin(((ThetaC(i)-PhiC(i))*pi)/180)^2))/(Hx*(Hx+Hy)));
    OneOverCC(i) = 1/C(i);
end

%% Gilbert Fit

F1 = @(x,xdata)x(1)*2*xdata*2*f/gamma;
x=0.003;                                % error estimate
[beta,r,J,COVB,mse] = nlinfit(OneOverC,LineWidth,F1,x);

Alpha = beta*10^-3*sqrt(3)              %damping parameter
AlphaErr=nlparci(beta,r,J);             %Upper and lower bounds for error
AlphaError=abs(AlphaErr(1)-AlphaErr(2))*10^-3*sqrt(3)

% dH values from fit
for j=1:length(OneOverCC)
    dH(j)=OneOverCC(j)*2*0.0056*f*10^3/gamma;
end

%% Plots
figure(1)
plot(ThetaC,dH,'g-');
title('Linewidth vs Theta');
hold on;
plot(Theta,LineWidth,'o');

figure(2)
plot(OneOverC,F1(beta,OneOverC),'g-');
title('Linewidth vs 1/C');
hold on;
plot(OneOverC,LineWidth,'o');

%% Outputs
file1=fopen('ThetaVsHmr.dat','w');
fprintf(file1,'%s\t%s\t%s\r\n','Theta','Hfmr','Phi');
fprintf(file1,'%f\t%f\t%f\r\n',[ThetaC;HfmrC;PhiC]);
fclose('all');

file2=fopen('ThetaVsDh.dat','w');
Header = {'Alpha','AlphaError','Theta';Alpha,AlphaError,'dH'};
fprintf(file2,'%s\t%s\r\n',Header{:,:});
fprintf(file2,'%f\t%f\r\n',[ThetaC;dH]);
fclose('all');