%OPAC   Driver script for the opacity program

%% Clear and set the workspace
clear all
close all
clc
global si params atm dust

%% Set up the model
physunits off
si=setUnits;
params=oParams('nBins',35,'maxTimeStep',1e6*si.s,'monomerSize',1e-6*si.m,...
    'binSpacingParameter',3,'coreMass',11.492*si.earth_mass,...
    'envMass',3.288E-01*si.earth_mass);
oSetup('atmos','atm2321.dat','index','ind.dat','source','tes2321.dat');
pause(1)

%% 2. Run model and save results
% % Create a file to hold output
% timeStamp=now;
% filename=['opacOut',num2str(timeStamp),'.mat'];
% save(filename,'si','params','atm','dust');
% 
% % Set a list of destination times
% run.markTimes=[0 3e9 1e10 3e10];
% 
% % Run the model and save results at specified times
% run.nd{1}=dust.nDensity;
% radt
% run.op{1}=atm.opacity;
% tt=run.markTimes;
% for k=2:length(tt)
%     jack(tt(k-1),tt(k))
%     pause(1)
%     radt
%     run.nd{k}=dust.nDensity;
%     run.op{k}=atm.opacity;
%     save('-append',filename,'run')
% end