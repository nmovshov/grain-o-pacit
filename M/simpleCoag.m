%% SIMPLECOAG
% SIMPLE COAGLUATION MODEL
% The purpose of this script is to calculate a charecteristic grain
% size in a steady state sedimentation scenario. Basically, I try to find 
% one grain size in each zone that equallizes the typical times of
% sedimentation and coagulation. Details are in the code. See also
% SIMPLECOAG.PDF.
% clear all
% close all
% clc
%% --- SETUP PHASE ---
%% physical constants
G=6.673e-11; % Universal gravitational constant (m^3 kg^-1 s^-2)
Navo=6.02e23; % Avogadro's number (1/mole)
KB=1.38e-23; % Boltzmann constant (m^2*kg*s^-2*K^-1)
EM=5.9742e24; % Earth mass (kg)
%% model parameters
%F0=4e-12; % Supply of dust into the atmosphere (kg/m^2/s)*
%F0=3e-14; % Supply of dust into the atmosphere (kg/m^2/s)**
F0=1.4e-18; % Supply of dust into the atmosphere (kg/m^2/s)
Mcore=6.9e25; % Mass of the planet's core (kg)
mu=2.3e-3; % mean molecular weight (kg/mole)
sigma=1e-19; % Molecular collison cross section (m^2)
rhograin=2800; % grain matter density (kg/m^3)
a=logspace(-8,-3); % grain sizes
%* Based on 1e-8 Earth masses a year
%** Based on 0.01 of the gas mass in the upper layer
%% read in atmosphere data
% The atmosphere data is in an input file, usually 'atm.dat', with a Morris
% like format: The first 19 lines specify various parameters, used as
% configuration of a more complete coagulation scheme. Here I ignore all of
% these as header lines. The rest of the file is a table of atmospheric
% data in the form of [zone density pressure temperature gam vconv]. Of
% these I need zone density and temperature. Values in Morris files are in
% cgs.
[z d t]=textread('atm.dat','%f%f%*f%f%*[^\n]',...
    'headerlines',19);
z=z/100;
d=d*1000;
Z=z(2:end);
D=d(2:end);
T=t(2:end);

%% --- CALCULATION PHASE ---
%% grain independent properties
l=abs(diff(z)); % layer thickness (m) l(1) uses values of Z(2)
g=G*Mcore./Z.^2; % acceleration of gravity (m/s^2)
ngas=D*Navo/mu; % gas number density (1/m^3)
vth=sqrt((8*KB*T)/(pi*mu/Navo)); % thermal speed of gas (m/s)
mfp=1./(sqrt(2)*ngas*sigma); % mean free path in gas (m)
%% drag forces (chose one)
Knud=[];
for k=1:length(a)
    Knud=[Knud mfp/a(k)]; % Knudsen number
end
A=1.249; B=0.42; C=0.87;
psi=1+Knud.*(A+B.*exp(-C./Knud));
EpsDrag=[];
for k=1:length(a)
    EpsDrag=[EpsDrag vth.*D*a(k)^2]; % Epstein drag per speed (N m^-1 s)
end
eta=8.6e-7*sqrt(T); % Viscosity (kg/m/s)
StokesDrag=[];
for k=1:length(a)
    StokesDrag=[StokesDrag 6*pi.*eta*a(k)]; % Stks drag per speed (N/m s)
end
psiDrag=StokesDrag./psi; % interpolated drag per speed (N/m s)
%% sedimentation speed and time
vsed=[];
for k=1:length(a)
    vsed=[vsed (rhograin*g*a(k))./(D.*vth)]; % sedimentation speed (m/s)*
end
%* Based on Epstein drag (Knud>>1)
tsed=[];
for k=1:length(a)
    tsed=[tsed l./vsed(:,k)]; % sedimentation time (s)
end
%% coagulation time
ngrain=[];
for k=1:length(a)
    ngrain=[ngrain (3*F0*Z(1)^2)./...
        (4*pi*rhograin*Z.^2.*vsed(:,k)*a(k)^3)]; % grain number density (m^-3)
end
mfpgrain=[];
for k=1:length(a)
    mfpgrain=[mfpgrain (1)./... % mean distance between grains (m)
        (sqrt(2)*4*pi*a(k)^2*ngrain(:,k))];
end
mgrain=4/3*pi*rhograin*a.^3; % grain mass (kg)
vBrown=[];
for k=1:length(a)       % speed of Brownian motion of grains (m/s)
    vBrown=[vBrown sqrt((8*KB*T)/(pi*mgrain(k)))];
end
tcoag=mfpgrain./vBrown; % coagulation time (s)
%% analytic solution
% Assuming Epstein drag only equal times of sedimentation and coagulation
% can be solved to obtain a grain size.
astar=((l.*D.^2.*(vth.^2)*(F0*4*pi*Z(1)^2))./...
    (rhograin^(7/2)*g.^2.*Z.^2).*...
    (1.5*sqrt(3*KB*T))/(pi^2)).^(2/9);
vsedstar=(rhograin*g.*astar)./(D.*vth);
tsedstar=l./vsedstar;
ngrainstar=(3*F0*Z(1)^2)./(4*pi*rhograin*vsedstar.*astar.^3.*Z.^2);
mfpgrainstar=(1)./(2^(2.5)*pi*ngrainstar.*astar.^2);
mgrainstar=4/3*pi*rhograin*astar.^3;
vBrownstar=sqrt((8*KB*T)./(pi*mgrainstar));
tcoagstar=mfpgrainstar./vBrownstar;
%% --- VISUALISATION PHASE ---
%% grain size vs. height plot (analytic)
f1=figure;
l1=plot(Z,astar);
ax1=gca;
xlabel(ax1,'ZHeight [m]')
ylabel(ax1,'grain size [m]')
set(l1,'color','k','linewidth',1)
strTitle=sprintf('F0=%g, M_{core}=%g',F0,Mcore);
title(strTitle)

f2=figure;
for k=1:4
ax2=subplot(2,2,k);
hold on
l1=plot(a,tsed(10*k-9,:),'r+-','displayname','t_{sed}');
l2=plot(a,tcoag(10*k-9,:),'b+-','displayname','t_{coag}');
l3=plot([astar(10*k-9) astar(10*k-9)],...
    [tsed(10*k-9,1) tsed(10*k-9,end)],'k:','displayname','a^{*}');
[c ci]=min(abs(tcoag(10*k-9,:)-tsed(10*k-9,:)));
xlim(ax2,[a(ci-3) a(ci+3)]);
ylim(ax2,[tcoag(10*k-9,ci-2) tcoag(10*k-9,ci+2)]);
legend([l1 l2 l3])
title(['R = ',num2str(Z(10*k-9),'%2.1e'),' [m]'])
end
