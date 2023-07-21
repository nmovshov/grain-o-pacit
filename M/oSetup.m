function oSetup(varargin)
% This is the setup module of program OPAC. It is responsible for creating
% the static variables that will be used throughout the program. These
% include the properties of the model atmosphere that are read from a file.
% The setup routine also creates all variables that are derived via
% one-time calculations on the input and parameters. Most important are the
% collision kernel and mass redistribution coefficients (explained further
% in the code). All the variables are saved in a small number of global
% structures that are accessible to all program modules.

global si params atm dust

%% File names and non-default configuration parameters
propertyArgIn=varargin;
while length(propertyArgIn)>=2
    prop=propertyArgIn{1};
    value=propertyArgIn{2};
    propertyArgIn=propertyArgIn(3:end);
        switch prop
            case 'kernel'
                kernelFlag=value;
            case 'atmos'
                atmFile=value;
            case 'index'
                indexFile=value;
            case 'source'
                sourceFile=value;
            case 'initial'
                initDistFile=value;
            case 'massbin'
                massBinFlag=value;
            otherwise
                error('Unrecognized parameter input.')
        end
end

%% The model atmosphere 
% The model atmosphere files have the format:
% [radius] [density] [pressure] [temperature] [...]
% This program requires the first, second, and fourth columns. Note: all
% values are in cgs.
fid=fopen(atmFile);
while(1) % Skip header lines
    C=fgetl(fid);
    if strcmp(strtrim(C),'--- BEGIN DATA ---'), break, end
end
C=textscan(fid,'%f%f%*f%f%*[^\n]');
fclose(fid);
% Input properties
atm.R=C{1}*si.cm; % Radius [m]
atm.Dd=C{2}*si.g/si.cm^3; % Gas density [kg/m^3]
atm.Tt=C{3}*si.kelvin; % Gas temperature [K]
% Derived properties
atm.Z=(atm.R(1:end-1)+atm.R(2:end))/2; % Layer height [m]*
% * Note. There are length(R)-1 'Zones'. Zone(i) is the volume between R(i)
% and R(i+1). All the properties of Zone(i), e.g. temperature, will be in
% the ith place of the appropriate vectors, and averaged from the values
% on the boundaries.
atm.dR=abs(diff(atm.R)); % Zone thickness [m]
atm.dV=4/3*pi*(atm.R(1:end-1).^3-atm.R(2:end).^3); % Zone volume [m^3]
atm.D=(atm.Dd(1:end-1)+atm.Dd(2:end))/2; % Gas density [kg/m^3]
atm.T=(atm.Tt(1:end-1)+atm.Tt(2:end))/2; % Gas temperature [K]
atm.g=grav(); % Acceleration of gravity [m/s^2]
    function y=grav()
        bigEm=params.coreMass+params.envMass-sum(atm.D.*atm.dV);
        diEm=atm.D.*atm.dV;
        for k=1:length(atm.Z)-1
            y(k)=si.gravity*(bigEm+sum(diEm(k+1:end)))/atm.Z(k)^2;
        end
        y(length(atm.Z))=si.gravity*bigEm/atm.Z(end)^2;
    end % of nested function grav
atm.vThermal=sqrt((8*si.boltzmann.*atm.T)/(pi*params.gasMass)); % Thermal velocity of the gas [m/s]
atm.nDensity=atm.D/params.gasMass; % Gas number density [m^-3]
atm.meanFreePath=1./(sqrt(2)*params.sigmaGas*atm.nDensity); % Mean free path between gas molecules [m]
atm.viscosity=8.6e-6*sqrt(atm.T)*si.g/si.cm/si.second/si.K^0.5; % Gas viscosity [kg m^-1 s^-1]

%% The dust structure
%% Create the size bins and corresponding mass bins
if ~exist('massBinFlag','var'), massBinFlag='geometric'; end
[dust.sizeBin dust.massBin]=computeBins(massBinFlag);
%% Calculate the largest allowed active bin
% We use the last few bins as "recycle bins", to make sure no mass is lost
% due to coagulation into too large particles. We need the largest "active"
% bin to be half the mass of the largest bin. (With a spacing parameter of
% 3 aBins=nBins-1.)
dust.aBins=find((2*dust.massBin-dust.massBin(params.nBins))>...
    dust.massBin(params.nBins)*1e-3,1,'first')-1;
%% Calculate sedimentation speed
dust.vSed=computeSedimentationSpeed;
    function y=computeSedimentationSpeed()
        % Create an array to hold the sedimentation speed of a grain in the
        % atmosphere. In this version the atmosphere is quiescent, i.e.,
        % there is no convection.
        y=preal(ones(length(atm.Z),params.nBins));
        for j=1:length(atm.Z)
            for k=1:params.nBins
                %y(j,k)=...
                %    (atm.g(j)*params.silicateDensity*dust.sizeBin(k))/...
                %    (atm.D(j)*atm.vThermal(j)); % Epstein style.
                kn=atm.meanFreePath(j)/dust.sizeBin(k); % Knudsen number
                psi=1+kn*(1.249+0.42*exp(-0.87/kn)); % interpolation factor
                y(j,k)=(dust.massBin(k)*atm.g(j)*psi)/...
                    (6*pi*dust.sizeBin(k)*atm.viscosity(j))+...
                    params.accretionRate/params.dust2gasRatio/atm.D(j);
            end
        end
    end % of nested function computeSedimentationSpeed
%% Initialize density array.
% The nDensity array holds the number density of dust grains in each bin
% and zone. nDensity(k,j) is the number, per unit volume, of j-bin grains
% in the kth zone.
dust.nDensity=preal(zeros(length(atm.Z),params.nBins));
if ~(exist('initDistFile','var')==1), initDistFile='new'; end
if strcmp(initDistFile,'new')
    dust.nDensity(:,1)=...
        params.dust2gasRatio*atm.D/dust.massBin(1);
    dust.nDensity(:,2:end)=0*si.m^-3;
else
    fid=fopen(initDistFile);
    while(1) % Skip header lines
        C=fgetl(fid);
        if strcmp(strtrim(C),'--- BEGIN DATA ---'), break, end
    end
    for k=1:length(atm.Z)
        C=fgetl(fid);
        dust.nDensity(k,:)=str2num(C)*si.cm^-3;
    end
    fclose(fid);
end
        
%% Create the collision kernel.
% The 'kernel' array holds the probabilities of collision, as a result of
% Brownian motion and 'overtaking'. kernel(i,j,k) is the probability of
% collision between an j-bin particle and a k-bin particle in the ith zone.
% This is unnormalized probability, i.e., kernel*density*density is the
% number of collisions per unit volume per unit time.
dust.kernel=computeKernel();
if exist('kernelFlag','var')
    if strcmp(kernelFlag,'fixed')
        dust.kernel=dust.kernel(1)*ones(size(dust.kernel));
    else
        error('Unknown kernel flag')
    end
end
%% Create coefficients of mass redistribution.
% The 'redist' array holds the probability of mass redistibution options
% after collision. redist(i,j,k) is the probability that after a
% collision of an i-bin particle with a j-bin particle, a k-bin particle
% will be created.
dust.redist=computeRedistProbs();
%% Fallout
% The fluxCapacitor field holds the mass that sediments out of the area of
% interest.
dust.fluxCapacitor=zeros(1,params.nBins)*si.m^-3;
dust.dFluxCapacitordt=zeros(1,params.nBins)*si.m^-3/si.s;
%% Index of refraction
% Read in a table of refractive indices for different wavelengths
fid=fopen(indexFile);
while(1) % Skip header lines
    C=fgetl(fid);
    if strcmp(strtrim(C),'--- BEGIN DATA ---'), break, end
end
C=textscan(fid,'%f%*f%*f%f%f');
dust.index=[C{:,1}*si.micron C{:,2} C{:,3}];
fclose(fid);
%% Source term (input of solids from planetesimal ablation)
% Read in a table of mass input at various heights
dust.source=preal(zeros(length(atm.Z),1));
if strcmp(sourceFile,'none')
    dust.source(1)=params.accretionRate/atm.dR(1);
    dust.source(2:end)=0*si.kg/si.m^3/si.s;
else
    fid=fopen(sourceFile);
    while(1) % Skip header lines
        C=fgetl(fid);
        if strcmp(strtrim(C),'--- BEGIN DATA ---'), break, end
    end
    C=textscan(fid,'%f%f');
    fclose(fid);
    table=[C{:,1}*si.cm C{:,2}*si.g/si.cm^3/si.s];
    for k=1:length(atm.Z)
        dust.source(k)=interp1(double(table(:,1)),double(table(:,2)),...
            double(atm.Z(k)),'linear',table(1,2))*si.kg/si.m^3/si.s;
    end
    dust.source(1)=dust.source(1)+params.accretionRate/atm.dR(1);
end
%% Clear unnecessary stuff
dust.massBin(end)=[]; % Used for sake of brevity in computeRedistProb.
dust.sizeBin(end)=[]; % --

end % End of function oSetup


%% --- Helper Functions ---

function y=computeRedistProbs()
% Calculates the probability of mass redistribution options after a
% collision between grains.

global params dust

n=params.nBins;
y=preal(zeros(n,n,n));
mass=(dust.massBin);
if 2*mass(1)<=mass(2)
    y(1,1,1)=(mass(2)-2*mass(1))/(mass(2)-mass(1));
end
for i=1:n
    for j=1:n
        for k=2:n
            if ((mass(i)+mass(j))>mass(k-1))&&((mass(i)+mass(j))<=mass(k))
                y(i,j,k)=(mass(i)+mass(j)-mass(k-1))/...
                    (mass(k)-mass(k-1));
                continue
            end
            if ((mass(i)+mass(j))>mass(k))&&((mass(i)+mass(j))<=mass(k+1))
                y(i,j,k)=(mass(k+1)-mass(i)-mass(j))/...
                    (mass(k+1)-mass(k));
                continue
            end
        end
    end
end
end % of function computeRedistProbs


function y=computeKernel() % WARNING: a "gray box" function
% Create the collision kernel.

global si params atm dust

bin=dust.sizeBin;
mass=dust.massBin;
mu=params.gasMass;
gama=params.stickCoef;
vth=atm.vThermal;
vSed=dust.vSed;
numZones=length(atm.Z);
difuse=preal(ones(numZones,params.nBins));
vDust=preal(ones(size(difuse)));
delta=preal(ones(size(difuse)));
for j=1:numZones;
    for k=1:params.nBins
        difuse(j,k)=(3*si.boltzmann*atm.T(j))/...
            (4*pi*bin(k)^2*atm.nDensity(j)*mu*vth(j));
        vDust(j,k)=vth(j)*sqrt(mu/mass(k));
        lB=(8*difuse(j,k))/(pi*vDust(j,k));
        delta(j,k)=(sqrt(2)/(6*bin(k)*lB))*...
            ((2*bin(k)+lB)^3-(4*bin(k)^2+lB^2)^(3/2))-...
            (2*bin(k));
    end
end
y=preal(ones(numZones,params.nBins,params.nBins));
for i=1:numZones
    for j=1:params.nBins
        for k=1:params.nBins
            D=(difuse(i,j)+difuse(i,k))/2;
            del=sqrt(delta(i,j)^2+delta(i,k)^2);
            v=sqrt(vDust(i,j)^2+vDust(i,k)^2);
            r=(bin(j)+bin(k))/2;
            P1=8*pi*r*D/...
                (r/(r+del/2)+(4*D)/(r*v));
            P2=4*pi*r^2*abs(vSed(i,j)-vSed(i,k));
            y(i,j,k)=gama*(P1+P2);
        end
    end
end
end % of function computeKernel


function [RB MB]=computeBins(flag)
% Create discrete bins for mass and radius.

global params

switch flag
    case 'geometric'
        RB=preal(ones(1,params.nBins+1));
        MB=preal(ones(1,params.nBins+1));
        for k=1:params.nBins+1 % (+1 to facilitate creation of redist)
            RB(k)=params.monomerSize*2^(k/params.binSpacingParameter);
            V=4*pi*RB(k)^3/3;
            MB(k)=V*params.silicateDensity;
        end
    case 'linear'
        for k=1:params.nBins+1
            RB(k)=params.monomerSize*k^(1/3);
            V=4*pi*RB(k)^3/3;
            MB(k)=V*params.silicateDensity;
        end
    otherwise
        error('Unknown mass bin flag')
end
end % of function computeBins