function jack(startTime,endTime,varargin)
%JACK   This is the part of program that tracks the time evolution of the
% grain number density (the dust.nDensity array), accounting for
% sedimentation of grains through the atmosphere and coagulation of grains
% into larger grains due to collisions.

global si params atm dust

%% Local variables
time=startTime;
icount=0;
nd=dust.nDensity;
fluxFlag=true; coagFlag=true; sourceFlag=true; msgFlag=true;
sourceEqualDistribution=false;

%% Set qualifiers
qualifiers=varargin;
while length(qualifiers)>=1
    qlf=qualifiers{1};
    qualifiers=qualifiers(2:end);
    switch qlf
        case 'noflux'
            fluxFlag=false;
        case 'nocoag'
            coagFlag=false;
        case 'nosource'
            sourceFlag=false;
            dust.source=0*si.kg/si.m^3/si.s;
        case 'equisource'
            sourceEqualDistribution=true;
        case 'nomessages'
            msgFlag=false;
        otherwise
            error('Unrecognized qualifier input.')
    end
end

%% Main computation loop
while time<endTime
%%  < Find the density rate of change >
    % First, the flux.
    if fluxFlag
        dndtFlux=flux();
    else
        dndtFlux=zeros(size(nd))*si.m^-3*si.s^-1;
    end
    % Then Coag.
    if coagFlag
        dndtCoag=coag();
    else
        dndtCoag=zeros(size(nd))*si.m^-3*si.s^-1;
    end
    % Next the source term.
    if sourceFlag
        dndtSource=zeros(size(nd))*si.m^-3*si.s^-1;
        if sourceEqualDistribution
            for k=1:length(atm.Z)
                dndtSource(k,1:dust.aBins)=(dust.source(k)/dust.aBins)./...
                    dust.massBin(1:dust.aBins);
            end
        else
            dndtSource(:,1)=dust.source/dust.massBin(1); % Put input mass into smallest bin
        end
    else
        dndtSource=zeros(size(nd))*si.m^-3*si.s^-1;
    end
    % And all together.
    dndt=dndtCoag+dndtFlux+dndtSource;

%%  < Determine an appropriate time step >
    ind=find(nd>repmat(atm.dV,1,params.nBins).^-1);
    if isempty(ind)
        dt=params.maxTimeStep;
    else
        x=max(abs(dndt(ind))./nd(ind));
        if (double(x))
            dt=0.01/x;
        else
        dt=params.maxTimeStep;
        end
    end
    if dt>params.maxTimeStep
        dt=params.maxTimeStep;
    end
    if time+dt>endTime, dt=endTime-time; end

%%  < Update the nDensity array >
    nd=nd+dndt*dt;
    nd(nd<0*si.m^-3)=0*si.m^-3;
    dust.nDensity=nd;

%%  < Update the flux capacitor >
    dust.fluxCapacitor=dust.fluxCapacitor+dust.dFluxCapacitordt*dt;

%%  < Advance the time variable and iteration counter >
    time=time+dt;
    icount=icount+1;
 
%%  < Update progress bar and print messages >
    if ~rem(icount,1000) && msgFlag
    fprintf('Time is %g seconds.\nTime step is %g.\n\n',...
        double(time),double(dt));
    end

end % of main loop

%% Finishing up messages.
if msgFlag
fprintf('END OF RUN (%g - %g s).\n',startTime,endTime)
end

end % of function jack


%% --- Helper Functions ---

function y=coag()
% Returns an array that holds the time rate of change of grain density in
% each size bin and zone, in [m^-3 s^-1].
% The rate of change in number density of the grains due to coagulation is
% given by (see podolak03.pdf, first two terms in Eq. 18)
% \frac{\partial{}n(m,t)}{\partial{}t}=...
% (1/2)*\int_{0}^{m}K(m',m-m')*n(m',t)*n(m-m',t)dm'-...
% n(m,t)*\int_{0}^{\infty}K(m,m')*n(m',t)dm'.
% This is approximated here for a discrete size distribution.

global si dust

% Local variables
nd=dust.nDensity;
ms=dust.massBin;
kernel=dust.kernel;
redist=dust.redist;
y=zeros(size(nd))*si.m^-3*si.s^-1;
nZones=size(nd,1);
nBins=size(nd,2);
aBins=dust.aBins;

for j=1:nZones
    for k=1:nBins
        % Calculate dndt for the jth zone and kth bin
        dndt=0*si.m^-3*si.s^-1;
        for m=1:k-1
            for n=1:m
                dndt=dndt+...
                    kernel(j,m,n)*redist(m,n,k)*nd(j,m)*nd(j,n)*...
                    (1-0.5*(m==n));
            end
        end
        for m=1:nBins
            dndt=dndt-...
                kernel(j,k,m)*nd(j,k)*nd(j,m)*(1-redist(k,m,k));
        end
        y(j,k)=dndt;
    end
    % Recycle "spill-out" back to active bins.
    for k=aBins+1:nBins
        y(j,aBins)=y(j,aBins)+y(j,k)*ms(k)/ms(aBins);
        y(j,k)=0*si.m^-3*si.s^-1;
    end
end

end % of function coag


function y=flux()
% Returns an array that holds the time rate of change of grain density, in
% each bin and zone, due to the net flux of grains of the same bin falling
% into and out of the zone.

global si atm dust

% Local variables
nd=dust.nDensity;
vs=dust.vSed;
R=atm.R;
y=zeros(size(nd))*si.m^-3*si.s^-1;
N=length(R);

% Uppermost zone - only flux out (flux in included in source term)
dndt=-nd(1,:).*vs(1,:)*3*R(2)^2/(R(1)^3-R(2)^3);
y(1,:)=dndt;
% The other zones - flux in and out
for k=2:N-1
    fluxIn=nd(k-1,:).*vs(k-1,:)*4*pi*R(k)^2;
    fluxOut=nd(k,:).*vs(k,:)*4*pi*R(k+1)^2;
    volume=atm.dV(k);
    dndt=(fluxIn-fluxOut)/volume;
    y(k,:)=dndt;
end
% Fallout zone - only inward flux*
% * The fallout "zone" collects all the mass that has sedimented out of
%   the last zone into a "flux capacitor", to allow a check on the
%   conservation of total mass.
fluxIn=nd(N-1,:).*vs(N-1,:)*4*pi*R(N)^2;
fluxOut=0*fluxIn;
volume=1*si.m^3;% A fictitious volume!
dndt=(fluxIn-fluxOut)/volume;
dust.dFluxCapacitordt=dndt;

end % of function flux.
