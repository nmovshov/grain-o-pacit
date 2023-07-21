function y=oParams(varargin)
% Returns a structure with model parameters.

global si 

% Default parameters
y.coreMass=10.0*si.earth_mass; % Mass of protoplanet's core (kg)
y.envMass=0.1*si.earth_mass; % Mass of protoplanet's atmosphere (kg)
y.gasMass=2.3*(si.g/si.mole)/si.avogadro; % Mass of a gas molecule (kg)
y.sigmaGas=1e-15*si.cm^2; % Collision cross section for gas (m^2)
y.silicateDensity=2800*si.kg/si.m^3; % Density of grain matter (kg/m^3)[1]
y.monomerSize=1e-6*si.m; % Radius (m)[2]
y.binSpacingParameter=3; % Dimensionless[3]
y.nBins=30; % How many size bins to use
y.stickCoef=1; % Sticking probability in grain-grain collision
y.maxTimeStep=1e6*si.s; % Maximum allowed time step (s)
y.accretionRate=0.0e-1*si.g/si.cm^2/si.s; % Mass flux of solids from outside the protoplanet
y.dust2gasRatio=0.01; % Ratio of dust to gas in solar nebula

% Override default parameters with input parameters
propertyArgIn=varargin;
while length(propertyArgIn)>=2
    prop=propertyArgIn{1};
    value=propertyArgIn{2};
    propertyArgIn=propertyArgIn(3:end);
    if any(strcmp(prop,fieldnames(y)))
        y.(prop)=value;
    else
        error('Unrecognized parameter input.')
    end
end

% [1] Morris always assumes a density of 2.8 gm/cm^3 as the typical matter
% density of silicate material.
%
% [2] We need to assume some monomer size in the aggregates that are the dust
% grains. This paramter should be experimented with. Suggested values are
% 0.01 to 10 microns. (ref:podolak03.pdf) (Morris calls this the "scale
% factor", or scfact, in old programs.)
%
% [3] Grains are put into discrete size bins characterized by radii
% a_n=a_0*2^(alpha*n). I call alpha the bin spacing parameter. We now use
% alpha=(1/3).

end