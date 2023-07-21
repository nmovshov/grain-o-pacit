function y=oTotalMass()
% Returns the total mass in the atmosphere, in kg.

global dust atm si params

nBins=params.nBins;
massDensity=dust.nDensity(:,1:nBins).*...
    repmat(dust.massBin(1:nBins),size(dust.nDensity,1),1);
volume=4/3*pi*(atm.R(1:end-1).^3-atm.R(2:end).^3);
mass=massDensity.*repmat(volume,1,nBins);
y{1}=mass;
mFallout=sum(dust.fluxCapacitor(1:nBins).*...
    dust.massBin(1:nBins)*(1*si.m^3)); % fictitious volume of 1
y{2}=mFallout;

end