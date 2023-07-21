function radt()
%RADT   This is the part that calculates the dust's contribution to the
% opacity of the atmosphere. The Rosseland mean opacity and the optical
% depth are calculated for each layer. This function uses Cuzzie's MIE.for
% subroutine for the computation of extinction efficiencies. Values for the
% refractive index of the grain material are taken from a table.

global si atm dust

%% Local variables
h=si.planck;
kB=si.boltzmann;
c=si.speed_of_light;
mnu=21; % Number of frequency (wave length) points for integration
table=double(dust.index);

%% For each zone in the atmosphere...
for z=1:length(atm.Z)
    
%% Define the range of relevant frequencies/wavelengths
    T=atm.T(z);
    nupeak=5.879e10*si.Hz/si.K*T; % Location of peak in Planck's function
    numin=nupeak/10;
    numax=nupeak*3.3;
    nu=linspace(numin,numax,mnu);
    lambda=c./nu;

%% Calculate the specific cross-section for all wave lengths and size bins
% (Using the prescription provided by Cuzzi.)
    crosec=preal(zeros(length(nu),length(dust.sizeBin)));
    for k=1:length(nu)
        cmplx=interp1(table(:,1),[table(:,2) table(:,3)],...
            double(lambda(k)),'linear');
        nr=cmplx(1); ni=cmplx(2);
        [x ind]=max(table(:,1));
        if double(lambda(k))>x
            nr=table(ind,2);
            ni=table(ind,3);
        end
        [x ind]=min(table(:,1));
        if double(lambda(k))<x
            nr=table(ind,2);
            ni=table(ind,3);
        end
        %nr=1.5; ni=0.01;
        crosec(k,:)=mie(nu(k),dust.sizeBin,nr,ni);
    end

%% Define the extinction coefficient and opacity
    K=crosec.*repmat(dust.nDensity(z,:),length(nu),1);
    kappa=K/atm.D(z);
    
%% Calculate the Rosseland mean opacity
    % Simpson-integrate inverse kappa weighted by dBdT
    dBdT=@(n)(2*h^2*c^-2*kB^-1*T^-2)*n.^4.*exp(h*n/kB/T).*(exp(h*n/kB/T)-1).^-2;
    weight=dBdT(nu);
    invkap=sum(kappa(1,:))^-1*weight(1);
    for k=2:2:mnu-1
        invkap=invkap+4*sum(kappa(k,:))^-1*weight(k);
    end
    for k=3:2:mnu-2
        invkap=invkap+2*sum(kappa(k,:))^-1*weight(k);
    end
    invkap=invkap+sum(kappa(mnu,:))^-1*weight(mnu);
    invkap=invkap*(nu(2)-nu(1))/3;    
    normal=quad(dBdT,numin,numax);
    invkap=invkap/normal;

    Ross(z)=invkap^-1;
    tau(z)=Ross(z)*atm.D(z)*atm.dR(z);

%% ... and continue to next zone
end

%% Update the atmosphere's opacity data
atm.opacity=Ross;
atm.opticalDepth=tau;