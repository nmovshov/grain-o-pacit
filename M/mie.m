function sigma=mie(nu,a,nr,ni)
%MIE    Cuzzie's prescription.

si=setUnits;
c=si.speed_of_light;
lambda=c/nu;
x=2*pi*a/lambda;

if ni<3
    cosbar=0.2*ones(size(x));
    cosbar(x>2.5)=0.8;
else
    cosbar=-0.2*ones(size(x));
    cosbar(x>2.5)=0.5;
end
Qs=ones(size(x));
Qs(x<1.3)=8/3*x(x<1.3).^4.*((nr^2-1)/(nr^2+2)).^2;
Qs(x>=1.3)=2*x(x>=1.3).^2.*((nr-1)^2+ni^2);
Qs(Qs>1)=1;
Qa=12*x.*(2*nr.*ni)./((nr.^2+ni.^2+2).^2+...
    4*nr.^2.*ni.^2);
Qa(Qa>1)=1;
Qe=Qs+Qa;
Q=Qe.*(1-(Qs./Qe).*cosbar);
sigma=Q.*pi.*a.^2;