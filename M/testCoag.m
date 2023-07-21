% testCoag  A verification of our coagulation algorythm

%% A. Constant kernel and all grains initially of equal size
% The soluion in this case is described in [1]. (This analytical solution
% applies to mass bins that are integer multiples of some m0.)
%
% [1] Comparison of Analytical and Physical Modeling of Planetesimal
% Accumulation, G. W. Whetherill, Icarus 88, 336-354, 1990.

%% Clear and set a workspace
clear all
close all
clc
global si params atm dust
physunits off
si=setUnits;

%% Run the model with a constant kernel
params=oParams('nBins',20,'maxTimeStep',1e6*si.s);
oSetup('atmos','atm2321.dat','index','ind.dat','source','tes2321.dat',...
    'kernel','fixed','massbin','geometric')
tt=(0:100)*1e8*si.s;
nd{1}=dust.nDensity(1,:);
h=waitbar(0);
for k=2:length(tt)
    jack(tt(k-1),tt(k),'noflux','nosource')
    nd{k}=dust.nDensity(1,:);
    waitbar(k/(length(tt)))
end
close(h)
pause(1)
save conker
%% Find total number of grains, from model and analyticaly
ndTotal=zeros(1,length(tt));
for k=1:length(nd)
    ndTotal(k)=sum(nd{k});
end
nd0=ndTotal(1);
alfa=dust.kernel(1);
f=@(t)(1+0.5*alfa*nd0.*t).^(-1);
ndTotalA=nd0*f(tt);

%% Find the distribution analyticaly
for j=1:length(nd)
    for k=1:params.nBins
        nk(k)=nd0*f(tt(j))^2*(1-f(tt(j)))^(k-1);
    end
    ndA{j}=nk;
end

%% Compare the model and analytic solutions
subplot(1,2,1)
plot(tt,ndTotal,'r+',tt,ndTotalA,'b-')
subplot(1,2,2)
semilogy(dust.sizeBin,nd{end},'r+',dust.sizeBin,ndA{end},'b-')
bla=mean(abs(ndTotal-ndTotalA)./ndTotal);

%% B. Linear kernel and all grains initially of equal size
% This case is also described in [1]. (This analytical solution applies
% to mass bins that are integer multiples of some m0.)

% %% Run the model with a linear kernel
% params=oParams('nBins',20,'maxTimeStep',1e6*si.s);
% oSetup('atmos','atm2321.dat','index','ind.dat','source','tes2321.dat',...
%     'kernel','linear','massbin','linear')
% tt=(0:100)*1e8*si.s;
% nd{1}=dust.nDensity(1,:);
% h=waitbar(0);
% for k=2:length(tt)
%     jack(tt(k-1),tt(k),'noflux','nosource')
%     nd{k}=dust.nDensity(1,:);
%     waitbar(k/(length(tt)))
% end
% close(h)
% pause(1)
% save linker
% %% Find total number of grains, from model and analyticaly
% ndTotal=zeros(1,length(tt));
% for k=1:length(nd)
%     ndTotal(k)=sum(nd{k});
% end
% nd0=ndTotal(1);
% beta=dust.kernel(1)/2;
% f=@(t)exp(-beta*nd0*t);
% ndTotalA=nd0*f(tt);
% 
% %% Find the distribution analyticaly
% for j=1:length(nd)
%     for k=1:params.nBins
%         nk(k)=nd0*(k^(k-1))/(factorial(k))*f(tt(j))*(1-f(tt(j)))^(k-1)*...
%             exp(-k*(1-f(tt(j))));
%     end
%     ndA{j}=nk;
% end
% 
% %% Compare
% subplot(1,2,1)
% plot(tt,ndTotal,'r+',tt,ndTotalA,'b-')
% subplot(1,2,2)
% semilogy(dust.sizeBin,nd{end},'r+',dust.sizeBin,ndA{end},'b-')