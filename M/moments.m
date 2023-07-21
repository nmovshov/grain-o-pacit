% MOMENTS   Some moments of the size distribution in steady state.
clear all
close all
clc
physunits off
[fileName pathName]=uigetfile();
load([pathName fileName])
nd=run.nd{end};
md=nd.*repmat(dust.massBin,length(atm.Z),1);
%% first moment - mean radius
for j=1:length(atm.Z)
    m1(j)=sum(dust.sizeBin.*nd(j,:))/sum(nd(j,:));
end
%% second moment - mean suface area
for j=1:length(atm.Z)
    m2(j)=sum((dust.sizeBin.^2).*nd(j,:))/sum(nd(j,:));
end
%% second moment corrected for scattering
for j=1:length(atm.Z)
    nupeak=5.879e10*si.Hz/si.K*atm.T(j); % Location of peak in Planck's function
    nr=1.5; ni=0.01;
    sigma=mie(nupeak,dust.sizeBin,nr,ni);
    q2(j)=sum(sigma.*nd(j,:))/sum(nd(j,:));
end
%% first moment - mean radius by mass
for j=1:length(atm.Z)
    M1(j)=sum(dust.sizeBin.*md(j,:))/sum(md(j,:));
end
%% second moment - mean suface area by mass
for j=1:length(atm.Z)
    M2(j)=sum((dust.sizeBin.^2).*md(j,:))/sum(md(j,:));
end
%% second moment corrected for scattering
for j=1:length(atm.Z)
    nupeak=5.879e10*si.Hz/si.K*atm.T(j); % Location of peak in Planck's function
    nr=1.5; ni=0.01;
    sigma=mie(nupeak,dust.sizeBin,nr,ni);
    Q2(j)=sum(sigma.*md(j,:))/sum(md(j,:));
end