% function res=singleSimEIModel(beta,sigma,thetaE,totSim)
%
% Call mex function eiModel to run a simulation of totSim seconds with
% parameters beta,sigma and thetaE.
%
% Example call for Fig 5C:
% singleSimEIModel(0.7,3.5,4.8,20)
%
% Created by Daniel Jercog


function res=singleSimEIModel(beta,sigma,thetaE,totSim)

run('./eiModelParams.m');
param.thetaE=param.Edesp-thetaE;
param.sigmaN=sigma;
param.betaE=beta;

paramS.t_start=0;
paramS.t_end=totSim*1000; %milliseconds
paramS.dt=0.2;
paramS.x0=0;
paramS.tauN=1;
paramS.sigmaN=param.sigmaN;

fParam=[param.Jee,param.Jei,param.Jie,param.Jii,param.thetaE,param.thetaI,param.Eslope,param.Edesp,param.Islope,param.Idesp,param.betaE,param.tauE,param.tauI,param.tauAdapt,paramS.tauN,paramS.sigmaN];
sParam=[paramS.t_start,paramS.t_end,paramS.dt];

%initial conditions for rateE, rateI, adaptation
iniC=[0,0,0];
res=eiModel(sParam,fParam,iniC);

figure;hold on;
plot(res(2,:));
plot(res(1,:),'r');
plot(res(3,:),'k');

