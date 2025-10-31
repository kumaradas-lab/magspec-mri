function [sigma,lamda,kappa]=get_DWGslk(epsilon,s,doPlot)
if nargin<3; doPlot=0;end;
%%
% epsilon=[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5;
% s=[0,1,1,0];

S=cumsum([0,diff(epsilon).*(s(1:end-1)+s(2:end))./2],2);
epsilonN=unique([linspace(epsilon(1),epsilon(end),1000001),epsilon]);
sN=interp1(epsilon,s,epsilonN,'linear');
sN=abs(sN);
SN=cumsum([0,diff(epsilonN).*(sN(1:end-1)+sN(2:end))./2],2);
sigma=S(end);
lamda=1./sigma.*cumsum([0,diff(epsilonN).*(SN(1:end-1)+SN(2:end))./2],2);
kappa=1./sigma.^2.*cumsum([0,diff(epsilonN).*((SN(1:end-1).^2)+(SN(2:end).^2))./2],2);
S2N=cumsum([0,diff(epsilonN).*(SN(1:end-1).^2+SN(2:end).^2)./2],2);
b=cumsum([0,diff(epsilonN).*(S2N(1:end-1)+S2N(2:end))./2],2);

if doPlot
    figure(11)
    clf
    plot(epsilon,s,epsilonN,sN,epsilonN,SN,epsilonN,lamda,epsilonN,kappa,epsilonN,S2N,epsilonN,b)
    legend('s', 'sN', 'SN', 'lamda', 'kappa', 'S2N', 'b')
end

lamda=lamda(end);
kappa=kappa(end);

