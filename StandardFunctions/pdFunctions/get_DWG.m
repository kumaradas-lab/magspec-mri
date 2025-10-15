function bNormEnd=get_DWG(Time,Amp,P)
if nargin<3; P=ones(size(Amp)); end

%%
% epsilon=[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5;
% Time=[ [-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5...
%         ,[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5+2...

% Time=[ [-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]+Seq.tEcho/2-Seq.DWG.DeltaTime/2 ...
%         ,[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]+Seq.tEcho/2+Seq.DWG.DeltaTime/2 ...
%         ,[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]+Seq.tEcho*3/2-Seq.DWG.DeltaTime/2 ...
%         ,[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]+Seq.tEcho*3/2+Seq.DWG.DeltaTime/2];
% Amp=[0,1,1,0,0,-1,-1,0,0,1,1,0,0,-1,-1,0];
% P=[0,1,1,0,0,-1,-1,0,0,-1,-1,0,0,1,1,0];

AmpTime=cumsum([0,diff(Time).*(Amp(1:end-1).*P(1:end-1)+Amp(2:end).*P(2:end))./2],2);
TimeN=unique([linspace(Time(1),Time(end),1000001),Time]);
AmpN=interp1(Time,Amp,TimeN,'linear');
PN=sign(interp1(Time,P,TimeN,'linear')).*sign(AmpN);
AmpN=abs(AmpN);
AmpTimeN=cumsum([0,diff(TimeN).*(AmpN(1:end-1).*PN(1:end-1)+AmpN(2:end).*PN(2:end))./2],2);
sigma=AmpTime(end);
AmpTime2N=cumsum([0,diff(TimeN).*((AmpTimeN(1:end-1).^2).*PN(1:end-1)+(AmpTimeN(2:end).^2).*PN(2:end))./2],2);
bNorm=cumsum([0,diff(TimeN).*(AmpTime2N(1:end-1)+AmpTime2N(2:end))./2],2);
bNormEnd=bNorm(end);

figure(11)
ax(1)=subplot(4,1,1);
plot(ax(1),Time,Amp,TimeN,AmpN.*PN)
legend('Amp', 'AmpN')
ax(2)=subplot(4,1,2);
plot(ax(2),TimeN,AmpTimeN)
legend('AmpTimeN')
ax(3)=subplot(4,1,3);
plot(ax(3),TimeN,AmpTime2N)
legend('AmpTime2N')
ax(4)=subplot(4,1,4);
plot(ax(4),TimeN,bNorm)
legend('bNorm')
linkaxes(ax,'x');
ylabel(['bNorm = ' num2str(bNorm(end))])

%%
% % epsilon=[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5;
% epsilon=[ [-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5...
%         ,[-Seq.DWG.Duration/2,-Seq.DWG.Duration/2+Seq.DWG.tRamp,Seq.DWG.Duration/2-Seq.DWG.tRamp,Seq.DWG.Duration/2]./Seq.DWG.Duration+0.5+2];
% s=[0,1,1,0,0,1,1,0];
% P=[0,1,1,0,0,-1,-1,0];
% 
% S=cumsum([0,diff(epsilon).*(s(1:end-1).*P(1:end-1)+s(2:end).*P(2:end))./2],2);
% epsilonN=unique([linspace(0,epsilon(end),100001),epsilon]);
% sN=interp1(epsilon,s,epsilonN,'linear');
% PN=sign(interp1(epsilon,P,epsilonN,'linear')).*sign(sN);
% sN=abs(sN);
% SN=cumsum([0,diff(epsilonN).*(sN(1:end-1).*PN(1:end-1)+sN(2:end).*PN(2:end))./2],2);
% sigma=S(end);
% lamda=1./sigma.*cumsum([0,diff(epsilonN).*(SN(1:end-1)+SN(2:end))./2],2);
% kappa=1./sigma.^2.*cumsum([0,diff(epsilonN).*((SN(1:end-1).^2).*PN(1:end-1)+(SN(2:end).^2).*PN(2:end))./2],2);
% S2N=cumsum([0,diff(epsilonN).*((SN(1:end-1).^2).*PN(1:end-1)+(SN(2:end).^2).*PN(2:end))./2],2);
% b=cumsum([0,diff(epsilonN).*(S2N(1:end-1)+S2N(2:end))./2],2);
% 
% figure(11)
% plot(epsilon,s,epsilonN,sN.*PN,epsilonN,SN,epsilonN,lamda,epsilonN,kappa,epsilonN,S2N,epsilonN,b)
% legend('s', 'sN', 'SN', 'lamda', 'kappa', 'S2N', 'b')