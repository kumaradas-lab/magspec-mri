function  window=LorentzToGauss(windowSize,LB,Wmax,fSample,doPlot)
%%
if nargin<1; windowSize=[]; end
if nargin<2; LB=[]; end
if nargin<3; Wmax=[]; end
if nargin<4; fSample=[]; end
if nargin<5; doPlot=0; end
if isempty(windowSize); windowSize=2^10;end;
if isempty(fSample); fSample=windowSize;end;
if isempty(LB); LB=-(windowSize/fSample)*2;end;
if isempty(Wmax); Wmax=1/3;end;
alfa=-LB*pi/2/Wmax/(windowSize/fSample);
k=(0:windowSize-1).';
window=exp(-pi*LB*k./fSample)  .*  exp(-alfa.*(k./fSample).^2);%
if(doPlot);
    figure(1); 
    plot(k./fSample,window); 
    xlabel('time / s')
    ylabel('window amplitude')
end