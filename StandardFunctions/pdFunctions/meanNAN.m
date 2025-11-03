function m=meanNAN(mi,dim,confidence)
% do a mean over a matrix even if it is including nan.

[s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10)]=size(mi);
if nargin==1;
    dim=find(s>1,1,'first');
    if isempty(dim); dim=1;end
end
if nargin<=2
    confidence=0;
end
nNAN=sum(isnan(mi),dim);
confidenceNAN=((s(dim)-nNAN)<=s(dim)*confidence);
mi(isnan(mi))=0;
m=sum(mi,dim)./(s(dim)-nNAN);
m(confidenceNAN)=nan;


%% ------------------------------------------------------------------------    
% (C) Copyright 2012 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%------------------------------------------------------------------------    
