function ComplexDataAnglewMean=get_MeanPhaseWeighted(ComplexData,dim)
if nargin==1; dim=find(size(ComplexData)>1,1,'first');end
    ComplexDataw=abs(ComplexData);
    ComplexDataAngle=unwrap(angle(ComplexData),[],dim);
    ComplexDataAnglew=ComplexDataAngle.*ComplexDataw;
    ComplexDataAnglewMean=sum(ComplexDataAnglew,dim)./sum(ComplexDataw,dim);
end
