function w=RaisedCosine(L,sflag)
if nargin == 1
  sflag='symmetric';
end
L=double(L);

switch sflag
  case 'periodic'
    toAdd=zeros(numel(L),1);
  case 'symmetric'
    toAdd=0.5+zeros(numel(L),1);
  case 'fftGrid'
    toAdd=0.5*mod(L(:),2);
  otherwise
    error('second argument: periodic, symmetric or fftGrid')
end

r2=zeros([L(:).',1]);
for t=1:numel(L)
  if L(t)>1
    r2=bsxfun(@plus,r2,shiftdim( (toAdd(t)+(0:L(t)-1).') /L(t) - 0.5, -(t-1) ).^2);
  end
end
r = r2.^0.5;
r(r>0.5) = 0.5;
w = 0.5+0.5*cosd(360.*r);
end
