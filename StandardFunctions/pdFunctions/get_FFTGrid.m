function FFTGrid=get_FFTGrid(span,n)
% FFTGrid=get_FFTGrid(span,n)
% span        60  60  60  60  60  60
% n            0   1   2   3   4   5
% FFTGrid(1)  []   0 -30 -20 -30 -24
% FFTGrid(2)           0   0 -15 -12
% FFTGrid(3)              20   0   0
% FFTGrid(4)                  15  12
% FFTGrid(5)                      24
%%

%  span=3;
%  n=4;
%  span=permute(cat(3,[20:22;30:32],[20:22;30:32]),[4,1,2,3]);
%  n=permute(cat(3,[2:3:9;0:2],[2:3:9;0:2]+5),[4,1,2,3]);
if isscalar(span) && isscalar(n)
    FFTGrid=((0:1:(n-1))-(floor(n/2))).'*span/n;
else
    n1=max(n(:));
    s0=(floor(n./2));
    FFTGrid=bsxfun(@times,bsxfun(@minus,(0:1:(n1-1)).',s0),span./n);
    FFTGrid(isnan(FFTGrid)|isinf(FFTGrid))=0;
end

%%
% x = [((0:span/n:(span-span/n))-span/n*(floor(n/2))).',...
%       (-span/n*(floor(n/2)):span/n:(span-span/n)-span/n*(floor(n/2))).',...
%       ((0:1/n:(1-1/n))*span-span/n*(floor(n/2))).',...
%       ((0:1/n:(1-1/n))-1/n*(floor(n/2))).'*span,...
%       ((0:1:(n-1))-(floor(n/2))).'*span/n]


