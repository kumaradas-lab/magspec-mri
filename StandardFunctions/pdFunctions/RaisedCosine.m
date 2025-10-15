function w=RaisedCosine(L,sflag)
if nargin == 1
    sflag='symmetric';
end   
L=double(L);
if strcmp(sflag,'periodic');L=L+1;end;
wp=(0:L-1).';
w=  +0.5...
    -0.5 .*cosd(360.*wp./(L-1));
if isnan(w(1));w(1)=1;end;
if strcmp(sflag,'periodic');w=w(1:end-1);end;
end
