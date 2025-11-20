function w = SincWindow(L, Z, sflag)
if nargin <= 1
  Z = 1;
end
if nargin <= 2
  sflag = 'symmetric';
end
L = double(L);
if strcmp(sflag, 'periodic')
  L = L + 1;
end
wp = ((0:L-1).'-(L-1)/2)./(L-1);
w = sinc(2.*Z.*wp);
if isnan(w(1))
  w = 1;
end
if strcmp(sflag, 'periodic')
  w = w(1:end-1);
end
if isscalar(w)
  w = 1;
end

end
