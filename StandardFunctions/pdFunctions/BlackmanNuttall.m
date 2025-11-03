function w = BlackmanNuttall(L, sflag)
%% Create Blackman-Nuttall filter window
%Das Blackman-Nuttall-Fenster ist bis auf die vier fast identischen Koeffizienten identisch mit dem Blackman-Harris-Fenster, was den Einfluss der notwendigen Genauigkeit bei der Implementierung der Koeffizienten bei dieser Klasse von Fensterfunktionen verdeutlicht.
%   w = BlackmanNuttall(L, sflag)
%
% INPUT:
%   L       number of samples
%   sflag   'symmetric' or 'periodic' (default: 'symmetric')
%
%
% ------------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

if nargin == 1
    sflag = 'symmetric';
end   

L = double(L);
if strcmp(sflag, 'periodic'), L = L+1; end;
wp = (0:L-1).';
a0=0.3635819;
a1=0.4891775;
a2=0.1365995;
a3=0.0106411;
w=  + a0 ...
    - a1 .* cos(2.*pi.*wp./(L-1))...
    + a2 .* cos(4.*pi.*wp./(L-1))...
    - a3 .* cos(6.*pi.*wp./(L-1));
    
  
if strcmp(sflag, 'periodic'), w = w(1:end-1); end;
if numel(w) == 1, w = 1; end

end