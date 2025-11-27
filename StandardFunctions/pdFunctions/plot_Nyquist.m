function hf = plot_Nyquist(HW, Network, hax)
%% Plot Nyquist diagram of network
%
%     hf = plot_Nyquist(HW, Network, hax)
%
% INPUT:
%   HW
%           HW object or structure.
%   Network
%           Network structure as returned by sequence_Network.
%   hax
%           Optional handle to axes that is used for plotting the Nyquist
%           impedance diagram. If omitted, a figure is created that displays two
%           axes with the Nyquist impedance and admittance diagrams.
%
% OUTPUT:
%   hf
%           Only if called with two input argument or if hax is empty, the
%           handle to the figure that contains the Nyquist diagrams.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input
if nargin == 2, hax = []; end

%% create figure
if isempty(hax)
  hf = figure(122); clf(hf);
  hax = subplot(2,1,1, 'Parent', hf);
end

%% plot impedance diagram
Frequency = interpn(Network.Frequency, 4);
Reflection = interpn(Network.Reflection, 4, 'spline');
[MinReflection, iMinReflection] = min(Reflection);
fMinReflection = Frequency(iMinReflection);

Z0 = HW.RX.Rin;
z = (1+Reflection) ./ (1-Reflection);
Z = z*Z0;
lR = plot(hax, Z);
title(hax, 'Impedance');
set(lR, 'ZData', Frequency);
xlabel(hax, ['R in ' 937]);  % Ohm
ylabel(hax, ['X in ' 937]);  % Ohm
xlim(hax, [0, 500]);
ylim(hax, [-250, 250]);
axis(hax, 'square');
grid(hax, 'on');

% add labels to plot
fLarmorR = double(interp1(Frequency, real(Z), HW.fLarmor, 'nearest'));
fLarmorX = double(interp1(Frequency, imag(Z), HW.fLarmor, 'nearest'));
minReflR = double(interp1(Frequency, real(Z), fMinReflection, 'nearest'));
minReflX = double(interp1(Frequency, imag(Z), fMinReflection, 'nearest'));

if fLarmorX < minReflX
  arrowFLarmor = '\uparrow';
  alignmentFLarmor = 'top';
  arrowMin = '\downarrow';
  alignmentMin = 'bottom';
else
  arrowFLarmor = '\downarrow';
  alignmentFLarmor = 'bottom';
  arrowMin = '\uparrow';
  alignmentMin = 'top';
end

% Larmor frequency
hArrowFLarmor = text(fLarmorR, fLarmorX, ...
  arrowFLarmor, ...
  'FontSize', 14, ...
  'FontWeight', 'bold', ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', alignmentFLarmor, ...
  'Rotation', 0, 'Parent', hax);
arrowFLarmorExtent = get(hArrowFLarmor, 'Extent');
if strcmp(arrowFLarmor, '\uparrow')
  % move uparrow slightly upwards to touch the line
  set(hArrowFLarmor, 'Position', [fLarmorR, fLarmorX + 0.2*arrowFLarmorExtent(4)]);
end
text(fLarmorR + 0.7*arrowFLarmorExtent(3), ...
  fLarmorX + strcmp(arrowFLarmor, '\downarrow')*0.2*arrowFLarmorExtent(4), ...
  sprintf('%.1f dB @ f_L', 20*log10(abs(interp1(Frequency, Reflection, HW.fLarmor, 'nearest')))), ...
  'FontWeight', 'bold', ...
  'HorizontalAlignment', 'left', ...
  'VerticalAlignment', alignmentFLarmor, ...
  'Rotation', 0, 'Parent', hax, ...
  'Color', [0.106, 0.31, 0.208]);

% minimum reflection
hArrowMin = text(minReflR, minReflX, ...
  arrowMin, ...
  'FontSize', 14, ...
  'FontWeight', 'bold', ...
  'Margin', 1, ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', alignmentMin, ...
  'Rotation', 0, 'Parent', hax);
arrowMinExtent = get(hArrowMin, 'Extent');
if strcmp(arrowMin, '\uparrow')
  % move uparrow slightly upwards to touch the line
  set(hArrowMin, 'Position', [fLarmorR, fLarmorX + 0.2*arrowMinExtent(4)]);
end
text(minReflR + 0.7*arrowMinExtent(3), ...
  minReflX + strcmp(arrowMin, '\downarrow')*0.2*arrowMinExtent(4), ...
  sprintf('%.1f dB @ %.3f MHz', 20*log10(abs(MinReflection)), fMinReflection/1e6), ...
  'FontWeight', 'bold', ...
  'Margin', 1, ...
  'HorizontalAlignment', 'left', ...
  'VerticalAlignment', alignmentMin, ...
  'Rotation', 0, 'Parent', hax);


% hold(hax, 'on');
% plot(hax, double(interp1(Frequency, real(Z), HW.fLarmor, 'nearest')), ...
%   double(interp1(Frequency, imag(Z), HW.fLarmor, 'nearest')), '-xr');
% plot(hax, double(interp1(Frequency, real(Z), fMinReflection, 'nearest')), ...
%   double(interp1(Frequency, imag(Z), fMinReflection, 'nearest')), '-xb');
% hold(hax, 'off');

% zRaw=(1+ReflectionRaw)./(1-ReflectionRaw);
% ZRaw=zRaw*Z0;
% lRaw=plot(ZRaw,'r');
% set(lRaw,'ZData',Frequency)
%
% zRaw=(1+HW.NetworkCal.Leerlauf.ReflectionRaw)./(1-HW.NetworkCal.Leerlauf.ReflectionRaw);
% ZRaw=zRaw*Z0;
% lLeer=plot(ZRaw,'g');
% set(lLeer,'ZData',HW.NetworkCal.Leerlauf.Frequency)
%
% zRaw=(1+HW.NetworkCal.Kurzschluss.ReflectionRaw)./(1-HW.NetworkCal.Kurzschluss.ReflectionRaw);
% ZRaw=zRaw*Z0;
% lKS=plot(ZRaw,'y');
% set(lKS,'ZData',HW.NetworkCal.Kurzschluss.Frequency)
%
% zRaw=(1+HW.NetworkCal.Abschluss.ReflectionRaw)./(1-HW.NetworkCal.Abschluss.ReflectionRaw);
% ZRaw=zRaw*Z0;
% lAS=plot(ZRaw,'m');
% set(lAS,'ZData',HW.NetworkCal.Abschluss.Frequency)


%% plot admittance diagram
if exist('hf', 'var')
  hax2 = subplot(2,1,2, 'Parent', hf);
  Y = 1./(z.*Z0);
  lY = plot(hax2, Y);
  title(hax2, 'Nyquist Admittance');
  set(lY, 'ZData', Frequency);
  xlabel(hax2, '1/R');
  ylabel(hax2, '1/X');
  xlim(hax2, [0, 0.5])
  ylim(hax2, [-0.25, 0.25])
  axis(hax2, 'square');
  grid(hax2, 'on');

  text(double(interp1(Frequency, get(lY, 'XData'), HW.fLarmor, 'nearest')), ...
    double(interp1(Frequency, get(lY, 'YData'), HW.fLarmor, 'nearest')), ...
    'fL \downarrow', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', ...
    'Rotation', 0, 'Parent', hax2);

  text(interp1(Frequency, get(lY, 'XData'), fMinReflection, 'nearest'), ...
    interp1(Frequency, get(lY, 'YData'), fMinReflection, 'nearest'), ...
    ['\uparrow ' num2str(20*log10(abs(MinReflection)), '%.1f') ' dB @ ' num2str(fMinReflection/1e6, '%.3f') ' MHz'], ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'Rotation', 0, 'Parent', hax2);
elseif nargout > 0
  hf = [];
end

% hold off

end
