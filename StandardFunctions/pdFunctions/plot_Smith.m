function [hax, Network] = plot_Smith(HW, Network, hax, type)
%% Plot a Smith chart with the data from sequence_Network
%
%     [hax, Network] = plot_Smith(HW, Network, hax, type)
%
% Plot the measured admittance and impedance in a Smith chart. Additionally,
% calculate the loaded Q-factor of the coil with the method described in
% https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf
%
% INPUT:
%
%   HW
%       HW structure or PD.HWClass object.
%
%   Network
%       Structure with at least the following fields:
%     Frequency
%         Vector with frequencies in Hz.
%     Reflection
%         Complex reflection ratio.
%
%   hax
%       Optional. Handle to an axes to use for the Smith chart. If omitted or
%       empty, a new figure (201) is created for the chart.
%
%   type
%       Optional. Type of Smith chart ('z', 'y', 'yz'), default: 'z'.
%
%
% OUTPUT:
%
%   hax
%       Handle to the axes containing the Smith chart.
%
%   Network
%       Structure as in input with additional fields (see get_QFactor_Smith).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Check input parameters
if nargin < 4
  type = 'Z';
end
if nargin < 3 || isempty(hax)
  fh = figure(201);
  % scrsz = get(0,'ScreenSize');
  % set(fh, 'Position', [100 scrsz(4)-1300 1200 1200])
  set(fh, 'Name', 'Smith Chart', 'NumberTitle', 'off');
  hax = gca;
end


%% Create Smith chart
Z0 = HW.RX.Rin;
smithValue = PD.smith.createGrid({0.1:0.1:0.4, 0.2:0.2:0.8, 0.4:0.4:2.0, 1:1:5, 5:5:10, 10:10:20 });
smithMinorValue = PD.smith.createGrid({0.025:0.025:0.4, 0.05:0.05:0.8, 0.1:0.1:2.0, 0.2:0.2:5, 1:1:10, 2:2:20, 10:10:50 });
co = get(hax, 'ColorOrder');
impColor = co(1,:);
color = 1-0.9*(1-impColor);
minorColor = 1-0.4*(1-impColor);
admColor = [1 0 1];
subColor = 1-0.9*(1-admColor);
subMinorColor = 1-0.4*(1-admColor);
if strcmpi(type, 'Y')
  % swap colors
  c = color;  color = subColor;  subColor = c;
  c = minorColor;  minorColor = subMinorColor;  subMinorColor = c;
end

smithProps = struct('Type', type, 'Value', smithValue, 'Color', color, ...
  'MinorValue', smithMinorValue, 'MinorLineStyle', '-', ...
  'MinorColor', minorColor, 'MinorVisible', 'on', ...
  'SubLineStyle', '-', 'SubColor', subColor, ...
  'SubMinorLineStyle', '-', 'SubMinorColor', subMinorColor, ...
  'SubMinorVisible', 'on');

hsm = getappdata(hax, 'smithObject');
if isempty(hsm) || ~ishghandle(hsm, 'smith')
  cla(hax)
  hsm = PD.smith(hax, smithProps);
else
  set(hsm, smithProps);
  hsm.clear();
  hold(hax, 'on');
end
setappdata(hax, 'smithObject', hsm);


%% plot (interpolated) impedance and admittance
Frequency = interpn(Network.Frequency, 4);
Reflection = interpn(Network.Reflection, 4, 'spline');
[MinReflection, iMinReflection] = min(Reflection);
fMinReflection = Frequency(iMinReflection);
fMinFormatStr = '%.1f dB @ %.3f MHz';
fMinStr = sprintf(fMinFormatStr, 20*log10(abs(MinReflection)), fMinReflection/1e6);

smithline = plot(hax, Reflection, 'LineWidth', 2);
smithlineY = plot(hax, -Reflection, 'm', 'LineWidth', 2);
hl = legend([smithline,smithlineY], {'impedance', 'admittance'});
if ~verLessThan('Matlab', '9.2')
  set(hl, 'AutoUpdate', 'off');  % Do not automatically add new entries to legend
end

% smithl.XData=get(smithline,'XData');
% smithl.YData=get(smithline,'YData');
% smithl.ZData=get(smithline,'ZData');

xSmith = get(smithline, 'XData');
ySmith = get(smithline, 'YData');

xMinReflection = double(interp1(Frequency, xSmith, fMinReflection, 'nearest'));
yMinReflection = double(interp1(Frequency, ySmith, fMinReflection, 'nearest'));


%% Loaded Q-factor (https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf):
% "Draw line" from (0,0) trough fMinReflection
% Find frequencies (f1,f2) of points on S11 (Q-circle) at 45 degrees from the
% "outer edge" of the Q-circle (i.e. furthest from that line)
% --> Q = fMinReflection / |f1 - f2|
% That direct approach is rather unstable for noisy data. Instead using the
% following equivalent approach:

% Find circle that fits points around min. reflection.
Network.minRefFrequencySSInt  = NaN;
Network.minRefReflection_dB   = NaN;
Network.QLFrequencySSInt      = NaN;

Network = get_loadedQFactor_Smith(Network);

if isfinite(Network.minRefFrequencySSInt)
  plot(hax, Network.Reflection([Network.Q_iStart, Network.Q_iStop]), 'kd')
  plot(hax, -Network.Reflection([Network.Q_iStart, Network.Q_iStop]), 'kd')
  xB = real(Network.Q_circle_center);
  yB = imag(Network.Q_circle_center);
  rB = Network.Q_circle_radius;
  fplot(hax, @(t) rB*sin(t)+xB, @(t) rB*cos(t)+yB, [-pi, pi], ':k', 'MeshDensity', 360*1);
  xB = -xB;
  yB = -yB;
  fplot(hax, @(t) rB*sin(t)+xB, @(t) rB*cos(t)+yB, [-pi, pi], ':k', 'MeshDensity', 360*1);

  plot(hax, Network.Q_circle_minCenter, 'ko')
  centerCircle2minRefCircle = Network.Q_circle_center - Network.Q_circle_minCenter;
  line2plot = [Network.Q_circle_minCenter; Network.Q_circle_center; NaN; ...
    Network.Q_circle_center+centerCircle2minRefCircle.*exp(-1i*(+pi/2)); ...
    Network.Q_circle_center+centerCircle2minRefCircle.*exp(-1i*(-pi/2)); ...
    NaN];
  plot(hax, [line2plot; -line2plot], ':ko');  % plot diameter and line to minRefCycle

  fMinFormatStr = [fMinFormatStr ' (Q_L = %.1f)'];
  fMinStr = sprintf(fMinFormatStr, ...
    20*log10(abs(Network.Q_circle_minCenter)), Network.minRefFrequencySSInt/1e6, ...
    Network.QLFrequencySSInt);
  xMinReflection = real(Network.Q_circle_minCenter);
  yMinReflection = imag(Network.Q_circle_minCenter);
end


%% Mark at Larmor frequency and resonant frequency
x = interp1(Frequency, xSmith, HW.fLarmor, 'nearest');
y = interp1(Frequency, ySmith, HW.fLarmor, 'nearest');
xy = x+1i*y;
% xy=Reflection;
Zt = (1+xy)./(1-xy);
text(double(interp1(Frequency, xSmith, HW.fLarmor, 'nearest')), ...
  double(interp1(Frequency, ySmith, HW.fLarmor, 'nearest')), ...
  ['\color[rgb]{0 .5 0}\leftarrow fL ' num2str(Zt,'%10.2f') ], ...
  'HorizontalAlignment', 'left',...
  'FontWeight', 'bold', ...
  'Rotation', 90, 'FontSize', 12, 'Parent', hax);

text(xMinReflection, ...
  yMinReflection, ...
  ['\color[rgb]{0 .5 0}\leftarrow' fMinStr], ...
  'HorizontalAlignment', 'left', ...
  'FontWeight', 'bold', ...
  'Rotation', -90, 'FontSize', 12, 'Parent', hax);


%% Calculate parameters of series LC circuit
R = real(Zt)*Z0;
if imag(Zt)>=0
  X = imag(Zt)*Z0/(2*pi*HW.fLarmor);
  Xtext = 'L';
  Xeinheit = 'nH';
else
  X = -1/(imag(Zt)*Z0*(2*pi*HW.fLarmor));
  Xtext = 'C';
  Xeinheit = 'nF';
end
ZT = Zt*Z0;

title1 = sprintf('Series circuit X = %.1f Ohm, R = %.1f Ohm and %s = %.3f %s @ fLarmor (%.6f MHz)', ...
  ZT, R, Xtext, X*1e9, Xeinheit, HW.fLarmor/1e6);


%% Calculate parameters of parallel LC circuit
Y0 = 1/Z0;
x = interp1(Frequency, get(smithlineY, 'XData'), HW.fLarmor, 'nearest');
y = interp1(Frequency, get(smithlineY, 'YData'), HW.fLarmor, 'nearest');
xy = x + 1i*y;
% xy=Reflection;
Yt = (1+xy)./(1-xy);
Y = real(Yt)*Y0;
if imag(Yt) <= 0
  X = -1/(imag(Yt)*Y0*(2*pi*HW.fLarmor));
  Xtext = 'L';
  Xeinheit = 'nH';
  XaddText = 'C';
  XaddEinheit = 'nF';
else
  X = imag(Yt)*Y0 / (2*pi*HW.fLarmor);
  Xtext = 'C';
  Xeinheit = 'nF';
  XaddText = 'L';
  XaddEinheit = 'nH';
end
YT = Yt*Y0;
Xadd = 1/(X*4*pi^2*HW.fLarmor^2);

title2 = sprintf('Parallel circuit Y = %.3f 1/Ohm, R = %.1f Ohm and %s = %.3f %s @ fLarmor (%.6f MHz)', ...
  YT, 1/Y, Xtext, X*1e9, Xeinheit, HW.fLarmor/1e6);
title3 = sprintf('Parallel resonant circuit %s = %.3f %s; add %s = %.3f %s to get resonant @ fLarmor (%.6f MHz)', ...
  Xtext, X*1e9, Xeinheit, XaddText, Xadd*1e9, XaddEinheit, HW.fLarmor/1e6);


%% Mark at Larmor frequency and resonant frequency
text(double(interp1(Frequency, get(smithlineY, 'XData'), HW.fLarmor, 'nearest')), ...
  double(interp1(Frequency, get(smithlineY, 'YData'), HW.fLarmor, 'nearest')), ...
  ['\color[rgb]{0 .5 0}\leftarrow fL ' num2str(Yt, '%10.2f') ], ...
  'HorizontalAlignment', 'left',...
  'FontWeight', 'bold', ...
  'Rotation', 0, 'FontSize', 12, 'Parent', hax)

text(-xMinReflection, ...
  -yMinReflection, ...
  ['\color[rgb]{0 .5 0}' fMinStr ' \rightarrow'], ...
  'HorizontalAlignment', 'right', ...
  'FontWeight', 'bold', ...
  'Rotation', 0, 'FontSize', 12, 'Parent', hax)


%% Add title
set(hax, 'Unit', 'pixels')
if all(get(hax, 'Position') > [-inf,-inf,400,400])
  title(hax, {title1; title2; title3});
else
  title(hax, '');
end


% l1=plot(ReflectionRaw,'r');
% set(l1,'ZData',Frequency)
%
% l2=plot(HW.NetworkCal.Leerlauf.ReflectionRaw,'g');
% set(l2,'ZData',HW.NetworkCal.Leerlauf.Frequency)
%
% l3=plot(HW.NetworkCal.Kurzschluss.ReflectionRaw,'y');
% set(l3,'ZData',HW.NetworkCal.Kurzschluss.Frequency)
%
% l4=plot(HW.NetworkCal.Abschluss.ReflectionRaw,'m');
% set(l4,'ZData',HW.NetworkCal.Abschluss.Frequency)


hold(hax, 'off');

end
