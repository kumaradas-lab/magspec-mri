function hax = plot_Smith(HW, Network, hax, type)
%% Plot a Smith chart with the data from sequence_Network
%
%     hax = plot_Smith(HW, Network, hax, type)
%
% Plot the measured admittance and impedance in a Smith chart. Additionally,
% calculate the loaded Q-factor of the coil with the method described in
% https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf
%
% INPUT:
%   HW        HW structure or PD.HW object
%   Network   Structure with at least the following fields:
%     Frequency     Vector with frequencies in Hz
%     Reflection    Complex reflection ratio
%   hax       Optional. Handle to an axes to use for the Smith chart. If omitted
%             or empty, a new figure (201) is created for the chart.
%   type      Optional. Type of Smith chart ('z', 'y', 'yz'), default: 'z'.
%
% OUTPUT:
%   hax       Handle to the axes containing the Smith chart
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2019 Pure Devices GmbH, Wuerzburg, Germany
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
  hold(hax, 'all')
end
setappdata(hax, 'smithObject', hsm);

%% plot (interpolated) impedance and admittance
Frequency = interpn(Network.Frequency, 4);
Reflection = interpn(Network.Reflection, 4, 'spline');
[MinReflection, iMinReflection] = min(Reflection);
fMinReflection = Frequency(iMinReflection);

smithline = plot(hax, Reflection, 'LineWidth', 2);
smithlineY = plot(hax, -Reflection, 'm', 'LineWidth', 2);
legend(hax, 'impedance', 'admittance')

% smithl.XData=get(smithline,'XData');
% smithl.YData=get(smithline,'YData');
% smithl.ZData=get(smithline,'ZData');

xSmith = get(smithline, 'XData');
ySmith = get(smithline, 'YData');

%% Loaded Q-factor (https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf):
% "Draw line" from (0,0) trough fMinReflection
% Find frequencies (f1,f2) of points on S11 (Q-circle) at 45 degrees from the
% "outer edge" of the Q-circle (i.e. furthest from that line)
% --> Q = fMinReflection / |f1 - f2|
% That direct approach is rather unstable for noisy data. Instead using the
% following equivalent approach:

% Find circle that fits points around min. reflection.
tempRef = sort(abs(Reflection));
iStart = find(abs(Reflection) < tempRef(round(numel(Reflection)/10)), 1, 'first');
iStop  = find(abs(Reflection) < tempRef(round(numel(Reflection)/10)), 1, 'last');
Par = CircleFitByPratt([xSmith(iStart:iStop)', ySmith(iStart:iStop)']);
xCenterQ = Par(1);
yCenterQ = Par(2);

% find points with minimal distance to orthogonal diameter (diameter line that
% is orthogonal to the line from the center of the smith chart through the
% center of the Q-circle)
reflCenterQ = sqrt(xCenterQ^2 + yCenterQ^2);
distanceOrth = abs(xCenterQ*xSmith + yCenterQ*ySmith - reflCenterQ^2) / reflCenterQ;
minDist1 = min(distanceOrth(Frequency>fMinReflection));
if ~isempty(minDist1)
  if1 = find(distanceOrth == minDist1, 1, 'first');
end
minDist2 = min(distanceOrth(Frequency<fMinReflection));
if ~isempty(minDist2)
  if2 = find(distanceOrth == minDist2, 1, 'last');
end

% calculate loaded Q factor and construct string for plot
fMinFormatStr = '%.1f dB @ %.3f MHz';
if isempty(minDist1) || isempty(minDist2) || (if1 - if2 < 5) || ...
    (Par(3) > 0.9) || ((minDist1 > 0.1*Par(3)) && (minDist2 > 0.1*Par(3)))
  fMinStr = sprintf(fMinFormatStr, 20*log10(abs(MinReflection)), fMinReflection/1e6);
else
  if minDist1 > 0.1*Par(3)
    df = 2*abs(Frequency(if2) - fMinReflection);
    if1 = [];
  elseif minDist2 > 0.1*Par(3)
    df = 2*abs(Frequency(if1) - fMinReflection);
    if2 = [];
  else
    df = abs(Frequency(if1) - Frequency(if2));
  end
  QL = fMinReflection / df;
  fMinFormatStr = [fMinFormatStr ' (Q_L = %.1f)'];
  fMinStr = sprintf(fMinFormatStr, 20*log10(abs(MinReflection)), fMinReflection/1e6, QL);
  plot(hax, Reflection([if1, if2]), 'ro')
  plot(hax, -Reflection([if1, if2]), 'ro')
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

text(double(interp1(Frequency, xSmith, fMinReflection, 'nearest')), ...
  double(interp1(Frequency, ySmith, fMinReflection, 'nearest')), ...
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

text(double(interp1(Frequency, get(smithlineY, 'XData'), fMinReflection, 'nearest')), ...
  double(interp1(Frequency, get(smithlineY, 'YData'), fMinReflection, 'nearest')), ...
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

function Par = CircleFitByPratt(XY)
%% Circle fit by Pratt
%   V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%   Computer Graphics, Vol. 21, pages 145-152 (1987)
%
% Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
% Output: Par = [a b R] is the fitting circle:
%                       center (a,b) and radius R
%
% Note: This fit does not use built-in matrix functions (except "mean"),
%       so it can be easily programmed in any programming language

n = size(XY,1);      % number of data points

centroid = mean(XY);   % the centroid of the data set

% computing moments (note: all moments will be normed, i.e. divided by n)

Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;

for i=1:n
  Xi = XY(i,1) - centroid(1);  %  centering data
  Yi = XY(i,2) - centroid(2);  %  centering data
  Zi = Xi*Xi + Yi*Yi;
  Mxy = Mxy + Xi*Yi;
  Mxx = Mxx + Xi*Xi;
  Myy = Myy + Yi*Yi;
  Mxz = Mxz + Xi*Zi;
  Myz = Myz + Yi*Zi;
  Mzz = Mzz + Zi*Zi;
end

Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;

% computing the coefficients of the characteristic polynomial

Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Mxz2 = Mxz*Mxz;
Myz2 = Myz*Myz;

A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;

epsilon = sqrt(eps(class(XY)));
ynew=1e+20;
IterMax=20;
xnew = 0;

% Newton's method starting at x=0

for iter=1:IterMax
  yold = ynew;
  ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
  if (abs(ynew)>abs(yold))
    disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
    xnew = 0;
    break;
  end
  Dy = A1 + xnew*(A22 + 16*xnew*xnew);
  xold = xnew;
  xnew = xold - ynew/Dy;
  if (abs((xnew-xold)/xnew) < epsilon), break, end
  if (iter >= IterMax)
    disp('Newton-Pratt did not converge');
    xnew = 0;
  end
  if (xnew<0.)
    fprintf(1, 'Newton-Pratt negative root:  x=%f\n', xnew);
    xnew = 0;
  end
end

% computing the circle parameters

DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];

end
