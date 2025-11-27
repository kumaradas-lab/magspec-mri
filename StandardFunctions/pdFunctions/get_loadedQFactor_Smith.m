function Network = get_loadedQFactor_Smith(Network)
%% Calculate Q-Factor from relection measurement using Smith-Chart geometry
%
%     Network = get_loadedQFactor_Smith(Network)
%
% Calculate the loaded Q-factor of the coil with the method described in
% https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf
%
%
% INPUT:
%
%   Network
%       Structure with at least the following fields:
%     Frequency
%         Vector with frequencies in Hz.
%     Reflection
%         Complex reflection ratio.
%
%
% OUTPUT:
%
%   Network
%       Structure as in input with the following additional fields:
%     minRefFrequencySSInt
%         Frequency at the point in the Smith chart that is closest to the 50
%         Ohm point (determined from a circle fit) in Hertz. I.e., the resonant
%         frequency of the connected coil.
%     minRefReflection_dB
%         Reflection at the point in the Smith chart that is closest to the 50
%         Ohm point (determined from a circle fit) in Dezibel. I.e., the
%         reflection at the resonant frequency of the connected coil.
%     QLFrequencySSInt
%         Loaded Q-factor as determined from a circle fit.
%     Q_iStart
%         Start index into Network.Reflection for the range that was used to fit
%         the circle into the Smith chart.
%     Q_iStop
%         End index into Network.Reflection for the range that was used to fit
%         the circle into the Smith chart.
%     Q_circle_center
%         Best fit coordinate for the center of the circle in the Smith chart
%         (complex valued).
%     Q_circle_radius
%         Best fit radius for the circle in the Smith chart.
%     Q_circle_minCenter
%         Coordinate of the point on the best fit circle that is closest to the
%         center of the Smith chart (i.e., the 50 Ohm point).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Loaded Q-factor (https://engineering.olemiss.edu/~eedarko/experience/rfqmeas2b.pdf):
% "Draw line" from (0,0) trough fMinReflection
% Find frequencies (f1,f2) of points on S11 (Q-circle) at 45 degrees from the
% "outer edge" of the Q-circle (i.e. furthest from that line)
% --> Q = fMinReflection / |f1 - f2|
% That direct approach is rather unstable for noisy data. Instead using the
% following equivalent approach:

% Find circle that fits points around min. reflection.
% tempRef = sort(abs(Reflection));
% iStart = find(abs(Reflection) < tempRef(round(numel(Reflection)/10)), 1, 'first');
% iStop  = find(abs(Reflection) < tempRef(round(numel(Reflection)/10)), 1, 'last');
FrequencyQ = Network.Frequency;
ReflectionQ = Network.Reflection;
[minRef, minRefI] = min(abs(ReflectionQ));
maxRef = max(abs(ReflectionQ));

Network.minRefFrequencySSInt  = NaN;
Network.minRefReflection_dB   = NaN;
Network.QLFrequencySSInt      = NaN;

if maxRef/minRef<sqrt(2)
  warning('PD:get_QFactor_Smith:TooFewData', ...
    'The input data does likely not cover the coil resonance');
  return;
end

for t=1:2
  switch t
    case 1
      iStart = max([1,minRefI+1-find(abs(ReflectionQ(minRefI:-1:1))                 > max([sqrt(0.6),maxRef-(maxRef-minRef)/2]), 1, 'first')]);
      iStop  = min([numel(ReflectionQ),minRefI-1+find(abs(ReflectionQ(minRefI:end)) > max([sqrt(0.6),maxRef-(maxRef-minRef)/2]), 1, 'first')]);
      Par = CircleFitByPratt([reshape(real(ReflectionQ(iStart:iStop)),[],1), reshape(imag(ReflectionQ(iStart:iStop)),[],1)]);

    case 2
      iStart = max([1,minRefI+1-find(angle(divReflection(minRefI:-1:1)) < -pi*120/180, 1, 'first')]);                 % from angle -120 degrees
      iStop  = min([numel(ReflectionQ),minRefI-1+find(angle(divReflection(minRefI:end)) > pi*120/180, 1, 'first')]);  % to angle 120 degrees
      weigthing = round((pi-abs(angle(divReflection)))/pi*max([round(2^20/numel(iStart:iStop)),10]));  % linear angle dependent weight factor from 10 to 0 at RAD 0 to +/-pi
      wReflectionQ = NaN(sum(weigthing(iStart:iStop)), 1);
      tti = 1;
      for tt=iStart:iStop
        % angle dependent weight (add weight dependent copys of points)
        wReflectionQ(tti:tti+weigthing(tt)-1) = ReflectionQ(tt);
        tti = tti + weigthing(tt);
      end
      % Circle with angle dependent weight (<1e6 points < 10ms)
      Par = CircleFitByPratt([real(wReflectionQ), imag(wReflectionQ)]);
  end
  xCenterQ = Par(1);
  yCenterQ = Par(2);
  rCenterQ = Par(3);

  % find point with minimum distance to center of Smith chart
  % on line from the center of the Smith chart through the
  % center of the Q-circle

  reflCenterQ = sqrt(xCenterQ^2 + yCenterQ^2);
  xMinRefCircle = (reflCenterQ-rCenterQ)/reflCenterQ * xCenterQ;
  yMinRefCircle = (reflCenterQ-rCenterQ)/reflCenterQ * yCenterQ;
  minRefCircle = xMinRefCircle + 1i * yMinRefCircle;
  centerCircle = xCenterQ + 1i * yCenterQ;
  centerCircle2minRefCircle = minRefCircle - centerCircle;
  centerCircle2Reflection = ReflectionQ - centerCircle;
  divReflection = centerCircle2minRefCircle ./ centerCircle2Reflection;

  %   figure;plot(abs(d))
  %   figure;plot(angle(d))
  %   figure;plot(abs(divReflection(iStart:iStop)))
  %   figure;plot(angle(divReflection(iStart:iStop)))
end

FrequencySS = FrequencyQ(iStart:iStop);
angleAtFrequency = angle(divReflection(iStart:iStop));
[p, S, mu] = polyfit(FrequencySS, angleAtFrequency, ...
  min([floor(numel(angleAtFrequency)/2)*2-1, 7]));
% interpolate frequency about 1Hz resolution
FrequencySSInt = interpn(FrequencySS, ceil(log2((FrequencySS(2)-FrequencySS(1))/2)));
% get interpolated angle to interpolated frequencies
angleAtFrequencySSInt = polyval(p, FrequencySSInt, S, mu);
% figure(4); plot(FrequencySS,angleAtFrequency,FrequencySSInt,angleAtFrequencySSInt)
% find frequency of minRefCycle (0 RAD)
[~, minRefFrequencySSIntI] = min(abs(mod(angleAtFrequencySSInt      +pi, pi*2)-pi));
minRefFrequencySSInt = FrequencySSInt(minRefFrequencySSIntI);
% find frequency below minRefCycle (pi/2 RAD)
[~, minF1FrequencySSIntI] = min(abs(mod(angleAtFrequencySSInt+pi/2  +pi, pi*2)-pi));
minF1FrequencySSInt = FrequencySSInt(minF1FrequencySSIntI);
% find frequency above minRefCycle (+pi/2 RAD)
[~, minF2FrequencySSIntI] = min(abs(mod(angleAtFrequencySSInt-pi/2  +pi, pi*2)-pi));
minF2FrequencySSInt = FrequencySSInt(minF2FrequencySSIntI);

if (minRefFrequencySSIntI == 1 || minRefFrequencySSIntI == numel(FrequencySSInt))
  % repeat polyfit with reduced range for better extrapolation
  [p2, S2, mu2] = polyfit(FrequencySS, angleAtFrequency, ...
    min([floor(numel(angleAtFrequency)/2)*2-1, 3]));
  FrequencySS2 = linspace(...
    FrequencyQ(iStart)-(FrequencyQ(iStop)-FrequencyQ(iStart))/2, ...
    FrequencyQ(iStop)+(FrequencyQ(iStop)-FrequencyQ(iStart))/2, ...
    1000).';
  % interpolate frequency with about 1 Hz resolution
  FrequencySSInt2 = interpn(FrequencySS2, ...
    ceil(log2((FrequencySS2(2)-FrequencySS2(1))/2)));
  % get interpolated angle to interpolated frequencies
  angleAtFrequencySSInt2 = polyval(p2, FrequencySSInt2, S2, mu2);
  % find frequency of minRefCycle (0 RAD)
  [~, minRefFrequencySSIntI2] = min(abs(mod(angleAtFrequencySSInt2      +pi,pi*2)-pi));
  minRefFrequencySSInt=FrequencySSInt2(minRefFrequencySSIntI2);
  % figure(4);hold on; plot(FrequencySSInt2,angleAtFrequencySSInt2,'--');hold off;
end

if (minF1FrequencySSIntI == 1 || minF1FrequencySSIntI == numel(FrequencySSInt)) && ...
    (minF2FrequencySSIntI == 1 || minF2FrequencySSIntI == numel(FrequencySSInt))
  % repeated polyfit with reduced range for better extrapolation
  [p, S, mu] = polyfit(FrequencySS, angleAtFrequency, ...
    min([floor(numel(angleAtFrequency)/2)*2-1, 3]));
  FrequencySS = linspace(...
    FrequencyQ(iStart)-(FrequencyQ(iStop)-FrequencyQ(iStart))/2, ...
    FrequencyQ(iStop)+(FrequencyQ(iStop)-FrequencyQ(iStart))/2, ...
    1000).';
  % interpolate frequency with about 1 Hz resolution
  FrequencySSInt = interpn(FrequencySS, ...
    ceil(log2((FrequencySS(2)-FrequencySS(1))/2)));
  % get interpolated angle to interpolated frequencies
  angleAtFrequencySSInt = polyval(p, FrequencySSInt, S, mu);
  % find frequency below minRefCycle (pi/2 RAD)
  [~, minF1FrequencySSIntI] = min(abs(mod(angleAtFrequencySSInt+pi/2  +pi, pi*2)-pi));
  minF1FrequencySSInt = FrequencySSInt(minF1FrequencySSIntI);
  % find frequency above minRefCycle (+pi/2 RAD)
  [~, minF2FrequencySSIntI] = min(abs(mod(angleAtFrequencySSInt-pi/2  +pi, pi*2)-pi));
  minF2FrequencySSInt = FrequencySSInt(minF2FrequencySSIntI);
  % figure(4);hold on; plot(FrequencySSInt,angleAtFrequencySSInt,':');hold off;
end

if (minF1FrequencySSIntI == 1 || minF1FrequencySSIntI == numel(FrequencySSInt)) && ...
    (minF2FrequencySSIntI == 1 || minF2FrequencySSIntI == numel(FrequencySSInt))
  dfFrequencySSInt = [];
elseif minF1FrequencySSIntI == 1 || minF1FrequencySSIntI == numel(FrequencySSInt)
  dfFrequencySSInt = 2*abs(minF2FrequencySSInt - minRefFrequencySSInt);
elseif minF2FrequencySSIntI == 1 || minF2FrequencySSIntI == numel(FrequencySSInt)
  dfFrequencySSInt = 2*abs(minF1FrequencySSInt - minRefFrequencySSInt);
else
  dfFrequencySSInt = abs(minF1FrequencySSInt - minF2FrequencySSInt);
end

% calculate loaded Q factor and construct string for plot

if ~isempty(dfFrequencySSInt)
  QLFrequencySSInt = minRefFrequencySSInt / dfFrequencySSInt;
else
  QLFrequencySSInt = NaN;
end

fprintf('Frequency = %.0f Hz, Reflection = %.2f dB, Q = %.2f\n', ...
  minRefFrequencySSInt, ...
  20*log10(abs(minRefCircle)), ...
  QLFrequencySSInt);
% disp(['Frequency = ' num2str(minRefFrequencySSInt,'%10.0f') ' Hz, Reflection = ' num2str(20*log10(abs(minRefCircle)),'%6.2f') ' dB, Q = ' num2str(QLFrequencySSInt,'%6.2f')])

Network.minRefFrequencySSInt = minRefFrequencySSInt;
Network.minRefReflection_dB = 20*log10(abs(minRefCircle));
Network.QLFrequencySSInt = QLFrequencySSInt;
Network.Q_iStart = iStart;
Network.Q_iStop = iStop;
Network.Q_circle_radius = Par(3);
Network.Q_circle_minCenter = minRefCircle;
Network.Q_circle_center = centerCircle;

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

n = size(XY,1);  % number of data points

centroid = mean(XY);  % the centroid of the data set

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
Center = [Mxz*(Myy-xnew)-Myz*Mxy, Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

Par = [Center+centroid, sqrt(Center*Center'+Mz+2*xnew)];

end
