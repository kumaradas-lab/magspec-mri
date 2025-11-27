function [fCoil] = get_fCoil(HW)
%% Measure resonance frequency of rf coil using S11 reflection
%
%   [fCoil] = get_fCoil(HW)
%
%
% INPUT:
%
%   HW
%       HW object or structure
%
%
% OUTPUT:
%
%   fCoil
%       Resonance frequency of the rf coil (minimum reflection) in Hertz.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

Seq.plotSmith = 0;
Seq.plotNyquist = 0;
Seq.plotReflection = 1;
Seq.plotReflectionRaw = 0;
Seq.Cal = 0;

% measure with large bandwidth around current HW.fLarmor
Seq.fCenter = HW.fLarmor;
Seq.fSpan = 1e6;
Seq.fSteps = 101;

[Network, ~] = sequence_Network(HW, Seq);

if (20*log10(abs(Network.MinReflection)) < -6)
  % measure with smaller bandwidth around minimum reflection
  HW.fLarmor = Network.fMinReflection;

  Seq.fCenter = HW.fLarmor;
  Seq.fSpan = 0.05e6;
  Seq.fSteps = 301;
  [Network, ~] = sequence_Network(HW, Seq);

  if (20*log10(abs(Network.MinReflection)) < -6)
    % fCoil = Network.fMinReflection + Seq.fSpan/Seq.fSteps*rand(1);
    try
      myP = polyfit(...
        (Network.Frequency-Network.fMinReflection) ./ diff(Network.Frequency([1,end])), ...
        abs(Network.Reflection)-abs(Network.MinReflection), 12);
      fCoil = fminbnd(@(x) polyval(myP, x), -0.1, 0.1) ...
        * diff(Network.Frequency([1,end])) + Network.fMinReflection;
    catch
      fCoil = Network.fMinReflection + Seq.fSpan/Seq.fSteps*(rand(1)-0.5);
    end
  else
    warning('PD:get_fCoil:NoMinimumReflection', 'Resonance frequency not found.');
    fCoil = HW.fLarmor;
  end

else
  warning('PD:get_fCoil:NoMinimumReflection', 'Resonance frequency not found.');
  fCoil = HW.fLarmor;
end

end
