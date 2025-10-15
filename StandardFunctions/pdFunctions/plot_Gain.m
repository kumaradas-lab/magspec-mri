function hf = plot_Gain(HW, Network, Seq)
%% Plot results of gain measurement
%
%   hf = plot_Gain(HW, Network, Seq)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% check input
if nargin <= 2
  hf = figure(205);
  hf = clf(hf);
  set(hf, 'Name', 'Gain and Phase');
  ax1 = subplot(2,1,1, 'Parent', hf);
  ax2 = subplot(2,1,2, 'Parent', hf);
end

if nargin == 3
  aherror = 0;

  if ~isempty(Seq.plotGainHandle)
    ax1 = Seq.plotGainHandle;
  else
    aherror = 1;
  end
  if ~isempty(Seq.plotPhaseHandle)
    ax2 = Seq.plotPhaseHandle;
  else
    aherror = 1;
  end
  if aherror
    hf = figure(205);
    hf = clf(hf);
    set(hf, 'Name', 'Gain and Phase');
    ax1 = subplot(2,1,1, 'Parent', hf);
    ax2 = subplot(2,1,2, 'Parent', hf);
  end
end


%% plot amplitude
if min(Network.FrequencyGain) == max(Network.FrequencyGain)
  plot(ax1, Network.FrequencyGain, 20*log10(abs(Network.Gain)), '*');
  ylim(ax1, [floor(min(20*log10(abs(Network.Gain)))/5)*5, ceil(max(20*log10(abs(Network.Gain)))/5)*5]);
  xlim(ax1, [min(Network.FrequencyGain)-1, max(Network.FrequencyGain)+1]);
  set(ax1, 'XTick', min(Network.FrequencyGain));
  set(ax1, 'XTickLabel', {num2str(min(Network.FrequencyGain),'%12.1f')});
else
  plot(ax1, Network.FrequencyGain, 20*log10(abs(Network.Gain)));
  ylim(ax1, [floor(min(20*log10(abs(Network.Gain)))/5)*5, ceil(max(20*log10(abs(Network.Gain)))/5)*5]);
  xlim(ax1, [min(Network.FrequencyGain), max(Network.FrequencyGain)]);
  resonance_gain_dB = 20*log10(abs(interp1(Network.FrequencyGain, double(Network.Gain), HW.fLarmor)));
  ht = text(HW.fLarmor, resonance_gain_dB, ...
    ['\leftarrow fL ' num2str(resonance_gain_dB, '%.1f') ' dB @ ' num2str(HW.fLarmor/1e6, '%.3f') ' MHz'], ...
    'HorizontalAlignment', 'left', 'Parent', ax1);
  if (ceil(max(20*log10(abs(Network.Gain)))/5)*5 + floor(min(20*log10(abs(Network.Gain)))/5)*5)/2 ...
      > resonance_gain_dB
    set(ht, 'Rotation', 90);
  else
    set(ht, 'Rotation', -90);
  end
end
ylabel(ax1, {'Gain', 'in dB'});
xlabel(ax1, 'Frequency in Hz');
set(ax1, 'XGrid', 'on');
set(ax1, 'YGrid', 'on');


%% plot phase
if min(Network.FrequencyGain) == max(Network.FrequencyGain)
  plot(ax2, Network.FrequencyGain, angle(Network.Gain), '*');
  xlim(ax2, [min(Network.FrequencyGain)-1, max(Network.FrequencyGain)+1]);
  set(ax2, 'XTick', min(Network.FrequencyGain));
  set(ax2, 'XTickLabel', {num2str(min(Network.FrequencyGain),'%12.1f')});
else
  plot(ax2, Network.FrequencyGain, angle(Network.Gain));
  xlim(ax2, [min(Network.FrequencyGain), max(Network.FrequencyGain)]);
  text(HW.fLarmor, ...
    angle(interp1(Network.FrequencyGain, double(Network.Gain), HW.fLarmor)), ...
    '\leftarrow fLarmor', ...
    'HorizontalAlignment', 'left', ...
    'Rotation', 90, 'Parent', ax2);
end
ylabel(ax2, {'Angle', 'in rad'});
xlabel(ax2, 'Frequency in Hz');
set(ax2, 'XGrid', 'on');
set(ax2, 'YGrid', 'on');
% drawnow


end
