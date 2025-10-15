%% Load settings for a named gradient amplifier

if ~exist('talker', 'var')
  talker = PD.Talker.GetInstance();
end

for iDevice = 1:numel(talker)
  if isempty(HW.Grad(iDevice).PaName), continue; end

  if strcmp(HW.Grad(iDevice).PaName, '1-10 mm')
    % external power supply properties

    combinedChannels = false;
    LoadGradSystem_300mT_01;

    LoadPowerSupply_01;

    HW.Grad.Slice.channel = 4;
    HW.Grad(iDevice).CoilMaxDcCurrent(4) = 6.5 * 0.9; % A
    HW.Grad(iDevice).CoilCurrentSquareTime(4) = 0.9 * 0.9; % A^2*sec

    % invert channel that selects amplifier (relay)
    HW.DigitalIO(iDevice).InvertChannelOut = bitor(HW.DigitalIO(iDevice).InvertChannelOut, ...
      2^(HW.PowerSupply.SelectAmplifierChannel-1));

    HW.RecoveryCPMG.thicknessSliceRange = [1e-3, 10e-3];

    HW.RecoveryCPMG.PowerSupply.EnableAnalogChannel = HW.PowerSupply.EnableAnalogChannel;

  elseif strcmp(HW.Grad(iDevice).PaName, '2-100 mm')
    % DC-600 properties used for 300mT magnet in Rimpar
    % See also LoadGradAmp_DC600_SN_23024.m

    LoadGradAmp_DC600_SN_38;

    combinedChannels = true;
    LoadGradSystem_300mT_01;

    HW.Grad.Slice.channel = 2;
    HW.Grad(iDevice).CoilMaxDcCurrent([2 4]) = 6.5 * 0.9; % A
    HW.Grad(iDevice).CoilCurrentSquareTime([2 4]) = 0.9 * 0.9; % A^2*sec

    % do not invert channel that selects amplifier (relay)
    HW.DigitalIO(iDevice).InvertChannelOut = bitand(HW.DigitalIO(iDevice).InvertChannelOut, ...
      bitxor(2^8-1, 2^(HW.PowerSupply.SelectAmplifierChannel-1)));

    HW.RecoveryCPMG.thicknessSliceRange = [2e-3, 100e-3];

    HW.RecoveryCPMG.PowerSupply.EnableAnalogChannel = 0;

  else
    error('PD:LoadGradAmp:InvalidAmplifier', 'No settings for amplifier "%s".', ...
      HW.Grad(iDevice).PaName);
  end

end
