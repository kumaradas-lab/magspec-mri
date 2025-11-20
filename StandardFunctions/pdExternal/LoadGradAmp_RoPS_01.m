%% Load settings for a named gradient amplifier

if ~exist('talker', 'var')
  talker = PD.Talker.GetInstance();
end

for iDevice = 1:numel(talker)
  if isempty(HW.Grad(iDevice).PaName), continue; end

  if strcmp(HW.Grad(iDevice).PaName, '1-10 mm')
    % external power supply properties

    LoadPowerSupply_01;

    % invert channel that selects amplifier (relay)
    HW.DigitalIO(iDevice).InvertChannelOut = bitor(HW.DigitalIO(iDevice).InvertChannelOut, ...
      2^(HW.PowerSupply.SelectAmplifierChannel-1));

    HW.Grad(iDevice).CombineCurrentOutputs = [];  % channels, each row will be treated as a combined output

    HW.RecoveryCPMG.thicknessSliceRange = [1e-3, 10e-3];  % FIXME: Re-set to 1-10mm after commissioning

    HW.RecoveryCPMG.PowerSupply.EnableAnalogChannel = HW.PowerSupply.EnableAnalogChannel;

    % 145 W continuous power can be dissipated without reaching 60 degrees.
    HW.Grad(iDevice).CoilPowerDissipation(2) = 0.985*145*HW.PowerSupply.AnalogOverDriveFactor^2;  % power dissipation at maximum temperature in Watt (measured in Rimpar: 70 W)

    HW.Grad(iDevice).CoilMaxDcCurrent(4) = 8.2*HW.PowerSupply.AnalogOverDriveFactor*1.01;  % maximum (nominal) DC current of fuse in Ampere
    HW.Grad(iDevice).CoilCurrentSquareTime(4) = 40;  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec (datasheet: 20.23)

  elseif strcmp(HW.Grad(iDevice).PaName, '2-100 mm')
    % DC-600 properties used at BAM in Berlin
    % See also LoadGradAmp_DC600_SN_23024.m

    LoadGradAmp_DC600_SN_23024;

    if isa(HW.PowerSupply, 'PD.PowerSupplyEA9000')
      % If the power supply hasn't been loaded yet, the values for configuration
      % are still unknown. Leave this setting unchanged for now.
      HW.DigitalIO(iDevice).InvertChannelOut = bitand(HW.DigitalIO(iDevice).InvertChannelOut, ...
        bitxor(2^8-1, 2^(HW.PowerSupply.SelectAmplifierChannel-1)));
    end

    % do not invert channel that selects amplifier (relay)
    HW.Grad(iDevice).CombineCurrentOutputs = [];  % channels, each row will be treated as a combined output

    HW.RecoveryCPMG.thicknessSliceRange = [2e-3, 100e-3];

    HW.RecoveryCPMG.PowerSupply.EnableAnalogChannel = 0;

    HW.Grad(iDevice).CoilPowerDissipation(2) = 2*65;  % power dissipation at maximum temperature in Watt (measured in Rimpar: 70 W)

    HW.Grad(iDevice).CoilMaxDcCurrent(4) = 4;  % maximum (nominal) DC current of fuse in Ampere
    HW.Grad(iDevice).CoilCurrentSquareTime(4) = 29.165 * 0.95;  % time-lag fuse parameter (maximum "accumulated heat") in A^2*sec

  else
    error('PD:LoadGradAmp:InvalidAmplifier', 'No settings for amplifier "%s".', ...
      HW.Grad(iDevice).PaName);
  end

end
