%% Configure and initialize power supply with S/N 01

HW.PowerSupply.SN = 1;  % serial number (of settings file)

if isa(HW.PowerSupply, 'PD.PowerSupplyEA9000') && isvalid(HW.PowerSupply)
  % change to USB remote control to reset settings
  HW.PowerSupply.SetLock(true);
end

HW.PowerSupply.PaUoutMax = 155;  % maximum output voltage of amplifier in Volts
HW.PowerSupply.PaPoutMax = 1400;  % maximum output power of amplifier in Watts

HW.PowerSupply.PaIoutMaxAllowed = 10;  % maximum allowed output current of amplifier in Amperes
HW.PowerSupply.PaUoutMaxAllowed = 155;  % maximum allowed output voltage of amplifier in Amperes
HW.PowerSupply.PaPoutMaxAllowed = 1400;  % maximum allowed output power of amplifier in Amperes

% protection settings of Power Supply (set manually!):
% OVP: 150 V
% OCP: 11 A
% OPP: 1300 W

HW.PowerSupply.writeCommandDelay = 0.1;  % time between sending a command and checking the corresponding state in seconds

HW.PowerSupply.EnableAnalogChannel = 1;  % number of the digital output channel for analog remote control
HW.PowerSupply.EnableAnalogTOffset = 0.1;  % offset time in seconds when the enable analog control signal should be switched on
HW.PowerSupply.EnableAnalogTPostset = 0.1;  % postset time in seconds when the enable analog control signal should be switched off

% Reduce the coil efficiency by this factor to cause an "over-drive" of the
% output current.
% Limit the output current digitally by the same factor to get a correct (and
% steady) current output.
HW.PowerSupply.AnalogOverDriveFactor = 1.1;  % over-drive by 10%

HW.PowerSupply.SelectAmplifierChannel = 5;  % number of the digital output channel that is used to switch amplifiers
% FIXME: Should we invert this channel already when the power supply is loaded?
%        Currently this is done in LoadGradAmp_RoPS_01.

if ~exist('PD.PowerSupplyEA9000', 'class')
  return;
end


if ~exist('iDevice', 'var'), iDevice = 1; end

% invert digital output that is used to enable analog control
if isemptyfield(HW.DigitalIO(iDevice), 'InvertChannelOut')
  oldInvertChannelOut = 0;
else
  oldInvertChannelOut = HW.DigitalIO(iDevice).InvertChannelOut;
end
HW.DigitalIO(iDevice).InvertChannelOut = bitor(oldInvertChannelOut, HW.PowerSupply.EnableAnalogChannel);

% Set the inverted channels *before* programming the power supply if necessary
if ~isequal(oldInvertChannelOut, HW.DigitalIO(iDevice).InvertChannelOut)
  talker = PD.Talker.GetInstance();
  if talker(iDevice).myMon.IsConfigured ~= 0
    talker(iDevice).CreateAbortCommands(HW, iDevice);
    talker(iDevice).myTX.sendCommand(talker(iDevice).abortCommandArray);
  end
end

if (~isa(HW, 'PD.HWClass') && isemptyfield(HW, 'PowerSupply')) || ~isa(HW.PowerSupply, 'PD.PowerSupplyEA9000') || ...
    isempty(HW.PowerSupply) || (isa(HW.PowerSupply, 'PD.PowerSupplyEA9000') && ~isvalid(HW.PowerSupply))
  % only re-create object if it was explicitly deleted before
  HW.PowerSupply = PD.PowerSupplyEA9000.GetInstance(HW.PowerSupply);
end


% settings for gradient output channel
iGradChannel = 4;  % gradient channel that controls the power supply

HW.Grad(iDevice).ExtGradSN = 999;
% HW.Grad(iDevice).PaEnable = 1;
HW.Grad(iDevice).PaCurrentControlled(iGradChannel) = 1;  % If current controlled, set to 1. If voltage controlled, set to 0
HW.Grad(iDevice).PaRin(iGradChannel) = 70e3;  % NT input impedance
HW.Grad(iDevice).PaRout(iGradChannel) = 1e6;  % power supply output impedance


% Set offset such that 0 outputs a slightly negative voltage.
% That is done to avoid that the power supply switches on and off repeatedly due
% to noise on the analog channel.
HW.Grad(iDevice).PaOffsetI(iGradChannel) = +0.030;

% amplification of power supply
HW.Grad(iDevice).PaUin2PaIout(iGradChannel) = 2.5033; % 2.5000;  % input voltage to output current ratio in A/V

% reduce coil efficiency (see HW.PowerSupply.AnalogOverDriveFactor above)
HW.Grad(iDevice).LoadIin2Amp(iGradChannel) = ...
  HW.Grad(iDevice).LoadIin2Amp(iGradChannel) / HW.PowerSupply.AnalogOverDriveFactor;

HW.Grad(iDevice).PaUoutMax(iGradChannel) = HW.PowerSupply.PaPoutMax;  % maximum output voltage (limited by openMatlab)
HW.Grad(iDevice).PaUoutMin(iGradChannel) = -HW.Grad(iDevice).PaUoutMax(iGradChannel);  % "maximum" negative output voltage (limited by openMatlab)
