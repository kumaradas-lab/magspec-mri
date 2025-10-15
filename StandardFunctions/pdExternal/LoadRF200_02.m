%% Settings for rf amplifier RF-200

HW.TX.ExtRFSN = 2;
HW.TX.ExtRFType = 'RF200';  % string with type identifier

HW.TX.ChannelDef = 2;  % Default TX Channel set to Tx2

% Amplifier
HW.TX.Uout2PaUout(2) = 10^((-15.0+53.0-3.0)/20);  % 53.0 dB amplification @ 24.0 MHz, 1.0 V amp only
HW.TX.Uout2PaUout(2) = 10^((-15.0+52.2-3.0)/20);  % 52.2 dB amplification @ 24.0 MHz, 1.0 V with Switch

HW.TX.Max.PaUout(2) = 250;  % maximum peak output voltage in V
HW.TX.Def.PaUout(2) = 80;  % default PaUout(2)=50 V => 25 W

% switch HW.UserName
%   case {'10mm_single', '10mm_both'}
%     HW.TX.Def.PaUout(2) = 50;  % 10 mm Coil: 50V
%   case {'15mm_single', '15mm_both', '14.8mm_light_single', '15mm_light_single', ...
%       '15mm_light_single_b1', '15mm_light_single_b1_DC600', 'teach'}
%     HW.TX.Def.PaUout(2) = 75;  % 15 mm Coil: 75V
%   otherwise
%     error('PD:LoadRF200:UnknownUser', ...
%       'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
% end

HW.TX.Max.Amplitude = [20,20]*1e-3;  % maximum B1+ in T
HW.TX.Def.Amplitude = [20,20]*1e-3;  % default B1+ in T

%%
UseExtRFAmpSwitch = 2; % 1: use switch of RF-200; 2: external switch in probe

LoadExtRFAmp_Cal;  % new cal Uout and 6 A FET

%% Coil 15mm (long)
HW.RX2TXdeadTime = 2.4e-6;  % receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 10e-6;  % receiver deadtime after TX pulse in s, e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 2000e-9;  % unblank of RF-200 before TX pulse
HW.TX.BlankPostset = 1e-6;  % blank of RF-200 after TX pulse

% rf coil damping Settings
HW.TX.DampCoil.Enable = true;  % enable coil damping
HW.TX.DampCoil.DigitalOutputLatency = 0.6e-6;  % latency of damping circuit in seconds
HW.TX.DampCoil.DigitalOutputChannel = 1;  % digital output channel for coil damping signal

% switch HW.UserName
%   case {'14.8mm_light_single'}
%     HW.TX.DampCoil.DigitalOutputDuration = 2.0e-6; % duration of signal on digital out in seconds
%     HW.TX.DampCoil.DampingDuration = 2.0e-6; % effective duration of damping the coil in seconds
%   case {'15mm_light_single'}
%     HW.TX.DampCoil.DigitalOutputDuration = 2.4e-6; % duration of signal on digital out in seconds
%     HW.TX.DampCoil.DampingDuration = HW.TX.DampCoil.DigitalOutputDuration+2e-6; % effective duration of damping the coil in seconds
%   case {'15mm_light_single_b1', '15mm_light_single_b1_DC600', 'teach'}
%     HW.TX.DampCoil.DigitalOutputDuration = 1.4e-6; % duration of signal on digital out in seconds
%     HW.TX.DampCoil.DampingDuration = HW.TX.DampCoil.DigitalOutputDuration+2e-6; % effective duration of damping the coil in seconds
%   otherwise
%     HW.TX.DampCoil.DigitalOutputDuration = 1.4e-6; % duration of signal on digital out in seconds
%     HW.TX.DampCoil.DampingDuration = 1.4e-6; % effective duration of damping the coil in seconds
% end

HW.TX.DampCoil.DigitalOutputDuration = HW.TX.BlankPostset + 1e-6;  % duration of signal on digital out in seconds
HW.TX.DampCoil.DampingDuration = HW.TX.DampCoil.DigitalOutputDuration + 2e-6;  % effective duration of damping the coil in seconds

HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 3e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 2000e-9;  % Blank of receiver before TX pulse
if HW.TX.DampCoil.Enable
  HW.TX.BlankPostsetAQ = HW.TX.DampCoil.DigitalOutputDuration + 1.0e-6;  % blank of receiver after TX pulse
else
  HW.TX.BlankPostsetAQ = 6000e-9;
end
HW.TX.BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% LNA Settings
% HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % reduce VGA gain to avoid saturation
%
% switch HW.UserName
%   case {'10mm_single', '10mm_both'}
%     HW.RX.LNAGain=10^(((-51.7175)-(-73.75))/20);  % 22.0325 dB gain @ 24300494.5011 MHz F=1.4714 dB (-73.75 dBm cal)
%   case {'15mm_single', '15mm_both', '14.8mm_light_single'}
%     HW.RX.LNAGain=10^(((-51.3859)-(-73.54))/20);  % 22.1541 dB gain @ 21500436.5131 MHz F=1.6048 dB (-73.54 dBm cal) 21.50+22.15+1.60
%   case {'15mm_light_single'}
%     HW.RX.LNAGain=10^(((-49.3502)-(-70.8655))/20);  % 21.5153 dB gain @ 21490363.0437 MHz F=1.4506 dB (-70.8655 dBm cal) with L/4 switch % in LoadTxRxLNA_50;
%   case {'15mm_light_single_b1', '15mm_light_single_b1_DC600', 'teach'}
%     HW.RX.LNAGain=10^(((-49.3502)-(-70.8655))/20);  % 21.5153 dB gain @ 21490363.0437 MHz F=1.4506 dB (-70.8655 dBm cal) with L/4 switch % in LoadTxRxLNA_50;
%   otherwise
%     error('PD:LoadRF200:UnknownUser', ...
%       'No settings for user "%s" specified in "%s"', HW.UserName, mfilename());
%     % center frequency: HW.RX.LNAGain=10^(((-51.5027)-(-73.67))/20); % 22.1673 dB gain @ 23200472.8774 MHz F=1.5127 dB (-73.67 dBm cal) 23.20+22.17+1.51
% end
