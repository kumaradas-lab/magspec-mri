%% Settings for rf amplifier RF-200

if ~exist('iDevice', 'var')
  iDevice = 1;
end

HW.TX(iDevice).ExtRFSN = 3;
HW.TX(iDevice).ExtRFType = 'RF200';  % string with type identifier

HW.TX(iDevice).ChannelDef = 2;  % default TX channel set to Tx2

% Amplifier
% HW.TX(iDevice).Uout2PaUout(2) = 10^((-15.0+53.0-3.0)/20);  % 53.0 dB amplification @ 24.0 MHz, 1.0 V without switch
HW.TX(iDevice).Uout2PaUout(2) = 10^((-15.0+52.2-3.0)/20);  % 52.2 dB amplification @ 24.0 MHz, 1.0 V with Switch

HW.TX(iDevice).Max.PaUout(2) = 250;  % maximum peak output voltage in V
HW.TX(iDevice).Def.PaUout(2) = 12.5;  % default PaUout(2)=50 V => 25 W

HW.TX(iDevice).Max.Amplitude = [20, 20]*1e-3;  % maximum B1+ in T
HW.TX(iDevice).Def.Amplitude = [20, 20]*1e-3;  % default B1+ in T

%% Coil 38mm
HW.RX2TXdeadTime = 3e-6;  % receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 66e-6;  % receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
HW.TX(iDevice).BlankOffset = 2e-6;  % unblank of RF-200 before TX pulse
HW.TX(iDevice).BlankPostset = 2e-6;  % blank of RF-200 after TX pulse

%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 40e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1;  % reduce VGA gain to avoid saturation

%%
UseExtRFAmpSwitch = 1;  % use switch of RF-200
LoadExtRFAmp_Cal;  % new cal Uout and 6 A FET
