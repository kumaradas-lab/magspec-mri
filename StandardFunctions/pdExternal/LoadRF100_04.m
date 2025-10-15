%% Settings for rf amplifier RF-100 with SN 04

HW.TX.ExtRFSN = 4;
HW.TX.ExtRFType = 'RF100';  % string with type identifier

HW.TX.ChannelDef = 2;             % default TX channel set to Tx2

%% Amplifier settings
HW.TX.Uout2PaUout(2) = 10^(40/20);  % 36 dB to 40 dB amplification

HW.TX.Max.PaUout(2) = 250;  % maximum peak output voltage in V
HW.TX.Def.PaUout(2) = 20;  % default peak output voltage in V (might be reduced by other HW.TX.Def settings)

HW.TX.Max.Amplitude = [20, 20]*1e-3;  % maximum B1+ in T
HW.TX.Def.Amplitude = [20, 20]*1e-3;  % sensible default B1+ in T

%% Coil in magspec
%% Coil 10mm
HW.RX2TXdeadTime = 1e-6;         % receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 50e-6;        % receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 800e-9;      % unblank of RF-100 before TX pulse
HW.TX.BlankPostset = 400e-9;     % blank of RF-100 after TX pulse
%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 1000e-9;   % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 4000e-9;  % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;               % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
HW.RX.VGAGainDef = HW.RX.VGAGainMax/2;  % reduce VGA gain to avoid saturation

%%
UseExtRFAmpSwitch = 1;  % use switch of RF-100
LoadExtRFAmp_Cal;
