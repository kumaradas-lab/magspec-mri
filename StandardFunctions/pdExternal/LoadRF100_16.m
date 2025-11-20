%% Settings for rf amplifier RF-100 with SN 16

HW.TX.ExtRFSN = 16;
HW.TX.ExtRFType = 'RF100';  % string with type identifier

HW.TX.ChannelDef = 2;         % Default TX Channel set to Tx2

% Amplifier
HW.TX.Uout2PaUout(2) = 50;  % 50x amplification
HW.TX.Max.PaUout(2) = 100;  % maximum peak output voltage in V

HW.TX.Def.PaUout(2) = 10;  % default peak output voltage in V (might be reduced by other HW.TX.Def settings)

HW.TX.Max.Amplitude = [20, 20]*1e-3;  % maximum B1+ in T
HW.TX.Def.Amplitude = [20, 20]*1e-3;  % default B1+ in T

%% Coil 10mm
HW.RX2TXdeadTime = 1e-6;         % Receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 50e-6;        % Receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 800e-9;      % Unblank of RF-100 before TX pulse
HW.TX.BlankPostset = 400e-9;     % Blank of RF-100 after TX pulse
%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 1000e-9;   % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 4000e-9;  % Blank of receiver after TX pulse
HW.TX.BlankAQ = 1;               % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
HW.RX.VGAGainDef = HW.RX.VGAGainMax/2; % reduce VGA gain to avoid saturation

%%
UseExtRFAmpSwitch = 1;  % use switch of RF-100
LoadExtRFAmp_Cal;
