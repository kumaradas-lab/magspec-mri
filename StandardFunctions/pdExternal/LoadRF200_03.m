%% Settings for rf amplifier RF-200

HW.TX.ExtRFSN = 3;
HW.TX.ExtRFType = 'RF200';  % string with type identifier

HW.TX.ChannelDef = 2;         % Default TX Channel set to Tx2

% Amplifier
% HW.TX.Uout2PaUout(2) = 10^((-15.0+53.0-3.0)/20);  % 53.0 dB amplification @ 24.0 MHz, 1.0 V without switch
HW.TX.Uout2PaUout(2) = 10^((-15.0+52.2-3.0)/20);  % 52.2 dB amplification @ 24.0 MHz, 1.0 V with Switch

HW.TX.Max.PaUout(2) = 250;  % maximum peak output voltage in V
HW.TX.Def.PaUout(2) = 12.5;  % default PaUout(2)=50 V => 25 W

HW.TX.Max.Amplitude = [20,20]*1e-3;   % maximum B1+ in T
HW.TX.Def.Amplitude = [20,20]*1e-3;   % default B1+ in T

%% Coil 38mm
HW.RX2TXdeadTime = 3e-6;         % Receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 66e-6;        % Receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 2e-6;      % Unblank of RF-100 before TX pulse
HW.TX.BlankPostset = 2e-6;     % Blank of RF-100 after TX pulse

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 2e-6;   % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 40e-6;  % Blank of receiver after TX pulse
HW.TX.BlankAQ = 1;               % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
HW.RX.VGAGainDef = HW.RX.VGAGainMax/1; % reduce VGA gain to avoid saturation

%%
UseExtRFAmpSwitch = 1;  % use switch of RF-200
LoadExtRFAmp_Cal;  % new cal Uout and 6 A FET
