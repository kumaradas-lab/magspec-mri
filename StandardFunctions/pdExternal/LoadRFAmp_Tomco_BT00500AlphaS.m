%% RF Amplifier Tomco BT00500-AlphaS (with 20 dB input attenuation)
% Some settings might be overwritten by values in LoadLNA or switch files.
% Must be loaded before LoadLNA and LoadTxRx files.

HW.TX.ChannelDef = 2;                     % default channel: TX2
HW.TX.Uout2PaUout(2) = 10^((57-20)/20);   % 57 dB gain - 20 dB input attenuation
% HW.TX.Max.PaUout(2) = 10^(57/20)*sqrt(0.001*50*2);   % dBm to peak voltage, e.g. 57 dBm => 500 W
HW.TX.Max.PaUout(2) = sqrt(500*50*2);     % maximum peak voltage PaUout(2) = 100 V => 100 W, 224 V => 500W

HW.TX.Def.PaUout(2) = 1;                  % default peak voltage at TX2 10 V => 1 W, 100 V => 100 W, 224 V => 500 W

HW.TX.Max.Amplitude(2) = 1/1e-6/(HW.Gamma.H1/2/pi); % maximum B1+ amplitude in T, e.g. 1 microsecond for 360 degrees pulse => 23 mT
HW.TX.Def.Amplitude(2) = 1;               % default B1+ amplitude in T

%% Timing of TX (all in seconds)
HW.RX2TXdeadTime = 1e-6;                  % dead time RX to TX
HW.TX2RXdeadTime = 50e-6;                 % dead time TX to RX
HW.TX.BlankOffset = 800e-9;               % Offset of Blank before TX pulse
HW.TX.BlankPostset = 160e-9;              % delay of Blank after TX pulse

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 400e-9;             % Blank of receiver before TX pulse in seconds
HW.TX.BlankPostsetAQ = 400e-9;            % Blank of receiver after TX pulse in seconds
HW.TX.BlankAQ = 1;                        % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
