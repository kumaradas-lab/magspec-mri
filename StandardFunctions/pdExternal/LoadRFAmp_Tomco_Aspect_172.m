%% settings for Tomco rf amplifier used with Aspect system
% Some settings might be overwritten by values in LoadLNA or switch files.
% Must be loaded before LoadLNA and LoadTxRx files.

HW.TX(1).ExtRFSN = 10172;  % dummy serial number

HW.TX(1).ChannelDef = 2;  % default channel: TX2
% HW.TX(1).Uout2PaUout(2) = 10^((57-20)/20);  % 57 dB gain - 20 dB input attenuation
HW.TX(1).Uout2PaUout(2) = 10^((57-10)/20);  % 57 dB gain - 10 dB input attenuation
% HW.TX(1).Max.PaUout(2) = 10^(57/20)*sqrt(0.001*50*2);   % dBm to peak voltage, e.g. 57 dBm => 500 W
HW.TX(1).Max.PaUout(2) = sqrt(500*50*2);  % maximum peak voltage PaUout(2) = 100 V => 100 W, 224 V => 500W

HW.TX(1).Def.PaUout(2) = 135;  % default peak voltage at TX2

HW.TX(1).Max.Amplitude(2) = 1/1e-6/(HW.Gamma.H1/2/pi);  % maximum B1+ amplitude in T, e.g. 1 microsecond for 360 degrees pulse => 23 mT
HW.TX(1).Def.Amplitude(2) = 1;  % default B1+ amplitude in T  (high dummy value, so PaUout limit is used)

%% Timing of TX (all in seconds)
HW.RX2TXdeadTime = 1e-6;                  % dead time RX to TX
HW.TX2RXdeadTime = 50e-6;                 % dead time TX to RX
HW.TX(1).BlankOffset = 2000e-9;           % offset of Blank before TX pulse
HW.TX(1).BlankPostset = 160e-9;           % delay of Blank after TX pulse

%% TRx switch during transmit at Tx2
HW.TX(1).BlankOffsetAQ = 400e-9;          % blank of receiver before TX pulse in seconds
HW.TX(1).BlankPostsetAQ = 400e-9;         % blank of receiver after TX pulse in seconds
HW.TX(1).BlankAQ = 1;                     % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
