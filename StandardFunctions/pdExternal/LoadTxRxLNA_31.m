% Load settings for LNA #31 21.5 MHz
HW.RX.LnaSN = 31;

% calibrated with RF-100_06 HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;               % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 40e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 5e-6;               % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;              % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 5e-6;             % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 10e-6;           % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/4;  % default receiver gain; for denominators approx < 3 -> receiver saturation
HW.RX.LNAGain=10^(((-35.9298)-(-60))/20);  % 24.0702 dB gain @ 21500325.633 MHz F=1.755 dB (-60 dBm cal) HW.RX.VGAGainDef = HW.RX.VGAGainMax/1; CalibrationRxNoise 0.26997 dB gain

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 35);  % def transmit voltage
