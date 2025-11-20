%% Settings for LNA #31 21.5 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 31;

% calibrated with RF-100_06 HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 40e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 5e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 5e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 10e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/4;  % default receiver gain; for denominators approx < 3 -> receiver saturation
HW.RX(iDevice).LNAGain = 10^(((-35.9298)-(-60))/20);  % 24.0702 dB gain @ 21500325.633 MHz F=1.755 dB (-60 dBm cal) HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1; CalibrationRxNoise 0.26997 dB gain

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 35);  % def transmit voltage
