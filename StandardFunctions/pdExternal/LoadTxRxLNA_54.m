%% Settings for high power L/4 (13.1 MHz) LNA #53
% built Sep 2023
% Built-in into RF-200 #04

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 54;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/7;  % default receiver gain

HW.RX(iDevice).LNAGain = 10^(((-48.3376)-(-70.144))/20);  % 21.8064 dB gain @ 13100246.4924 MHz F=2.2135 dB (-70.144 dBm cal)
 
HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 40);  % def transmit voltage

% LoadExtLNA_Cal;
% Do not load when used with RF-200. Calibration done with LowPower-script.
