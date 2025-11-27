%% Settings for LNA #50 21.5 MHz

% Low power broadband LNA built Mar 2023 combined with passive L/4 switch

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 50;

% HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 3e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 10e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 2e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 1e-6;  % blank of transmit after TX pulse

% HW.TX.BlankOffsetAQ = 2.4e-6;  % blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ = 2.4e-6;  % blank of receiver after TX pulse
% HW.TX.BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

% HW.RX(iDevice).LNAGain = 10^(((-49.5269)-(-71.24))/20);  % 21.7131 dB gain @ 21500395.1886 MHz F=1.1788 dB (-71.24 dBm cal) without L/4 switch
HW.RX(iDevice).LNAGain = 10^(((-49.3502)-(-70.8655))/20);  % 21.5153 dB gain @ 21490363.0437 MHz F=1.4506 dB (-70.8655 dBm cal) with L/4 switch

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = 75;  % def transmit voltage

% LoadExtLNA_Cal;  % only direct to TX2
