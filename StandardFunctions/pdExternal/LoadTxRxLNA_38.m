% Load settings for LNA #38 22.9 MHz und 24.3 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 38;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.3/20);  % -0.3 dB Gain (TX to Coil) @ 24.15 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

% HW.RX(iDevice).LNAGain = 10^(((-44.9731)-(-67.633))/20);  % 22.6599 dB gain @ 22900331.5471 MHz F=0.78729 dB (-67.633 dBm cal)
% HW.RX(iDevice).LNAGain = 10^(((-45.2458)-(-67.633))/20);  % 22.3872 dB gain @ 24300351.3135 MHz F=0.71811 dB (-67.633 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-44.8925)-(-67.633))/20);  % 22.7405 dB gain @ 22900330.2545 MHz F=0.66069 dB (-67.633 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-45.1655)-(-67.633))/20);  % 22.4675 dB gain @ 24300350.7268 MHz F=0.70906 dB (-67.633 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 72);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 50);  % def transmit voltage
