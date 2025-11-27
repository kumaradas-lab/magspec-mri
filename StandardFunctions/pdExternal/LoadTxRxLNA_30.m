%% Settings for LNA #30 2-40 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 30;

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
HW.RX(iDevice).LNAGain = 10^(((-37.8162)-(-60))/20);  % 22.1838 dB gain @ 2300025.7339 MHz F=1.9949 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-37.6362)-(-60))/20);  % 22.3638 dB gain @ 5300059.0745 MHz F=1.3271 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-38.4329)-(-60))/20);  % 21.5671 dB gain @ 12300136.7997 MHz F=1.4805 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-40.1127)-(-60))/20);  % 19.8873 dB gain @ 24300269.5874 MHz F=1.4805 dB (-60 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 10);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 4);  % def transmit voltage
