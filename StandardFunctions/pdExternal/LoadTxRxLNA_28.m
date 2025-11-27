%% Settings for LNA #28 2-40 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 28;

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
HW.RX(iDevice).LNAGain = 10^(((-37.8982)-(-60))/20);  % 22.1018 dB gain @ 2300024.6129 MHz F=1.9084 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-37.6936)-(-60))/20);   % 22.3064 dB gain @ 5300056.7849 MHz F=1.3648 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-38.4415)-(-60))/20); % 21.5585 dB gain @ 12300131.8212 MHz F=1.4811 dB (-60 dBm cal)
HW.RX(iDevice).LNAGain = 10^(((-40.1317)-(-60))/20);  % 19.8683 dB gain @ 24300260.1229 MHz F=1.5357 dB (-60 dBm cal)

% when using internal switch of LNA (instead of switch in Rf-100 #15):
% HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 10);  % max transmit voltage
% HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 4);  % def transmit voltage
