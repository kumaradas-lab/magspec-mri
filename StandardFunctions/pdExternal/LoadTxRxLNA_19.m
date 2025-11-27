%% Settings for LNA #19 23.5 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 19;

% HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.4/20);  % -0.4 dB Gain @ 23.5 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 1e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 4e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 1e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 480e-9;  % blank of transmit after TX pulse

if 0
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % reduce VGA gain to avoid saturation
  % HW.RX(iDevice).LNAGain = 10^(((-35.8684)-(-60))/20);  % 24.1316 dB gain @ 23550303.9158 MHz F=0.85376 dB (-60 dBm cal) (without RF-100)
  HW.RX(iDevice).LNAGain = 10^(((-36.8717)-(-60))/20);  % 23.1283 dB gain @ 23500381.9362 MHz F=2.0899 dB (-60 dBm cal)(with RF-100)
else
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % reduce VGA gain to avoid saturation
  HW.RX(iDevice).LNAGain = 10^(((-36.9085)-(-60))/20);  % 23.0915 dB gain @ 23500382.707 MHz F=2.2518 dB (-60 dBm cal)
end

HW.TX(iDevice).BlankOffsetAQ = 400e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 400e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 50);  % def transmit voltage
