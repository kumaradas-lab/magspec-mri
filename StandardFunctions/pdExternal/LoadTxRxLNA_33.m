%% Settings for LNA #33 - COPY FROM LNA #32 24.5 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 33;

% calibrated with RF-100_06 HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 1e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 1e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 1e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 1e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

if 0
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain
  HW.RX(iDevice).LNAGain = 10^(((-34.1243)-(-60))/20);  % 25.8757 dB gain @ 24500397.7985 MHz F=0.82651 dB (-60 dBm cal)
else
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % default receiver gain
  HW.RX(iDevice).LNAGain = 10^(((-34.1668)-(-60))/20);  % 25.8332 dB gain @ 24500396.855 MHz F=0.92902 dB (-60 dBm cal)
end

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 101);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 30);  % def transmit voltage

LoadExtLNA_Cal;
