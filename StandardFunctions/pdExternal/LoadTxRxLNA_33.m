% Load settings COPY FROM LNA #32 24.5 MHz
HW.RX.LnaSN = 33;

% calibrated with RF-100_06 HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 1e-6; % >7µs         % blank of transmit before TX pulse
HW.TX.BlankPostset = 1e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 1e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 1e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

if 0
  HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain
  HW.RX.LNAGain=10^(((-34.1243)-(-60))/20); % 25.8757 dB gain @ 24500397.7985 MHz F=0.82651 dB (-60 dBm cal)
else
  HW.RX.VGAGainDef = HW.RX.VGAGainMax/6;  % default receiver gain
  HW.RX.LNAGain=10^(((-34.1668)-(-60))/20); % 25.8332 dB gain @ 24500396.855 MHz F=0.92902 dB (-60 dBm cal)
end
 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 101);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 30);  % def transmit voltage

LoadExtLNA_Cal