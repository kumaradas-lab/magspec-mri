% High power L/4 (13.1 MHz) LNA built Sep 2023
% Built-in into RF-200 #04
HW.RX.LnaSN = 54;

HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 9e-6; % >7µs         % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain
HW.RX.VGAGainDef = HW.RX.VGAGainMax/7;  % default receiver gain

HW.RX.LNAGain=10^(((-48.3376)-(-70.144))/20); % 21.8064 dB gain @ 13100246.4924 MHz F=2.2135 dB (-70.144 dBm cal)
 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 40);  % def transmit voltage

% LoadExtLNA_Cal;
% Don't Load when used with RRF-200. Calibration done with LowPower-Script.