% Low power broadband LNA modified to be used in conjuction with a
% L/4-switch at 22.4MHz
HW.RX.LnaSN = 65;

% Settings in calibration file for Rf-200 #07
% HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
% HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
% HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
% HW.TX.BlankOffset = 9e-6; % >7µs         % blank of transmit before TX pulse
% HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse
% 
% HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
% HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain

HW.RX.LNAGain = 10^(((-50.2719)-(-71.14))/20);  % 20.8681 dB gain @ 22400423.3518 MHz F=1.2861 dB (-71.14 dBm cal)

% HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 10);  % max transmit voltage
% HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 4);  % def transmit voltage

% LoadExtLNA_Cal;