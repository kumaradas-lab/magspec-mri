% Load settings for LNA #24 24.377452 MHz
HW.RX.LnaSN = 24;

% HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.4/20); % -0.4 dB Gain @ 23.5 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 1e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 4e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 1e-6;               % blank of transmit before TX pulse
HW.TX.BlankPostset = 480e-9;            % blank of transmit after TX pulse
HW.RX.VGAGainDef = HW.RX.VGAGainMax/6;  % default receiver gain

HW.TX.BlankOffsetAQ = 400e-9;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 400e-9;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.LNAGain=10^(((-36.7699)-(-60))/20); % 23.2301 dB gain @ 24377838.756 MHz F=1.0667 dB 50 Ohm (-60 dBm cal)

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 50);  % def transmit voltage
