% Load settings for LNA #23
HW.RX.LnaSN = 23;

HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.4/20); % -0.4 dB Gain @ 23.5 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 3e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 13e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 2e-6;               % blank of transmit before TX pulse
HW.TX.BlankPostset = 480e-9;            % blank of transmit after TX pulse
HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain

HW.TX.BlankOffsetAQ = 2000e-9;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 4000e-9;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.LNAGain = 10^(21.3/20);           % 21.3 dB LNA receive gain @ 23.5 MHz

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 100);  % def transmit voltage
