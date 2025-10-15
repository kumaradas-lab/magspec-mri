% Load settings for LNA #25 22.813840  MHz
HW.RX.LnaSN = 25;

% HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.4/20); % -0.4 dB Gain @ 23.5 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 1e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 5e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 1e-6;               % blank of transmit before TX pulse
HW.TX.BlankPostset = 1.4e-6;            % blank of transmit after TX pulse
HW.RX.VGAGainDef = HW.RX.VGAGainMax/6;  % default receiver gain

HW.TX.BlankOffsetAQ = 1.4e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 1e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.LNAGain=10^(((-35.599)-(-60))/20); % 24.401 dB gain @ 22814173.6183 MHz F=0.85267 dB (-60 dBm cal) %HW.RX.VGAGainDef = HW.RX.VGAGainMax/6;  % default receiver gain

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 101);  % max transmit voltage
% HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 50);  % def transmit voltage
