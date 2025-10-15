% Load settings for LNA #30 2-40 MHz
HW.RX.LnaSN = 30;

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
HW.RX.LNAGain=10^(((-37.8162)-(-60))/20); % 22.1838 dB gain @ 2300025.7339 MHz F=1.9949 dB (-60 dBm cal)
HW.RX.LNAGain=10^(((-37.6362)-(-60))/20); % 22.3638 dB gain @ 5300059.0745 MHz F=1.3271 dB (-60 dBm cal)
HW.RX.LNAGain=10^(((-38.4329)-(-60))/20); % 21.5671 dB gain @ 12300136.7997 MHz F=1.4805 dB (-60 dBm cal)
HW.RX.LNAGain=10^(((-40.1127)-(-60))/20); % 19.8873 dB gain @ 24300269.5874 MHz F=1.4805 dB (-60 dBm cal)

 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 10);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 4);  % def transmit voltage