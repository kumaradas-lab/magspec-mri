% Load settings for LNA #0 PHA-13HLN+ 22.813840  MHz
HW.RX.LnaSN = 0;

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

% HW.RX.LNAGain=10^(((-33.564)-(-60))/20); % 26.436 dB gain @ 22814201.9968 MHz F=1.63 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-31.5957)-(-60))/20); % 28.4043 dB gain @ 5200082.2538 MHz F=1.824 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-31.5925)-(-60))/20); % 28.4075 dB gain @ 1200019.3851 MHz F=4.1544 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-31.299)-(-60))/20); % 28.701 dB gain @ 2200033.9818 MHz F=2.951 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-31.3398)-(-60))/20); % 28.6602 dB gain @ 3600055.3662 MHz F=2.076 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-33.7134)-(-60))/20); % 26.2866 dB gain @ 30814251.397 MHz F=1.7331 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-33.5526)-(-60))/20); % 26.4474 dB gain @ 22814189.4836 MHz F=1.5213 dB (-60 dBm cal)
% 
% HW.RX.LNAGain=10^(((-41.364)-(-60))/20); % 18.636 dB gain @ 22814229.7189 MHz F=1.6449 dB (-60 dBm cal) 22.5R bias

%3V
% HW.RX.LNAGain=10^(((-36.9297)-(-60))/20); % 23.0703 dB gain @ 22902123.1764 MHz F=1.2427 dB (-60 dBm cal)
% HW.RX.LNAGain=10^(((-38.6551)-(-60))/20); % 21.3449 dB gain @ 22902189.6138 MHz F=1.8196 dB (-60 dBm cal) mit RF-100
HW.RX.LNAGain=10^(((-35.9059)-(-60))/20); % 24.0941 dB gain @ 8065135.296 MHz F=0.92167 dB (-60 dBm cal)


HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 50);  % def transmit voltage
