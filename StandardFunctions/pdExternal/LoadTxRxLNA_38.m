% Load settings for LNA #38 22.9 MHz und 24.3 MHz
HW.RX.LnaSN = 38;

HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.3/20); % -0.3 dB Gain (TX to Coil) @ 24.15 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 9e-6; % >7µs         % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain

% HW.RX.LNAGain=10^(((-44.9731)-(-67.633))/20); % 22.6599 dB gain @ 22900331.5471 MHz F=0.78729 dB (-67.633 dBm cal)
% HW.RX.LNAGain=10^(((-45.2458)-(-67.633))/20); % 22.3872 dB gain @ 24300351.3135 MHz F=0.71811 dB (-67.633 dBm cal)
HW.RX.LNAGain=10^(((-44.8925)-(-67.633))/20); % 22.7405 dB gain @ 22900330.2545 MHz F=0.66069 dB (-67.633 dBm cal)
HW.RX.LNAGain=10^(((-45.1655)-(-67.633))/20); % 22.4675 dB gain @ 24300350.7268 MHz F=0.70906 dB (-67.633 dBm cal)
 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 72);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 50);  % def transmit voltage