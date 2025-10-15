%% Settings for LNA S/N 18

HW.RX.LnaSN = 18;


HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 9e-6; % >7µs         % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

% HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain
% HW.RX.LNAGain=10^(((-35.0552)-(-60))/20); % 24.9448 dB gain @ 24300413.2391 MHz F=0.82022 dB (-60 dBm cal)

HW.RX.VGAGainDef = HW.RX.VGAGainMax/6;  % default receiver gain
HW.RX.LNAGain=10^(((-34.9642)-(-60))/20); % 25.0358 dB gain @ 24300410.4328 MHz F=0.9 dB (-60 dBm cal)

% HW.RX.VGAGainDef = HW.RX.VGAGainMax/10;  % default receiver gain
% HW.RX.LNAGain=10^(((-34.944)-(-60))/20); % 25.056 dB gain @ 24300413.9119 MHz F=1.0946 dB (-60 dBm cal)

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 10);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 4);  % def transmit voltage
