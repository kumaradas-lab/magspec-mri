% Load settings for LNA #19 23.5 MHz
HW.RX.LnaSN = 19;

% HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.4/20); % -0.4 dB Gain @ 23.5 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 1e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 4e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 1e-6;               % blank of transmit before TX pulse
HW.TX.BlankPostset = 480e-9;            % blank of transmit after TX pulse

if 0
  HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
  % HW.RX.LNAGain=10^(((-35.8684)-(-60))/20); % 24.1316 dB gain @ 23550303.9158 MHz F=0.85376 dB (-60 dBm cal) (without RF-100)
  HW.RX.LNAGain=10^(((-36.8717)-(-60))/20); % 23.1283 dB gain @ 23500381.9362 MHz F=2.0899 dB (-60 dBm cal)(with RF-100)
else
  HW.RX.VGAGainDef=HW.RX.VGAGainMax/6;    % reduce VGA gain to avoid saturation
  HW.RX.LNAGain=10^(((-36.9085)-(-60))/20); % 23.0915 dB gain @ 23500382.707 MHz F=2.2518 dB (-60 dBm cal)
end

HW.TX.BlankOffsetAQ = 400e-9;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 400e-9;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 50);  % def transmit voltage
