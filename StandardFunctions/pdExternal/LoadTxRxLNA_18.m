%% Settings for LNA S/N 18

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 18;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB gain (TX to coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;   % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain
% HW.RX(iDevice).LNAGain = 10^(((-35.0552)-(-60))/20);  % 24.9448 dB gain @ 24300413.2391 MHz F=0.82022 dB (-60 dBm cal)

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % default receiver gain
HW.RX(iDevice).LNAGain = 10^(((-34.9642)-(-60))/20);  % 25.0358 dB gain @ 24300410.4328 MHz F=0.9 dB (-60 dBm cal)

% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/10;  % default receiver gain
% HW.RX(iDevice).LNAGain = 10^(((-34.944)-(-60))/20);  % 25.056 dB gain @ 24300413.9119 MHz F=1.0946 dB (-60 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 10);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 4);  % def transmit voltage
