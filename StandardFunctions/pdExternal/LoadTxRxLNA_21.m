%% Settings for LNA S/N 21

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 21;

% LNA only (no switch)
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.TX2RXdeadTime = 30e-6;  % receiver dead-time after TX pulse >10 us, e.g., ringing of coil ~40 us

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % reduce VGA gain to avoid saturation
HW.RX(iDevice).LNAGain = 10^(((-37.41)-(-60))/20);  % 22.59 dB gain @ 21.7 MHz F=0.8 dB (-60 dBm cal pico 50 %FS) LNA 21 MosFET BSN20BKR

%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 1000e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 4000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
