%% Settings for LNA S/N 21

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 21;

% LNA only (no switch)
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.TX2RXdeadTime = 5e-6;  % receiver dead-time after TX pulse >10 us, e.g., ringing of Coil ~40 us

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % reduce VGA gain to avoid saturation
HW.RX(iDevice).LNAGain = 10^(((-34.63)-(-60))/20);  % 25.37 dB gain @ 8.95 MHz F=1.0 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR
HW.RX(iDevice).LNAGain = 10^(((-35.57)-(-60))/20);  % 2xlambda/4 24.43 dB gain @ 8.95 MHz F=1.2 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR
HW.RX(iDevice).LNAGain = 10^(((-37.42)-(-60))/20);  % 2xlambda/4 22.58 dB gain @ 8.95 MHz F=1.2 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR

%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 1000e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 4000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
