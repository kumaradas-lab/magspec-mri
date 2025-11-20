%% Settings for LNA S/N 17

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 17;

% LNA only (no switch)
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel

% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/10;  % reduce VGA gain to avoid saturation
% HW.RX(iDevice).LNAGain = 10^(((-55.93)-(-80))/20);  % 21.55 dB gain @ 24.0 MHz F=0.9 dB

% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % reduce VGA gain to avoid saturation
% HW.RX(iDevice).LNAGain = 10^(((-56.08)-(-80))/20);  % 21.55 dB gain @ 24.0 MHz F=0.9 dB


% HW.RX(iDevice).LNAGain = 10^(((-33.60)-(-60))/20);  % 26.40 dB gain @ 8.95 MHz F=0.8 dB ; groﬂer C am Eingang test
% HW.RX(iDevice).LNAGain = 10^(((-30.97-5.76+1.52)-(-60))/20);  % 26.40 dB gain @ 8.95 MHz F=0.8 dB ;
HW.TX2RXdeadTime = 5e-6;  % receiver deadtime after TX pulse, e.g., ringing of TRx coil ~40 us

% HW.RX(iDevice).LNAGain = 10^(((-39.89)-(-60))/20);  % 20.02 dB gain @ 8.95 MHz F=0.8 dB ;
% HW.RX(iDevice).LNAGain = 10^(((-39.49)-(-60))/20);  % 20.02 dB gain @ 8.95 MHz F=2.5 dB ; 20µH 220 Ohm
HW.RX(iDevice).LNAGain = 10^(((-40.05)-(-60))/20);  % 19.95 dB gain @ 8.95 MHz F=2.8 dB ; 2x20µH 220 Ohm

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % reduce VGA gain to avoid saturation

% %% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 1000e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
