%% Settings for LNA S/N 17

HW.RX.LnaSN = 17;

% LNA only
HW.TX.ChannelDef=2;         % default TX rf channel

% HW.RX.VGAGainDef=HW.RX.VGAGainMax/10;    % reduce VGA gain to avoid saturation
% HW.RX.LNAGain=10^(((-55.93)-(-80))/20); % 21.55 dB gain @ 24.0 MHz F=0.9 dB

% HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
% HW.RX.LNAGain=10^(((-56.08)-(-80))/20); % 21.55 dB gain @ 24.0 MHz F=0.9 dB


% HW.RX.LNAGain=10^(((-33.60)-(-60))/20); % 26.40 dB gain @ 8.95 MHz F=0.8 dB ; groﬂer C am Eingang test
% HW.RX.LNAGain=10^(((-30.97-5.76+1.52)-(-60))/20); % 26.40 dB gain @ 8.95 MHz F=0.8 dB ; 
HW.TX2RXdeadTime=5e-6;         % Receiver deadtime after TX pulse eg. ringing of Coil ~40 µs

% HW.RX.LNAGain=10^(((-39.89)-(-60))/20); % 20.02 dB gain @ 8.95 MHz F=0.8 dB ; 
% HW.RX.LNAGain=10^(((-39.49)-(-60))/20); % 20.02 dB gain @ 8.95 MHz F=2.5 dB ; 20µH 220 Ohm
HW.RX.LNAGain=10^(((-40.05)-(-60))/20);   % 19.95 dB gain @ 8.95 MHz F=2.8 dB ; 2x20µH 220 Ohm

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation

% %% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=5000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
