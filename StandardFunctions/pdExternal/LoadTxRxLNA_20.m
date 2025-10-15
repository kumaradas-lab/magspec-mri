%% Settings for LNA S/N 20

HW.RX.LnaSN = 20;

% LNA only
HW.TX.ChannelDef=2;         % default TX rf channel
HW.TX2RXdeadTime=20e-6;     % Receiver deadtime after TX pulse >10 µs eg. ringing of Coil ~40 µs

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
HW.RX.LNAGain=10^(((-39.113)-(-60))/20); % 19.93 dB gain @ 45 MHz F=1.0 dB (-60 dBm cal pico 50 %FS) LNA20 MosFET BSN20BKR

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
