%% Settings for LNA S/N 21

HW.RX.LnaSN = 21;

% LNA only
HW.TX.ChannelDef=2;         % default TX rf channel
HW.TX2RXdeadTime=30e-6;         % Receiver deadtime after TX pulse >10 µs eg. ringing of Coil ~40 µs

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
HW.RX.LNAGain=10^(((-37.41)-(-60))/20); % 22.59 dB gain @ 21.7 MHz F=0.8 dB (-60 dBm cal pico 50 %FS) LNA 21 MosFET BSN20BKR

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
