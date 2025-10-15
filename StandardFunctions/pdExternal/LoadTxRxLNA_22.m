% LNA only
HW.TX.ChannelDef=2;         % Def RF Channal TX2
HW.TX2RXdeadTime=5e-6;         % Receiver deadtime after TX pulse >10 µs eg. ringing of Coil ~40 µs

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
HW.RX.LNAGain=10^(((-34.63)-(-60))/20); % 25.37 dB gain @ 8.95 MHz F=1.0 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR
HW.RX.LNAGain=10^(((-35.57)-(-60))/20); % 2xlambda/4 24.43 dB gain @ 8.95 MHz F=1.2 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR
HW.RX.LNAGain=10^(((-37.42)-(-60))/20); % 2xlambda/4 22.58 dB gain @ 8.95 MHz F=1.2 dB (-60 dBm cal pico 50 %FS) LNA 22 MosFET BSN20BKR

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
