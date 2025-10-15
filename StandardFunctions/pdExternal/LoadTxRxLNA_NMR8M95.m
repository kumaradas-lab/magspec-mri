% LNA only
HW.TX.ChannelDef=2;         % Def RF Channal TX2
HW.TX2RXdeadTime=4e-6;         % Receiver deadtime after TX pulse eg. ringing of Coil ~40 µs

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
HW.RX.LNAGain=10^(((-57.55)-(-80))/20); % 22.45 dB gain @8.95 MHz F=0.6 dB (-80 dBm cal)

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=1400e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
