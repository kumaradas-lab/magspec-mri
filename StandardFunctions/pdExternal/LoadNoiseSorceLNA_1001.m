%% Settings for LoadNoiseSorceLNA_1001

HW.RX.LnaSN = 1001;

% LNA only
HW.TX.ChannelDef=2;         % default TX rf channel
HW.TX2RXdeadTime=30e-6;         % Receiver deadtime after TX pulse eg. ringing of Coil ~40 µs

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
% HW.RX.LNAGain=10^(((-55.78+0.32)-(-80))/20); % 24.22 dB gain @ 24.799470 MHz F=1.3 dB (-80 dBm cal)
% HW.RX.LNAGain=10^(((-35.19)-(-60))/20); % 24.81 dB gain @ 24.799470 MHz F=1.3 dB (-80 dBm cal)
HW.RX.LNAGain=10^(((-42.7992)-(-60))/20); % 17.2008 dB gain @ 13100201.351 MHz F=3.1909 dB (-60 dBm cal)


%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=0;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
