%% Settings for LNA S/N 16

HW.RX.LnaSN = 16;

% LNA only
HW.TX.ChannelDef=2;         % default TX rf channel

HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % reduce VGA gain to avoid saturation
HW.RX.LNAGain=10^(((-58.45)-(-80))/20); % 21.55 dB gain @ 24.0 MHz F=0.9 dB

% %% TRx switch during transmit at Tx2
% HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
% HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
