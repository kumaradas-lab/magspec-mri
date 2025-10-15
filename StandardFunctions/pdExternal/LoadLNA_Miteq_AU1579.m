%% LNA Miteq AU-1579
% This file overwrites some settings from LoadRF file.
% Must be loaded after LoadRF file and before LoadTxRx file.

HW.TX.ChannelDef = 2;                   % default channel: TX2
HW.RX2TXdeadTime = 1e-6;                % deadtime of receiver before transmission (after acquisition)
HW.TX2RXdeadTime = 50e-6;               % deadtime of receiver after transmission
HW.TX.BlankOffset = 160e-9;             % offset of Blank before TX pulse
HW.TX.BlankPostset = 160e-9;            % delay of Blank after TX pulse
HW.RX.VGAGainDef = HW.RX.VGAGainMax/10; % default gain of variable gain amplifier in drive-l

HW.RX.LNAGain = 10^(37/20); % 37 dB gain @ 1-100 MHz with noise figure F=1 dB

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 800e-9;           % Blank of receiver before TX pulse in seconds
HW.TX.BlankPostsetAQ = 2000e-9;         % Blank of receiver after TX pulse in seconds
HW.TX.BlankAQ = 1;                      % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
