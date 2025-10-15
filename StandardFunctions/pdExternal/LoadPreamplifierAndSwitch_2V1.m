% Load this file in LoadMySystem.m if you use the external Preamplifier & Switch

HW.TX.Uout2PaUout=[1,10^(-0.4/20)]; % -0.4 dB Gain @ 23.5 MHz
HW.TX.ChannelDef=2;                 % Default TX HF Channel
HW.RX2TXdeadTime=3e-6;              % dead-time of receiver before TX pulse
HW.TX2RXdeadTime=13e-6;             % dead-time of receiver after TX pulse
HW.TX.BlankOffset=2e-6;             % Blank of transmit before TX pulse
HW.TX.BlankPostset=480e-9;          % Blank of transmit after TX pulse
HW.RX.VGAGainDef=HW.RX.VGAGainMax/1;% Default receiver gain

HW.TX.BlankOffsetAQ=2000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.

HW.RX.LNAGain=10^(21.3/20);         % 21.3 dB LNA receive gain @ 23.5 MHz  

HW.TX.Max.PaUout(2)=100;               % Max transmit voltage
HW.TX.Def.PaUout(2)=100;               % def transmit voltage

