%% settings for LNA in Aspect system @ approx. 42 MHz
HW.RX(1).LnaSN = 10172;  % dummy serial number

HW.TX(1).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(1).BlankOffset = 1e-6;  % blank of transmit before TX pulse
HW.TX(1).BlankPostset = 1e-6;  % blank of transmit after TX pulse

HW.TX(1).BlankOffsetAQ = 1e-6;  % blank of receiver before TX pulse
HW.TX(1).BlankPostsetAQ = 20e-6;  % blank of receiver after TX pulse
HW.TX(1).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.


HW.RX(1).VGAGainDef = HW.RX(1).VGAGainMax/6;  % default receiver gain range divisor (3 to 6)
HW.RX(1).LNAGain = 10^(11.7/20);  % 30 dB gain - estimated value


% HW.TX.PaUout2Amplitude(2) = 0.00004;  % 5 µs estimated pulse @ 60 V
