%% Settings for LNA S/N 13

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 13;

% Resonant coil 10 mm
HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.TX(iDevice).Def.Uout(2) = 3.8;  % default TX Uout in V

HW.RX2TXdeadTime = 4e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % 50 us (slow PIN diodes in LNA #12) - dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 3e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 1e-6;  % blank of transmit after TX pulse

switch 4
  case 1
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1;  % max receiver gain -0 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.34)-(-80))/20);  % 24.66 dB gain @ 24.35 MHz F=1.3 dB
  case 4
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/4;  % max receiver gain -12 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.39)-(-80))/20);  % 24.61 dB gain @ 24.35 MHz F=1.4 dB
  case 10
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/10;  % max receiver gain -20 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.36)-(-80))/20);  % 24.64 dB gain @ 24.35 MHz F=1.6 dB
  case 20
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/20;  % max receiver gain -20 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.30)-(-80))/20);  % 24.70 dB gain @ 24.35 MHz F=2.2 dB
  case 40
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/40;  % max receiver gain -20 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.27)-(-80))/20);  % 24.73 dB gain @ 24.35 MHz F=4 dB
  case 100
    HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/100;  % max receiver gain -20 dB
    HW.RX(iDevice).LNAGain = 10^(((-55.46)-(-80))/20);  % 24.54 dB gain @ 24.35 MHz F=9.4 dB
end
HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);


%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 1e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 4e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
