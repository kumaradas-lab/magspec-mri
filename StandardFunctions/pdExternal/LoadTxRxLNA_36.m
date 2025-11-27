%% Settings for LNA #36

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 36;

% calibrated with RF-100 HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain
HW.RX(iDevice).LNAGain = 10^(((-42.2888)-(-67.6))/20);  % 25.3112 dB gain @ 24500352.6402 MHz F=1.065 dB (-67.6 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 101);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 30);  % def transmit voltage

switch HW.UserName
  case 'magnet_01_probe_1H'
    % Magnet 1: 1H Wechselprobenkopf
    % Magnet 165
  case 'magnet_01_probe_31P'
    % Magnet 1: 31P Wechselprobenkopf
    % Magnet 165
    warning('LNA frequency range: 24-25 MHz; 31P probe selected @ 9.8 MHz');
  case 'magnet_02'
    % Magnet 2: 10mm Magnet KEIN Wechselprobenkopf
    % Magnet 166
end
