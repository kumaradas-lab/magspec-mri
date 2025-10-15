% Load settings for LNA #36
HW.RX.LnaSN = 36;

% calibrated with RF-100HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX.ChannelDef = 2;                   % default TX rf channel
HW.RX2TXdeadTime = 10e-6;                % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;               % dead-time of receiver after TX pulse
HW.TX.BlankOffset = 9e-6; % >7µs         % blank of transmit before TX pulse
HW.TX.BlankPostset = 5e-6;            % blank of transmit after TX pulse

HW.TX.BlankOffsetAQ = 2e-6;          % blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 5e-6;         % blank of receiver after TX pulse
HW.TX.BlankAQ = 1;                      % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain
HW.RX.LNAGain=10^(((-42.2888)-(-67.6))/20); % 25.3112 dB gain @ 24500352.6402 MHz F=1.065 dB (-67.6 dBm cal)
 
HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 101);  % max transmit voltage
HW.TX.Def.PaUout(2) = min(HW.TX.Def.PaUout(2), 30);  % def transmit voltage

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
