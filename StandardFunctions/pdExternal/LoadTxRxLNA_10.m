%% Settings for LNA S/N 10

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 10;

% Resonant coil 15 mm
% HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
% HW.RX2TXdeadTime = 5e-6;  % dead-time of receiver before TX pulse
% HW.TX2RXdeadTime = 3e-6;  % 1 us -0.2 dB 5us ok (fast PIN diodes in LNA #10) - dead-time of receiver after TX pulse
% HW.TX(iDevice).BlankOffset = 160e-9;  % blank of transmit before TX pulse
% HW.TX(iDevice).BlankPostset = 160e-9;  % blank of transmit after TX pulse

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

% HW.RX(iDevice).LNAGain = 10^((22.7-1.2)/20);  % 22.7 dB gain @ 24.71 MHz F=1 dB
HW.RX(iDevice).LNAGain = 10^(((-37.8286)-(-60))/20);  % 22.1714 dB gain @ 24790905.9089 MHz F=0.65698 dB (-60 dBm cal) (HW.RX.VGAGainDef=HW.RX.VGAGainMax/3)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);

% %% TRx switch during transmit at Tx2
% HW.TX(iDevice).BlankOffsetAQ = 80e-9;  % blank of receiver before TX pulse
% HW.TX(iDevice).BlankPostsetAQ = 240e-9;  % blank of receiver after TX pulse
% HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

% %% TRx switch during transmit at Tx2
% HW.TX(iDevice).BlankOffsetAQ = 1000e-9;  % blank of receiver before TX pulse
% HW.TX(iDevice).BlankPostsetAQ = 2000e-9;  % blank of receiver after TX pulse
% HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
% if ~HW.TX(iDevice).ExtRFSN && isfield(HW.TX(iDevice).CalibrationRfAmp, 'TriScatteredAmp')
%   HW.TX(iDevice).CalibrationRfAmp = rmfield(HW.TX(iDevice).CalibrationRfAmp,'TriScatteredAmp');
% end
if ~HW.TX(iDevice).ExtRFSN
  HW.TX(iDevice).CalibrationRfAmp = [];
end
