%% Settings for LNA S/N 10

HW.RX.LnaSN = 10;

% Resonant coil 15 mm
% HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.2/20);
HW.TX.ChannelDef=2;         % default TX rf channel
% HW.RX2TXdeadTime=5e-6;      % Totzeit des Empfängers vor dem Senden
% HW.TX2RXdeadTime=3e-6;      % 1µs -0.2 dB 5µs ok (schnelle Pindioden im LNA 10) Totzeit des Empfängers nach dem Senden 
% HW.TX.BlankOffset=160e-9;     % Offset des Blankingsignals vor dem TX Puls
% HW.TX.BlankPostset=160e-9;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverstärkung

% HW.RX.LNAGain=10^((22.7-1.2)/20); % 22.7 dB gain @ 24.71 MHz F=1 dB
HW.RX.LNAGain=10^(((-37.8286)-(-60))/20); % 22.1714 dB gain @ 24790905.9089 MHz F=0.65698 dB (-60 dBm cal) (HW.RX.VGAGainDef=HW.RX.VGAGainMax/3)

HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

% %% TRx switch during transmit at Tx2
% HW.TX.BlankOffsetAQ=80e-9;        % Blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ=240e-9;       % Blank of receiver after TX pulse
% HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.

% %% TRx switch during transmit at Tx2
% HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ=2000e-9;       % Blank of receiver after TX pulse
% HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.