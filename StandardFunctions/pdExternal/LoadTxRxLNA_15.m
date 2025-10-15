%% Settings for LNA S/N 15

HW.RX.LnaSN = 15;

% Resonant coil 15 mm
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.2/20);
HW.TX.ChannelDef=2;         % default TX rf channel
HW.RX2TXdeadTime=5e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=40e-6;     % 40µs (langsame Pindioden im LNA 12) Totzeit des Empfängers nach dem Senden 
HW.TX.BlankOffset=160e-9;   % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=160e-9;  % Zeitzusatz des Blankingsignals nach dem TX Puls

if 0
  HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverstärkung
  HW.RX.LNAGain=10^(((-54.4)-(-80))/20); % 26.0 dB gain @ 13.0 MHz F=1 dB
else
  HW.RX.VGAGainDef=HW.RX.VGAGainMax/1; % Normale Empfangsverstärkung
  HW.RX.LNAGain=10^(((-54.28)-(-80))/20); % 26.0 dB gain @ 13.0 MHz F=0.6 dB
end
HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=800e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=2000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
