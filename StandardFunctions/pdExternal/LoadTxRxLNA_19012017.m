% LNA 1901017
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.3/20);
HW.TX.ChannelDef=2;         % Normaler TX HF Channal
HW.RX2TXdeadTime=5e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=25e-6;      % Totzeit des Empfängers nach dem Senden
HW.TX.BlankOffset=9e-6;     % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/5; % Normale Empfangsverstärkung

% HW.RX.LNAGain=10^(23.9/20); % 24.8 MHz NF=1.8 dB  (26.07.2017 VNWA 20.16 dB)
HW.RX.LNAGain=10^(19/20); % 24.8 MHz NF=1.8 dB  (26.07.2017 VNWA 20.16 dB)

HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=HW.TX.BlankOffset+1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=HW.TX.BlankPostset+05000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                    % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
