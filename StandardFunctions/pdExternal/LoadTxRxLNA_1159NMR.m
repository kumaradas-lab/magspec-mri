% Resonanter Spule 15mmm
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.3/20);
HW.TX.ChannelDef=2;         % Normaler TX HF Channal
HW.RX2TXdeadTime=5e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=40e-6;      % Totzeit des Empfängers nach dem Senden
HW.TX.BlankOffset=9e-6;     % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverstärkung

HW.RX.LNAGain=22.5*10^(1.8/20); % ? dB gain @ 13.1 MHz F=0.5 dB?

HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

