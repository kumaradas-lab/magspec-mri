% Resonanter Spule 11mmm?
HW.TX.Uout2PaUout=[1,10^(-0.1/20)];
HW.TX.ChannelDef=2;         % Normaler TX HF Channal
HW.RX2TXdeadTime=1e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=9e-6;      % Totzeit des Empfängers nach dem Senden
HW.TX.BlankOffset=9e-6;     % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverstärkung

HW.RX.LNAGain=10^(23.8/20); % 23,7 MHz
HW.RX.LNAGain=10^(23.3/20); % 23,7 MHz

HW.TX.Max.PaUout(2)=100;
HW.TX.Def.PaUout(2)=100;

