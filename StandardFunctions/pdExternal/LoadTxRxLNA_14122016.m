% Resonanter Spule 15mmm
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.2/20);
HW.TX.ChannelDef=2;         % Normaler TX HF Channal
HW.RX2TXdeadTime=5e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=30e-6;      % Totzeit des Empfängers nach dem Senden
HW.TX.BlankOffset=9e-6;     % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3; % Normale Empfangsverstärkung

% HW.RX.LNAGain=10^(23.4/20); % 23.4 dB gain @ 13.1 MHz F=2.5 dB
% HW.RX.LNAGain=10^(22.9/20); % 23.4 dB gain @ 13.1 MHz F=2.5 dB
% HW.RX.LNAGain=10^(24.3/20); % 23.4 dB gain @ 13.1 MHz F=2.5 dB
HW.RX.LNAGain=10^(26.3/20)*990/1414; % 23.4 dB gain @ 13.1 MHz 
HW.RX.LNAGain=10^(26.3/20)*990/1414*925/1000; % 23.4 dB gain @ 13.1 MHz -3db
% HW.RX.LNAGain=1; % 23.4 dB gain @ 13.1 MHz Ohne LNA
% HW.RX.Uin2VgaUin=1.41*1594/1414;%@ 13.1 MHz peak
HW.RX.Uin2VgaUin=1.41*1594/1414;%@ 13.1 MHz eff

HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

