% GRADHF für 63 MHz umgebaut (Gradient defekt)
HW.TX.ChannelDef=1;         % Channel 1 kommt auf TX2 raus. Normaler TX HF Channal
HW.TX.Max.Uout=[0.333,4];    % Maximale Ausgangsspannung
HW.TX.Def.Uout=[0.3,3.0];   % Ende des linearen Bereiches der HF Verstärkung
HW.TX.BlankOffset=1e-6;    % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls
HW.TX.MmrtUout2Uout=[1.724,2.5];   % 1.724 bei 63 MHz Verstärkung der Sendeemfangskarte bei 50 Ohm Last
HW.TX.Rout=50;

HW.TX.BlankOffsetAQ=1000e-9;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=1000e-9;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=0;                    % nicht verwenden da kein bestückt. Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.

HW.TX2RXdeadTime=2e-6;      % Totzeit des Empfängers nach dem Senden
HW.RX2TXdeadTime=2e-6;      % Totzeit des Empfängers nach dem Senden

HW.RX.Uin2VgaUin=0.75;    %0.8 bei 63 MHz Volt am Eingang ( 50 Ohm ) zu VGA Eingang (50 Ohm)
HW.RX.Gain2VGAVolts=@AD8331Gain2Volts; % VGAVolts=HW.RXGain2VGAVolts(Gain, HW);
[~,HW.RX.VGAGainMin,HW.RX.VGAGainMax]=AD8331Gain2Volts(1,HW);
HW.RX.Rin=50;
HW.RX.VgaUout2MmrtUin=0.909;   %20 Rout /200 Rin ohm

%LNA
HW.RX.VGAGainDef=HW.RX.VGAGainMax/5; % Normale Emfangsverstärkung
