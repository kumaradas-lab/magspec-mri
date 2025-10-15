%% LoadTxRxLNA_06102015 
% Add filename to LoadMySystem.m if you use the TX/RX switch with LNA

% connected to high Q resonant coil with 10mmm
HW.TX.Uout2PaUout=[1,10^(-0.1/20)];     % linear Gain of TX channel 1 and 2
HW.TX.ChannelDef=2;                     % normal TX channel (1 => TRx, 2 => TX2)
HW.RX2TXdeadTime=1e-6;                  % dead time between receive and transmit
HW.TX2RXdeadTime=10e-6;                 % dead time between transmit and receive (you can reduce this value until you see pulse artefacts in the AQ)
HW.TX.BlankOffset=2e-6;                 % time from Blk and Sw changes to TX pulse
HW.TX.BlankPostset=0.32e-6;             % time from TX pulse to Blk and Sw changes  
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % Normale Empfangsverstärkung
 
% HW.RX.LNAGain=10^(23.8/20); % Gain of preamp at 24 MHz
% HW.RX.LNAGain=10^(14.3/20); % Gain of preamp  at 22 MHz 
% HW.RX.LNAGain=10^(19.6/20); % Gain of preamp  at 4.7 MHz MAR-6SM+
% HW.RX.LNAGain=10^(19.75/20); % Gain of preamp  at 9 MHz MAR-6SM+

% HW.RX.LNAGain=10^(27.3/20); % Gain of preamp  at 9 MHz PSA4-5043+ 5V 
% HW.RX.LNAGain=10^(28/20); % Gain of preamp  at 4.7 MHz PSA4-5043+ 5V
% HW.RX.LNAGain=10^(25.7/20); % Gain of preamp  at 9 MHz PSA4-5043+ 3V
HW.RX.LNAGain=10^(26.4/20); % Gain of preamp  at 4.7 MHz PSA4-5043+ 3V NF=2.0 dB
HW.RX.LNAGain=10^(25.4/20); % Gain of preamp  at 8.95 MHz PSA4-5043+ 3V NF<=1 dB
 
HW.TX.Max.PaUout(2)=100; % limit of input amplitude of the switch 
HW.TX.Def.PaUout(2)=100; % default input amplitude of the switch
 
%% TRx switch during transmit at Tx2 to prevent saturation effects of TRx
HW.TX.BlankOffsetAQ=HW.TX.BlankOffset+1000e-9;       % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=HW.TX.BlankPostset+3000e-9;     % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                                    % Switch TRx to 50 Ohm Resistor during TX2 pulse, to avoid saturation of TRx. (0 => do not switch TRx, 1 => Switch TRx to 50 Ohm during TX2 pulse)
