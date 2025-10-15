%% LoadTxRxLNA_63_4MHz_11112016
% Add filename to LoadMySystem.m if you use the TX/RX switch with LNA

% connected to high Q resonant coil with 10mmm
HW.TX.Uout2PaUout=[1,10^(-0.1/20)];     % linear Gain of TX channel 1 and 2
HW.TX.ChannelDef=2;                     % normal TX channel (1 => TRx, 2 => TX2)
HW.RX2TXdeadTime=1e-6;                  % dead time between receive and transmit
HW.TX2RXdeadTime=10e-6;                 % dead time between transmit and receive (you can reduce this value until you see pulse artefacts in the AQ)
HW.TX.BlankOffset=2e-6;                 % time from Blk and Sw changes to TX pulse
HW.TX.BlankPostset=0.32e-6;             % time from TX pulse to Blk and Sw changes  
HW.RX.VGAGainDef=HW.RX.VGAGainMax/3;    % Normale Empfangsverstärkung
 
HW.RX.LNAGain=10^(15.6/20); % Gain of preamp  at 63.4 MHz MAR-6SM+

HW.TX.Max.PaUout(2)=100; % limit of input amplitude of the switch 
HW.TX.Def.PaUout(2)=100; % default input amplitude of the switch
 
%% TRx switch during transmit at Tx2 to prevent saturation effects of TRx
HW.TX.BlankOffsetAQ=HW.TX.BlankOffset+1000e-9;       % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=HW.TX.BlankPostset+3000e-9;     % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                                    % Switch TRx to 50 Ohm Resistor during TX2 pulse, to avoid saturation of TRx. (0 => do not switch TRx, 1 => Switch TRx to 50 Ohm during TX2 pulse)
