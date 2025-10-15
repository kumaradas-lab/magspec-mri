%% Settings for TxRx switch
% Must be loaded after LoadRFAmp and LoadLNA files.

%% attenuation of external TxRx switch
HW.TX.Uout2PaUout(2) = HW.TX.Uout2PaUout(2)*10^(-0.2/20); % switch with TX attenuation of -0.2 dB
HW.RX.LNAGain = HW.RX.LNAGain*10^(-0.0/20);               % switch with RX attenuation of -0.0 dB
% HW.TX.Max.PaUout(2) = min(HW.TX.Max.PaUout(2), 100);    % maximum input voltage of external TRx switch
