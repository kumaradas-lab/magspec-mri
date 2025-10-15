if ~exist('iDevice', 'var'), iDevice = 1; end

%% rf audio amplifier AmpN2k SN 001
HW.TX(iDevice).ExtRFSN = 1;
HW.TX(iDevice).ExtRFType = 'RFHoellstern';  % string with type identifier

HW.TX(iDevice).ChannelDef = 1;  % Default TX Channel set to Tx1


TX_PaRin = 5e3;  % impetance of rf amplifier at input in Ohm (from datasheet)
TX_UoutRout = 50;  % impetance of TX2 port of MMRT in Ohm
matchFactor = 2 * TX_PaRin/(TX_PaRin+TX_UoutRout);  % factor compared to 50 Ohm termination
HW.TX(iDevice).MmrtUout2Uout(2) = HW.TX(iDevice).MmrtUout2Uout(2) * matchFactor;
HW.TX(iDevice).Uout2PaUout(1) = 1;  % amplification factor (gain)
HW.TX(iDevice).Uout2PaUout(2) = 10;  % amplification factor (gain) (data sheet: 19.5-20.5 dB)

HW.TX(iDevice).Max.Uout(1) = 10;  % maximum peak amplifier input voltage in V
HW.TX(iDevice).Max.Amplitude(1) = Inf;  % maximum peak output voltage in V
% HW.TX(iDevice).Max.NormCalibrated(1) = 2;

HW.TX(iDevice).Max.Uout(2) = 10;  % maximum peak amplifier input voltage in V
HW.TX(iDevice).Max.Amplitude(2) = Inf;  % maximum peak output voltage in V
% HW.TX(iDevice).Max.NormCalibrated(2) = 2;

HW.TX(iDevice).AmplPowerOutCW = 100000;
HW.TX(iDevice).AmplMaxPower = 100000;

% Offset U PowerAmplifer
AUXDAC2_Amp = round(-0.1e-3/3e-3*(2^10-1));  % AmpN2k SN002 - TODO!! 05-07-2024 after 30 min pre-heat
HW.TX.AuxDacOffset(2) =  HW.TX.AuxDacOffset(2) + AUXDAC2_Amp;  % DC offset at TX channels
