if ~exist('iDevice', 'var'), iDevice = 1; end

%% (gradient) amplifier DC-600 used on rf channels
HW.TX(iDevice).ExtRFSN = 41;
HW.TX(iDevice).ExtRFType = 'RFDC600';  % string with type identifier

HW.TX(iDevice).ChannelDef = 2;  % default TX channel

HW.TX(iDevice).Uout2PaUout(1) = 1;  % active gradiometer
TX_UoutRout = 50;  % impedance of TX2 port of MMRT in Ohm
TX_PaRin = 5e3;  % impedance of differential converter at input in Ohm
matchFactor = 2 * TX_PaRin/(TX_PaRin+TX_UoutRout);  % factor compared to 50 Ohm termination
HW.TX(iDevice).MmrtUout2Uout(2) = HW.TX(iDevice).MmrtUout2Uout(2) * matchFactor;

HW.TX(iDevice).Max.Uout(1:2) = 15/2/2*0.81;  % maximum peak amplifier input voltage in V
HW.TX(iDevice).Max.Amplitude(1:2) = Inf;  % maximum peak output voltage in V

