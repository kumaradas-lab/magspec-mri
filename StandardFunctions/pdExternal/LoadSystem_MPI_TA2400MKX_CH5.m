if ~exist('iDevice', 'var'), iDevice = 1; end

%% rf audio amplifier TA2400 MK-X
HW.TX(iDevice).ExtRFSN = 0;
HW.TX(iDevice).ExtRFType = 'TA2400MKX';  % string with type identifier

HW.TX(iDevice).ChannelDef = 1;  % Default TX Channel set to Tx1

% factor 2 because without 50 Ohm termination
HW.TX(iDevice).Uout2PaUout(1:2) = 91.2*2;  % amplification factor (gain)

HW.TX(iDevice).Max.Uout(1:2) = 5/2;  % maximum peak amplifier input voltage in V
HW.TX(iDevice).Max.Amplitude(1:2) = Inf;  % maximum peak output voltage in V


%% rf coil (channel 5)
HW.TX(iDevice).AmplitudeUnit = 'V';
HW.TX(iDevice).AmplitudeName = '5:sol';
HW.TX(iDevice).PaUout2Amplitude(1:2) = 1;  % output in Volts  - FIXME: Is this correct for a resonant coil?
HW.TX(iDevice).AmplitudeUnitScale = 1;
