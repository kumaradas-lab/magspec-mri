%% Settings for rf amplifier RF-200 - no switch included

HW.TX.ExtRFSN = 9;
HW.TX.ExtRFType = 'RF200';  % string with type identifier

HW.TX.ChannelDef = 2;         % Default TX Channel set to Tx2

% Amplifier
HW.TX.Uout2PaUout(2) = 10^((-15.0+53.0-3.0)/20);  % 53.0 dB amplification @ 10.0 MHz, 1.0 V without switch
% HW.TX.Uout2PaUout(2) = 10^((-15.0+52.2-3.0)/20);  % 52.2 dB amplification @ 24.0 MHz, 1.0 V with Switch

HW.TX.Max.PaUout(2) = 250;  % maximum peak output voltage in V
HW.TX.Def.PaUout(2) = 80;  % default PaUout(2)=50 V => 25 W

HW.TX.Max.Amplitude = [20,20]*1e-3;   % maximum B1+ in T
HW.TX.Def.Amplitude = [20,20]*1e-3;   % default B1+ in T

%% no switch included
UseExtRFAmpSwitch = 2;  % 2: use external switch in probe (not yet calibrated)

LoadExtRFAmp_Cal;  % new cal Uout and 6 A FET

%% Coil 15mm (long)
HW.RX2TXdeadTime = 2.4e-6;  % receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 10e-6;  % receiver deadtime after TX pulse in s, e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 2000e-9;  % unblank of RF-200 before TX pulse
HW.TX.BlankPostset = 1e-6;  % blank of RF-200 after TX pulse

% rf coil damping Settings
HW.TX.DampCoil.Enable = true;  % enable coil damping
HW.TX.DampCoil.DigitalOutputLatency = 0.6e-6;  % latency of damping circuit in seconds
HW.TX.DampCoil.DigitalOutputChannel = 1;  % digital output channel for coil damping signal

HW.TX.DampCoil.DigitalOutputDuration = HW.TX.BlankPostset + 1e-6;  % duration of signal on digital out in seconds
HW.TX.DampCoil.DampingDuration = HW.TX.DampCoil.DigitalOutputDuration + 2e-6;  % effective duration of damping the coil in seconds

HW.TX.DampCoil.TX2RXdeadTime = HW.TX.DampCoil.DampingDuration + 3e-6;  % dead time between pulse and acquisition in seconds with enabled damping (additional time: ~3 * Q/pi/f0)

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 2000e-9;  % Blank of receiver before TX pulse
if HW.TX.DampCoil.Enable
  HW.TX.BlankPostsetAQ = HW.TX.DampCoil.DigitalOutputDuration + 1.0e-6;  % Blank of receiver after TX pulse
else
  HW.TX.BlankPostsetAQ = 6000e-9;
end
HW.TX.BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX settings
% HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % reduce VGA gain to avoid saturation

% HW.RX.LNAGain = 10^(((-50.2719)-(-71.14))/20);  % 20.8681 dB gain @ 22400423.3518 MHz F=1.2861 dB (-71.14 dBm cal)
% LoadTxRxLNA_65;
