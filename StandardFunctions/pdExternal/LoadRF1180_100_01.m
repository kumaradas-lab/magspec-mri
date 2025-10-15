%% Settings for rf amplifier RFPA RF1180-100 for PrimoGAIA

HW.TX(iDevice).ExtRFSN = 218826;
%HW.TX(2).ExtRFSN = 218826;
HW.TX(iDevice).ChannelDef = 2;  % default TX channel set to Tx2

% Amplifier
HW.TX(iDevice).Uout2PaUout(2) = 10^((51-12-6)/20);  % 56.9dB amplification @ 80 MHz, minimum attenuation at input: 17 dB
%HW.TX(iDevice).Uout2PaUout(2) = 10^((48-12-6)/20);  % 56.9dB amplification @ 150 MHz, minimum attenuation at input: 17 dB
%HW.TX(iDevice).Uout2PaUout(2) = 10^((52.3-12-6-5)/20);  % 56.9dB amplification @ 170 MHz, minimum attenuation at input: 12 dB

HW.TX(iDevice).Max.PaUout(2) = sqrt(2*101*50);  % maximum output voltage (peak amplitude), P=100 W, R=50 Ohm, (P=U^2/R)
HW.TX(iDevice).Def.PaUout(2) = sqrt(2*10*50);   % default output voltage (peak amplitude), P=10 W, R=50 Ohm, (P=U^2/R)

% HW.TX(iDevice).Max.Amplitude(2) = 20e-3;   % maximum output in B1 amplitude
% HW.TX(iDevice).Def.Amplitude(2) = 20e-3;   % default output in B1 amplitude


%% Coil 10mm
% HW.RX2TXdeadTime = 1e-6;   % Receiver deadtime before TX pulse in s
% HW.TX2RXdeadTime = 10e-6;  % Receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
% HW.TX(iDevice).BlankOffset = 10e-3;  % offset of un-blanking signal before rf pulse in s
% HW.TX(iDevice).BlankPostset = 1e-3;  % additional un-blanking time after rf pulse


%% Coil without resistor (big)
% HW.RX2TXdeadTime = 5e-6;    % dead-time of receiver before transmission
% HW.TX2RXdeadTime = 2e-6;    % dead-time of reciver after transmission
HW.TX(iDevice).BlankOffset = 4000e-9;  % offset of un-blanking signal before rf pulse
HW.TX(iDevice).BlankPostset = 160e-9;  % additional un-blanking time after rf pulse


%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 5000e-9;   % blank of receiver before TX pulse in s
HW.TX(iDevice).BlankPostsetAQ = 1000e-9;  % blank of receiver after TX pulse in s
HW.TX(iDevice).BlankAQ = 1;               % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.


%% power consumption and duty cycle properties of (external) rf amplifier
HW.TX(iDevice).AmplPowerOutCW = 100;  % maximum CW power in W
HW.TX(iDevice).AmplBasicPower = 250;  % power consumption of rf amplifier when un-blanked in W
HW.TX(iDevice).AmplMaxPower = 300;  % maximum power consumption of rf amplifier before over-heating in W
HW.TX(iDevice).AmplPowerTime = 0.2;  % sliding window duration for analysis in s


%% RX Settings
% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1;  % reduce VGA gain to avoid saturation

% HW.RX.LNAGain=10^(((-61.12)-(-60))/20); % 25.8757 dB gain @ 24500397.7985 MHz F=0.82651 dB (-60 dBm cal)


%%
% UseRF100Switch = 1; % use switch of RF-100
% LoadRF100_Cal; % new cal Uout and 6 A FET
