%% rf amplifier @ > ~60MHz

if ~exist('iDevice', 'var'),  iDevice = 1;  end

HW.TX(iDevice).ChannelDef = 2;  % Default TX Channel

% Amplifier
% HW.TX(iDevice).Uout2PaUout = [1, 10^((52-20)/20)];    % 20dB attenuation at input 52 +/- 1 dB Amplification
HW.TX(iDevice).Uout2PaUout = [1, 10^((52-18-5)/20)];  % 20dB attenuation at input 52 +/- 1 dB Amplification
% HW.TX(iDevice).Uout2PaUout = [1, 10^((52-16)/20)];    % 10dB attenuation at input 52 +/- 1 dB Amplification
% HW.TX(iDevice).Uout2PaUout = [1, 10^((52-0)/20)];     % no attenuation at input 52 +/- 1 dB Amplification

HW.TX(iDevice).Max.PaUout(2) = 110;             % 121 W
HW.TX(iDevice).Def.PaUout(2) = 5;               % Set Def output amplitude in Volt peak-peak

HW.TX(iDevice).Max.Amplitude = [20, 20]*1e-3;   % maximum output in B1 amplitude
HW.TX(iDevice).Def.Amplitude = [20, 20]*1e-3;   % default output in B1 amplitude

%% Coil without resistor (big)
% HW.RX2TXdeadTime = 5e-6;    % dead-time of receiver before transmission
% HW.TX2RXdeadTime = 2e-6;    % dead-time of reciver after transmission
HW.TX(iDevice).BlankOffset = 04000e-9;  % offset of un-blanking signal before rf pulse
HW.TX(iDevice).BlankPostset = 160e-9;   % additional un-blanking time after rf pulse
%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 5000e-9;      % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 1000e-9;     % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;                  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.


%% RX Settings
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1;  % default receiver amplification
% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % receiver amplification for shorter dead-time

%% 50 Ohm
% HW.RX2TXdeadTime = 1e-6;    % dead-time of receiver before transmission
% HW.TX2RXdeadTime = 400e-9;  % dead-time of reciver after transmission
% HW.TX(iDevice).BlankOffset = 200e-9;  % offset of un-blanking signal before rf pulse
% HW.TX(iDevice).BlankPostset = 80e-9;  % additional un-blanking time after rf pulse
% HW.TX(iDevice).BlankOffsetAQ = 300e-9;    % blank of receiver before TX pulse
% HW.TX(iDevice).BlankPostsetAQ = 1000e-9;  % blank of receiver after TX pulse
% HW.TX(iDevice).BlankAQ = 1;               % switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.


%%
% if exist('CalibrationPaUout2Amplitude.mat', 'file');
%   load('CalibrationPaUout2Amplitude.mat');
%   HW.TX(iDevice).CalibrationPaUout2Amplitude = CalibrationTX;
%   clear CalibrationPaUout2Amplitude
% end
% if exist('CalibrationUout.mat', 'file');
%   load('CalibrationUout.mat');
%   HW.TX(iDevice).CalibrationUout = CalibrationUout;
%   clear CalibrationUout
% end
% if exist('CalibrationNorm2MmrtUout.mat', 'file');
%   load('CalibrationNorm2MmrtUout.mat');
%   HW.TX(iDevice).CalibrationNorm2MmrtUout = CalibrationTX;
%   clear CalibrationNorm2MmrtUout
% end;
%
% if any(HW.TX(iDevice).Uout2PaUout~=1)
%   if exist('CalibrationRfAmp.mat', 'file');
%     load('CalibrationRfAmp.mat');
%     HW.TX(iDevice).CalibrationRfAmp = CalibrationRfAmp;
%     clear CalibrationRfAmp
%   end
% end
%
% if exist('CalibrationRx.mat', 'file');
%   load('CalibrationRx.mat');
%   HW.RX(iDevice).Calibration = CalibrationRx;
%   clear CalibrationRx
% end
