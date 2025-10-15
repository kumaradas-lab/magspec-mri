%% Settings for rf amplifier RF-100 with SN 6

HW.TX.ExtRFSN = 6;
HW.TX.ExtRFType = 'RF100';  % string with type identifier

HW.TX.ChannelDef = 2;         % Default TX Channel set to Tx2

% Amplifier
HW.TX.Uout2PaUout(2) = 50;  % 36 dB to 40 dB amplification
HW.TX.Max.PaUout(2) = 250;  % maximum peak output voltage in V

HW.TX.Def.PaUout(2) = 65;  % default peak output voltage in V (might be reduced by other HW.TX.Def settings)

HW.TX.Max.Amplitude = [20, 20]*1e-3;  % maximum B1+ in T
HW.TX.Def.Amplitude = [20, 20]*1e-3;  % default B1+ in T

%% Coil 10mm
HW.RX2TXdeadTime = 1e-6;       % Receiver deadtime before TX pulse in s
HW.TX2RXdeadTime = 10e-6;      % Receiver deadtime after TX pulse in s e.g. ringing of Coil ~40 µs
HW.TX.BlankOffset = 400e-9;    % Unblank of RF-100 before TX pulse in s
HW.TX.BlankPostset = 400e-9;    % Blank of RF-100 after TX pulse in s
%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ = 400e-9;  % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ = 2400e-9; % Blank of receiver after TX pulse
HW.TX.BlankAQ = 1;             % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%% RX Settings
% HW.RX.VGAGainDef = HW.RX.VGAGainMax/1; % default VGA gain
HW.RX.VGAGainDef = HW.RX.VGAGainMax/1; % reduce VGA gain to avoid saturation

%% 50 Ohm
% HW.RX2TXdeadTime = 1e-6;        % Receiver deadtime before TX pulse in s
% HW.TX2RXdeadTime = 400e-9;      % Receiver deadtime after TX pulse in s
% HW.TX.BlankOffset = 200e-9;     % Unblank of RF-100 before TX pulse in s
% HW.TX.BlankPostset = 80e-9;     % Blank of RF-100 after TX pulse in s
% HW.TX.BlankOffsetAQ = 300e-9;   % Blank of receiver before TX pulse
% HW.TX.BlankPostsetAQ = 1000e-9; % Blank of receiver after TX pulse
% HW.TX.BlankAQ = 1;              % Switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

%%
UseExtRFAmpSwitch = 1;  % 1: use switch of RF-100, 0: use active switch of LNA
LoadExtRFAmp_Cal;

%%
% if exist('CalibrationPaUout2Amplitude.mat','file');
%     load('CalibrationPaUout2Amplitude.mat');
%     HW.TX.CalibrationPaUout2Amplitude=CalibrationTX;
%     clear CalibrationPaUout2Amplitude
% end;
% if exist('CalibrationUout.mat','file');
%     load('CalibrationUout.mat');
%     HW.TX.CalibrationUout=CalibrationUout;
%     clear CalibrationUout
% end;
% if exist('CalibrationNorm2MmrtUout.mat','file');
%     load('CalibrationNorm2MmrtUout.mat');
%     HW.TX.CalibrationNorm2MmrtUout=CalibrationTX;
%     clear CalibrationNorm2MmrtUout
% end;
%
% if any(HW.TX.Uout2PaUout~=1)
%     if exist('CalibrationRfAmp.mat','file');
%         load('CalibrationRfAmp.mat');
%         HW.TX.CalibrationRfAmp=CalibrationRfAmp;
%         clear CalibrationRfAmp
%     end;
% end
%
% if exist('CalibrationRx.mat','file');
%     load('CalibrationRx.mat');
%     HW.RX.Calibration=CalibrationRx;
%     clear CalibrationRx
% end;
