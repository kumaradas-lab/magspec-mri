%% Reset the values of the structures Seq, TX, AQ, and Grad
% Previous values are cleared and the structures are re-created with the most
% commonly used fields.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Initialize blank structures
Seq = struct();
TX = struct();
AQ = struct();
Grad = struct();


%% Sequence parameters

Seq.tRep = [];            % array containing the lengths of the repetition times in seconds
% %         tRep1   tRep2   tRep3
% Seq.tRep=[ 0.5    0.2     2    ];
% TX.Start=[ 0.1    0.01    0.01 ;...   % Pulse 1
%            0.2    0.04    0.3  ;...   % Pulse 2
%            0.3    0.1     1    ];     % Pulse 3
% or
% %         tRep1   tRep2   tRep3
% Seq.tRep=[ 0.5    0.2     2    ];
% TX.Start=[ 0.1    nan    0.01 ;...   % Pulse 1
%            0.2    nan    0.3  ;...   % Pulse 2
%            0.3    nan    nan    ];   % Pulse 3

Seq.CLTime = [];          % array containing the time required to reload parameters
Seq.StartSequenceTime = [];  % time and date in seconds when sequence will be started
Seq.average = [];         % number of averages; set to 1 for no averaging
Seq.averageBreak = [];    % waiting time in seconds between two measurements for averaging

Seq.plotSeqTR = [];       % plot the running sequence; multiple TRs are plotted in different colors into the first TR. [] no Plot, 0 only TX and AQ, [1] X-Gradient TX and AQ, [1:3] XYZ Gradient TX and AQ
Seq.plotSeq = [];         % plot the running sequence; [] no Plot, 0 only TX and AQ, [1] X-Gradient TX and AQ, [1:3] XYZ Gradient TX and AQ


% Transmission parameters

TX.Channel = HW.TX.ChannelDef;  % use transmitting channel 1 or 2
TX.Start = [];            % matrix containing the starting time of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR (column)
TX.Frequency = [];        % matrix containing transmitting frequency for each pulse of each TR
TX.Duration = [];         % matrix containing the duration of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR
TX.Amplitude = [];        % matrix containing the amplitude for each pulse of each TR
TX.Phase = [];            % matrix containing the phase offset for each pulse of each TR

TX.BlankOffset = [];      % array of times to define, when the blanking signal is applied prior to the RF pulse; can be adjusted every TR
TX.BlankPostset = [];     % array of times to define how long blanking stays active after a RF pulse; can be adjusted every TR
TX.Repeat = [];           % if there is a ‘1’ in the array, the pulses of the prior TR are used. (reduces data traffic between the PC and the MRI device)
TX.Device = 1;            % multiple devices (1: master)


%% Acquisition parameter
% AQ = mriDevice.sequence.AQ;
AQ.Start = [];            % matrix containing the starting times (center of the first Sample – 0.5/AQ.fSample) of each acquisition window for each TR; up to 510 AQ windows are possible per TR
AQ.nSamples = [];         % matrix containing the number of samples to be acquired for each AQ window for each TR
AQ.fSample = [];          % matrix containing the sampling frequency of each acquisition window for each TR
AQ.Frequency = [];        % matrix containing the mixing frequency of each acquisition windows for each TR
AQ.Phase = [];            % matrix containing the phase offset for each acquisition window
AQ.SamplingFactor = [];   % matrix containing the oversampling factor for each acquisition window
AQ.ResetPhases = [];      % if set to 1 the TX and RX phase is matched at the beginning of the TR
AQ.Gain = [];             % relative gain of the acquisition path (1/AQ.Gain = max input Amplitude). Use HW.RX.GainDef for the best noise figure.
AQ.Repeat = [];           % if there is a ‘1’ in the array, the settings for the window of the prior TR are used. (reduces data traffic between the PC and the MRI device)
AQ.Device = 1;            % multiple devices (1: master)


%% Gradient parameters
% Grad(1) x
% Grad(2) y
% Grad(3) z
% Grad(4) B0

[Grad(1:sum([HW.Grad.n])).Time] = deal([]);   % time the values in Grad(t).Amp with the corresponding index are set
[Grad(1:sum([HW.Grad.n])).Amp] = deal([]);    % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
[Grad(1:sum([HW.Grad.n])).Shim] = deal([]);   % additional shim; magnet shim is already considered in HW.MagnetShim. Caution: Using high values over a long time will damage the gradient coils and amplifiers!
[Grad(1:sum([HW.Grad.n])).Repeat] = deal([]); % if there is a one in the array the gradients of the prior TR are used. (reduces data traffic between the PC and the MRI device)
