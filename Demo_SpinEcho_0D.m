%% Demo Sequence "Simple Spin Echo 0D"
% This demo sequence demonstrates how to create a basic Spin Echo.

%% Simple Spin Echo 0D
% % Preparations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadSystem;                                   % load system parameters

[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency


% % Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters required for the measurement

% Parameters used for timing calculations
Seq.tEcho     = 50e-3;                        % echo time in s
Seq.p90       = HW.tFlip90Def; % or 45e-6;    % duration of 1st TX pulse in s
Seq.p180      = HW.tFlip180Def; % or 90e-6;   % duration of 2nd TX pulse in s
Seq.plotSeq   = [];                           % plot sequence off

% Sequence parameters
Seq.tRep      = 100e-3;                       % repetition time in s


% RF transmission parameters
TX.Start      = [ -Seq.p90/2; ...             % start time of 1st TX pulse
                  Seq.tEcho/2-Seq.p180/2];    % Start Time of 2nd TX pulse
TX.Duration   = [ Seq.p90; ...                % duration of 1st TX pulse
                  Seq.p180];                  % duration of 2nd TX pulse
TX.Frequency  = [ HW.fLarmor; ...             % frequency of 1st TX pulse
                  HW.fLarmor];                % frequency of 2nd TX pulse
TX.Phase      = [ 0; ...                      % phase of 1st TX pulse
                  -90];                       % phase of 2nd TX pulse

% Acquisition parameters
AQ.fSample    = [ 30e3; ...                   % sampling rate of 1st AQ window
                  50e3];                      % sampling rate of 2nd AQ window
deadTime = get_DeadTimeTX2RX(HW, AQ.fSample); % dead time after TX pulse
AQ.Start      = [ TX.Start(1)+TX.Duration(1)+deadTime(1); ...  % acquisition start 1st AQ window
                  TX.Start(2)+TX.Duration(2)+deadTime(2)];     % acquisition start 2nd AQ window
AQ.nSamples   = [ 128; ...                    % number of samples in 1st AQ window
                  2048];                      % number of samples in 2nd AQ window
AQ.Frequency  = [ HW.fLarmor; ...             % frequency of 1st AQ window
                  HW.fLarmor];                % frequency of 2nd AQ window
AQ.Phase      = [ 0; ...                      % phase of 1st AQ window
                  0];                         % phase of 2nd AQ window


% % Start measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);


% % Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_data_1D(HW, data_1D);


%% -----------------------------------------------------------------------------
% (C) Copyright 2012-2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
