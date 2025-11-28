%% Demo Sequence "Simple Spin Echo 0D with multiple TRs"
% This script demonstrates how to create a Spin Echo pulse sequence using
% multiple tReps (repetition times). The echo time and repetition time are set
% to different values in each tRep.

%% Simple Spin Echo 0D with multiple tReps (repetition times)
% % Preparations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadSystem;                                     % load system parameters

[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency

% % Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters required for the measurement

% Parameters used for timing calculations
Seq.tEcho       = [20e-3, 50e-3, 100e-3];       % echo time; 1st TR: 20 ms, 2nd TR: 50 ms, 3rd TR: 100 ms
Seq.p90         = HW.tFlip90Def; % or 45e-6;    % duration of 1st TX pulse
Seq.p180        = HW.tFlip180Def; % or 90e-6;   % duration of 2nd TX pulse
Seq.plotSeq     = [];                           % plot sequence off

% Sequence parameters
% Defining an array for Seq.tRep will perform multiple measurements.
% Here: size(Seq.tRep) = 3 which means 3 measurements will be executed
Seq.tRep        = [500e-3, 400e-3, 300e-3];


% RF transmission parameters for all TRs
TX.Start        = [-Seq.p90/2*ones(size(Seq.tEcho)); ... % start time of 1st TX pulse
                    Seq.tEcho/2-Seq.p180/2];    % start time of 2nd TX pulse
TX.Duration     = [ Seq.p90; ...                % duration of 1st TX pulse
                    Seq.p180];                  % duration of 2nd TX pulse
TX.Frequency    = [ HW.fLarmor; ...             % frequency of 1st TX pulse
                    HW.fLarmor];                % frequency of 2nd TX pulse
TX.Phase        = [ 0; ...                      % phase of 1st TX pulse
                    90];                        % phase of 2nd TX pulse

% Acquisition parameters for all TRs
AQ.fSample      = 25e3;                         % sampling rate for all AQ windows
deadTime = get_DeadTimeTX2RX(HW, AQ.fSample);
AQ.Start        = [ TX.Start(1,:)+TX.Duration(1,:)+deadTime; ... % acquisition start of 1st AQ window
                    TX.Start(2,:)+TX.Duration(2,:)+deadTime];    % acquisition start of 2nd AQ window
AQ.nSamples     = [ 128; ...                    % number of samples in 1st AQ window
                    2048];                      % number of samples in 2nd AQ window


% % Start measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ Raw, SeqOut, data, data_1D ] = set_sequence(HW, Seq, AQ, TX, Grad);


% % Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_data_1D(HW, data_1D);


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
