%% Demo Sequence "Spin Echo 0D Parameter Search"
% This demo sequence demonstrates how to find parameters like Larmor-frequency
% fLarmor, pulse duration Seq.p90 and shim values Grad.Shim when e.g. customized
% hardware is used.

%% Simple Spin Echo 0D for parameter search
% % Preparations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadSystem;                                   % load system parameters


% % Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to be changed by the user (use the "increment value near cursor and
% evaluate cell" function provided by the Matlab editor)
Seq.p90         = 40e-6;                      % duration of 1st TX pulse in s
% HW.fLarmor      = 24631000;                   % Larmor/system frequency
Grad(1).Shim    = 10e-3;                      % x-shim gradient in T/m
Grad(2).Shim    = 0e-3;                       % y-shim gradient in T/m
Grad(3).Shim    = 0e-3;                       % z-shim gradient in T/m
Grad(4).Shim    = 0e-3;                       % B0 shift (optional)

% HW.B0 = 0.508;                              % search frequency around HW.B0 in Hz
% [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 100);  % find magnet frequency


% Define parameters required for the measurement

% Parameters used for timing calculations
Seq.tEcho       = 10e-3;                      % echo time in s
Seq.p180        = 2*Seq.p90;                  % duration of 2nd TX pulse in s

% Sequence parameters
Seq.tRep        = 100e-3;                     % repetition time in s


% RF transmission parameters
TX.Start        = [ -Seq.p90/2; ...           % start time of 1st TX pulse
                    Seq.tEcho/2-Seq.p180/2];  % start time of 2nd TX pulse
TX.Duration     = [ Seq.p90; ...              % duration of 1st TX pulse
                    Seq.p180];                % duration of 2nd TX pulse
TX.Frequency    = [ HW.fLarmor; ...           % frequency of 1st TX pulse
                    HW.fLarmor];              % frequency of 2ns TX pulse
TX.Phase        = [ 0; ...                    % phase of 1st TX pulse
                    90];                      % phase of 2nd TX pulse

% Acquisition parameters
AQ.Start        = [ 100e-6; ...               % acquisition start of 1st AQ window
                    5100e-6];                 % acquisition start of 2nd AQ window
AQ.fSample      = [ 30e3; ...                 % sampling rate of 1st AQ window
                    100e3];                   % sampling rate of 2nd AQ window
AQ.nSamples     = [ 128; ...                  % number of samples in 1st AQ window
                    1024];                    % number of samples in 2nd AQ window
AQ.Frequency    = [ HW.fLarmor; ...           % frequency of 1st AQ window
                    HW.fLarmor];              % frequency of 2nd AQ window
AQ.Phase        = [ 0; ...                    % phase of 1st AQ window
                    145];                     % phase of 2nd AQ window


% % Start measurement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);


% % Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_data_1D(HW, data_1D);


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
