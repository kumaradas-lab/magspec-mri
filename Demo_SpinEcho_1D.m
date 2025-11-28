%% Demo Sequence "Simple Spin Echo 1D"
% This demo sequence demonstrates how to create a simple spin echo with
% gradients. The FFT is applied to the encoded spin echo signal which
% shows a profile of the sample.

%% Simple Spin Echo 1D
% Preparations
    LoadSystem                                      % load system parameters
    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1); % find magnet frequency

    useGrad = 3;  % choose which gradient to be used (1==x, 2==y, 3==z)

% Define parameters required for the measurement

    % Parameters used for timing calculations
    Seq.tEcho       =   10e-3;                      % echo time
    Seq.p90         =   HW.tFlip90Def;              % duration of 1st TX pulse; use default value
    Seq.p180        =   HW.tFlip180Def;             % duration of 2nd TX pulse; use default value
    Seq.plotSeq     =   useGrad;                    % plot sequence on, plot RF, AQ and Grad(2)

    % Sequence parameters
    Seq.tRep        =   100e-3;                     % repetition time
    Seq.average     =   1;                          % averages

    % RF transmission parameters
    TX.Start        = [ -Seq.p90/2; ...             % start time of 1st TX pulse
                        Seq.tEcho/2-Seq.p180/2];    % start time of 2nd TX pulse
    TX.Duration     = [ Seq.p90; ...                % duration of 1st TX pulse
                        Seq.p180];                  % duration of 2nd TX pulse
    TX.Frequency    = [ HW.fLarmor; ...             % frequency of 1st TX pulse
                        HW.fLarmor];                % frequency of 2nd TX pulse
    TX.Phase        = [ 0; ...                      % phase of 1st TX pulse
                        90];                        % phase of 2nd TX pulse

    % Acquisition parameters
    AQ.Start        = [ 100e-6; ...                 % acquisition start of 1st AQ window
                        8500e-6];                   % acquisition start of 2nd AQ window
    AQ.fSample      = [ 30e3; ...                   % sampling rate of 1st AQ window
                        100e3];                     % sampling rate of 2nd AQ window
    AQ.nSamples     = [ 128; ...                    % number of samples in 1st AQ window
                        300];                       % number of samples in 2nd AQ window
    AQ.Frequency    = [ HW.fLarmor; ...             % frequency of 1st AQ window
                        HW.fLarmor];                % frequency of 2nd AQ window
    AQ.Phase        = [ 0; ...                      % phase of 1st AQ window
                        145];                       % phase of 2nd AQ window

    % Gradient along x, y or z (1==x, 2==y, 3==z)
    Grad(useGrad).Time = [ 0; 1; 2; 3; 4; 7.5; 8.5; 11.5; 12.5]*1e-3;   % time of gradient points in seconds
    Grad(useGrad).Amp  = [ 0; 0; 1; 1; 0; 0;   1;   1;    0;  ]*100e-3; % amplitude of gradient points in Tesla/meter

% Start measurement
    [Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

% Plot results
    plot_data_1D(HW, data_1D);

% Plot FFT
    figure(31);
    plot((data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/HW.Gamma.H1*2*pi/max(Grad(useGrad).Amp), ...
      abs(data(1).fft1_data(1:AQ.nSamples(2),2)));
    xlabel('distance/meter');
    ylabel('amplitude');

    figure(32);
    plot((data(1).f_fft1_data(1:AQ.nSamples(2),2))-HW.fLarmor, ...
      abs(data(1).fft1_data(1:AQ.nSamples(2),2)));
    xlabel('frequency');
    ylabel('amplitude');

    clear usegrad;

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
