%% Demo Sequence "Demo_T2"
% This demo sequence demonstrates how to use a pre-made sequence. In this
% case the sequence get_T2.m is used for T2 measurement.

%% Measure T2 values using get_T2
% Preparations
    LoadSystem                                      % load system parameters

% Define parameters required for the measurement

    Seq.T1Estimated     = 0.2;                      % estimated T1 in seconds, e.g. oil T1 ~200e-3 sec
    Seq.T2Estimated     = 0.1;                      % estimated T2 in seconds, e.g. oil T2 ~100e-3 sec
    Seq.T2EstimatedMin  = [];                       % shortest estimated T2 in seconds
    Seq.T2EstimatedMax  = [];                       % longest estimated T2 in seconds

    Seq.tEcho           = [];                       % echo time, e.g. 4e-3
    Seq.tEchoTrain      = [];                       % duration of echo train, e.g. 1       (T2*5)
    Seq.tAQFID          = -1;                       % duration of FID to acquire in seconds
    Seq.plot            = 1;                        % plot data
    Seq.RawData         = 0;
    Seq.plotSeq         = 1;                        % plot sequence (no gradients: 0)
    Seq.ConsoleOut      = 1;                        % display results in console
    Seq.SeqAverage.average = 2;                     % number of averages, e.g. 2 (>=2 and even for suppression of stimulated echoes)
    Seq.SubtractEmptyMeasurement = 0;               % Run the same measurement without sample to subtract it.
    Seq.preparationPulse = @Pulse_Rect;             % pulse shape used for flip pulse
    Seq.inversionPulse = @Pulse_Rect;               % pulse shape used for inversion pulse


    HW.FindFrequencyPause = Seq.T1Estimated*3;      % choose a time longer than 3*T1 to avoid saturation
    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0); % find magnet frequency

    Seq.fitExp.CorrectFrequencyOffset = 1;          % Correct frequency offset before averaging samples
    Seq.fitExp.CorrectFrequencyDrift = 0;           % Also correct linear frequency drift
    Seq.fitExp.FIDMean = 0;                         % Use average of FID. If false, each sample of the FID is used for the fit.

% Start measurement by calling the pre-made sequence
    [t2, T2, data, SeqOut] = get_T2(HW, Seq);

 %% Get the inverse Laplace 1D
    SeqOut{1}.iLaplace1D.SpectrumTimeStart = data{1}.DataTime(1)/1;
    SeqOut{1}.iLaplace1D.SpectrumTimeEnd = data{1}.DataTime(end)*2;
    SeqOut{1}.iLaplace1D.QualityFactor = 50;        % sharpness of peaks e.g. 50
    [data{1}, SeqOutIL] = get_iLaplace1D(data{1}, SeqOut{1});

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
