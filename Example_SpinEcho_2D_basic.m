%% Demo Sequence "Example_SpinEcho_2D_basic.m"
% This demo sequence demonstrates how to create a basic Spin Echo 2D sequence.
% See "sequence_Spin_Echo" for a more versatile Spin Echo imaging function.

%% Basic Spin Echo 2D

%% Preparations
LoadSystem;                                          % load system parameters
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % determine magnet frequency
% HW.fLarmor = 21.947e6;                               % uncomment for entering frequency


%% Define parameters for a basic Spin Echo 2D sequence
nPhase   = 64;                  % phase encoding steps
nSamples = 64;                  % number of samples per acquisition
echoTime = 5e-3;                % Echo time in seconds
fSample  = 100e3;               % sample frequency in Hz
tRep     = 0.2;                 % repetition time in seconds
averages = 1;                   % number of averages (1=none)

grad_read  = 3;                 % gradient channel for read gradient
grad_phase = 1;                 % gradient channel for phase gradient

grad_read_amp  = 150e-3;        % gradient amplitude read gradient in Tesla / meter
grad_phase_amp = 150e-3;        % gradient amplitude phase gradient in Tesla / meter
grad_ramptime  = 100e-6;        % ramp time in seconds for all gradients

p90  = HW.tFlip90Def;           % rf-pulse duration 90 degrees pulse in seconds
p180 = HW.tFlip180Def;          % rf-pulse duration 180 degrees pulse in seconds


%% Pulse program
% Sequence parameters
Seq.plotSeq = [1, 3];               % plot the sequence: RF, AQ, Grad(1), Grad(3)
Seq.tRep    = ones(1,nPhase)*tRep;  % assign repetition time to Seq structure
Seq.average = averages;             % assign averages to Seq structure

% RF transmitter parameters
TX.Start        = [ 0; ...          % start time of 90 degrees rf-pulse
                    echoTime/2];    % start time of 180 degrees rf-pulse
TX.Duration     = [ p90; ...        % duration of 90 degrees rf-pulse
                    p180];          % duration of 180 degrees rf-pulse
TX.Frequency    = [ HW.fLarmor; ... % frequency of 90 degrees rf-pulse
                    HW.fLarmor];    % frequency of 180 degrees rf-pulse
TX.Phase        = [ 0; ...          % phase of 90 degrees rf-pulse
                    0];             % phase of 180 degrees rf-pulse

% Acquisition parameters
AQ.fSample      = [ fSample ];                              % sample frequency of acquisition
AQ.Start        = [ echoTime - ((nSamples / 2)/fSample)];   % start time of acquisition
AQ.nSamples     = [ nSamples ];                             % samples to acquire
AQ.Frequency    = [ HW.fLarmor];                            % frequency of acquisition
AQ.Phase        = [ 0 ];                                    % phase of acquisition

% Help variables
zero_matrix = zeros(1, nPhase);
ones_matrix = ones(1, nPhase);


% Gradients time and amplitude calculations
% Read gradient
Grad(grad_read).Time = [(echoTime - ((nSamples / 2)/fSample) - 2 * grad_ramptime - 1 * grad_ramptime - ((nSamples / 2)/fSample) - 1 * grad_ramptime); ...
                          (echoTime - ((nSamples / 2)/fSample) - 2 * grad_ramptime - 1 * grad_ramptime - ((nSamples / 2)/fSample)); ...
                          (echoTime - ((nSamples / 2)/fSample) - 2 * grad_ramptime - 1 * grad_ramptime); ...
                        (echoTime - ((nSamples / 2)/fSample) - 2 * grad_ramptime); ...
                          (echoTime - ((nSamples / 2)/fSample)); ...
                          (echoTime + ((nSamples / 2)/fSample)); ...
                        (echoTime + ((nSamples / 2)/fSample) + 2 * grad_ramptime)];

Grad(grad_read).Amp  = [0; ...
                          -grad_read_amp; ...
                          -grad_read_amp; ...
                        0; ...
                          grad_read_amp; ...
                          grad_read_amp; ...
                        0];


% Phase gradient
Grad(grad_phase).Time = [Grad(grad_read).Time(1)*ones_matrix;...
                           Grad(grad_read).Time(2)*ones_matrix;...
                           Grad(grad_read).Time(3)*ones_matrix;...
                         Grad(grad_read).Time(4)*ones_matrix];

% calculation of the amplitudes for the phase encoding steps
phase_amps = cumsum(ones_matrix * -grad_phase_amp/nPhase)+grad_phase_amp/2;

Grad(grad_phase).Amp  = [0*zero_matrix; ...
                           phase_amps; ...
                           phase_amps; ...
                         0*zero_matrix];

% Check if placement of read-gradient is after 2nd rf-pulse
if Grad(grad_read).Time(1) < TX.Start(2)
  error('Echotime too short.');
end


%% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);


%% Plot results
plot_data_1D(HW, data_1D);

% Display the result
kspace = squeeze(data(1).data);
imagespace = fftshift(ifft2(ifftshift(kspace)));

hf = figure(4);
hax = subplot(2,1,1, 'Parent', hf);
imagesc(log(abs(kspace)), 'Parent', hax);
title(hax, 'k-space');
xlabel(hax, 'phase direction');
ylabel(hax, 'read direction');
colormap(hax, gray);
hax = subplot(2,1,2, 'Parent', hf);
imagesc(abs(imagespace), 'Parent', hax);
title(hax, 'Image (2D FFT of the k-space)');
xlabel(hax, 'phase direction');
ylabel(hax, 'read direction');


% cleanup
clear ones_matrix;
clear zero_matrix;


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
