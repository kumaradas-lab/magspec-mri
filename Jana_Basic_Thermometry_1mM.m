%% This code calculates the Temperature of a Mn-Sample at known concentration by measuring T1
% it uses only one Reference curve (the one measured by Khalid at a
% concentration of 1mM (botom/pink line of figure on slide 7)

% Preparations
LoadSystem;                                     % Load system parameters

%DEFINE EXPERIMENT PARAMETERS
%define concentration of the sample used
concentration = 1; %input concentration in mM

% Define parameters required for the measurement, e.g. oil T1 ~100e-3 sec
Seq.T1Estimated   = 100e-3;                     % estimated mean T1
% Seq.T1EstimatedMin = Seq.T1Estimated/3;         % minimum estimated T1,                         e.g. Seq.T1Estimated/3
% Seq.T1EstimatedMax = Seq.T1Estimated*3;         % maximum estimated T1,                         e.g. Seq.T1Estimated*3

Seq.tRelax        = Seq.T1Estimated*5;          % relaxation time from excitation pulse to next inversion pulse (0 for single-shot inversion recovery FLASH), e.g. 500e-3 (T1*5)
% Seq.tFlip         = 10e-3;                      % flip repetition time,                         e.g. 10e-3          (T1/10)
% Seq.excitationFlipAngle = 3;                    % flip angle of excitation pulse in degrees,    e.g. 3              (30 degrees/(T1/Seq.tFlip))
% Seq.nFids         = 50;                         % number of measured FIDs,                      e.g. 50             (T1*20/Seq.tFlip)
Seq.recovery      = 'Inversion';                % type of recovery,                             e.g. 'Inversion'    ('Inversion' or 'Saturation')
Seq.ConsoleOut    = 1;                          % display results in console,                   e.g. 1              (1 or 0)
Seq.average       = 1;                          % number of averages,                           e.g. 1              (>=1)
Seq.averageBreak  = Seq.tRelax;                 % time between two averages,                    e.g. 0.3            (T1*3)

Seq.tAQFID        = 200e-6;                     % acquisition time in seconds,                  e.g. 200e-6
Seq.tFlipLog      = 1;                          % increase wait time logarithmically (1) or linearly (0)
Seq.Spoil.UseCoordinate = 1;                    % Coordinate direction of the spoiler,          e.g. 1 or 2         (0 for no spoiler)

Seq.fitExp.CorrectFrequencyOffset = 1;          % Correct frequency offset before averaging samples
Seq.fitExp.CorrectFrequencyDrift = 0;           % Also correct linear frequency drift

% Start measurement by calling the pre-made sequence
HW.FindFrequencyPause = Seq.T1Estimated*10;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % Find magnet frequency.

[t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq);  % Actual measurement. See help for more options.

% use determined T1 to obtain temperature of the sample 

% calculate relaxivity r1
t1 = t1*0.001; % convert t1 from ms to s
r1 = 1/t1/concentration; % make relaxivity concentration independant

disp(['Relaxivity r1 = ' num2str(r1) ' Hz']);


syms x;
eqn = r1 == -0.0013048*x^2 + 0.027533*x + 5.2391;

S = solve(eqn, x);
Temperature = double(abs(S)); % to show S numerically
disp(['The Temperature is T = ' Temperature(2) ' deg. C']);