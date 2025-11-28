%% T2 spectrum using a solid echo followed by a CPMG echo train
%
% NOTE: Solid echo pulses and their ring-down time should be short.
%       Results might vary for systems that do not fulfill that requirement.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% preparation
LoadSystem;  % load system parameters (reset to default: HW Seq AQ TX Grad)

% reduce dead time to allow early acquisition of solid echo
HW.TX2RXdeadTime = 3e-6;  % receiver deadtime after TX pulse in s


%% adjustable settings
Seq.sampleName = 'Test_Measurement';

% measurement settings
Seq.tEcho = max(250e-6, HW.RecoveryCPMG.tEchoMin);  % echo time in s
% acquisition window duration at echoes in CPMG echo train in seconds
Seq.tAQEcho = 100e-6;
Seq.nEchoes = 8000;  % number of echoes in CPMG echo train

Seq.tauSolidEcho = max(10e-6, 2*HW.tFlip90Def);  % solid echo time in seconds
% acquisition duration of solid echo in seconds
Seq.tAQFID = min(100e-6, ...
  Seq.tEcho/2 ...
  - (HW.tFlip90Def/2 + get_DeadTimeTX2RX(HW, HW.RecoveryCPMG.fSampleFID) ...
     + HW.tFlip180Def/2 + get_DeadTimeRX2TX(HW, HW.RecoveryCPMG.fSampleFID)));

Seq.tRelax = 15;  % relaxation time after CPMG Echo train in s

Seq.nAverages = 1;  % number of averages

Seq.phaseCycling = true;  % phase cycling on or off

Seq.saveResults = false;  % save measurement data and results for each measurement step
% path that is used for storing results
pathName = fullfile(HW.RootPath, 'output', ...
  'Solid Echo', [datestr(now, 'yyyymmdd-HHMMSS') '_results']);

nLoop = 1;  % number of repetitions

tRepetitionLoop = 60;  % minimum repetition time of loop in seconds


%% measurement settings

Seq.Plot = 0;  % plot measured signal
Seq.plotSeq = [];

% duration of the refocusing pulses in seconds
Seq.tRefocus = HW.RecoveryCPMG.tFlip180Def;
if isempty(Seq.tRefocus)
  Seq.tRefocus = HW.tFlip180Def;
end
Seq.tEchoTrain = Seq.nEchoes * Seq.tEcho;  % total time of echo train in s

Seq.FixfLarmorTotEcho = false;  % adjust the Larmor frequency to the Echo time

% rf pulse settings
Seq.excitationPulse = @Pulse_Rect_SolidEcho;
Seq.tExcitation = Seq.tRefocus/2;
Seq.refocusingPulse = @Pulse_Rect;

Seq.tAQEchoDelay = HW.RecoveryCPMG.tAQEchoDelay;  % delays the acquisition windows in s

% optional T1 preparation steps
Seq.nTau1 = 0;  % number of recovery times. Set to 0 for CPMG only
Seq.nTau1SteadyState = 1; % number of pre-shots with 1st Tau1
% Seq.Tau1Start = 0.20e-3;  % shortest recovery time for automatic spacing in s
% Seq.Tau1End = 1000e-3;  % longest recovery time for automatic spacing in s
% Seq.Tau1Log = 1;  % Boolean to select logarithmic spacing of recovery times (linear spacing otherwise)
% Seq.Recovery = 'Inversion';  % type of recovery experiment, 'Inversion' or 'Saturation'

if Seq.phaseCycling
  Seq.SeqAverage.average = 2*Seq.nAverages;  % set to 2 for phase cycling
  % Cycle phase of excitation and acquisition instead of refocusing pulses
  Seq.SeqAverage.TXPhaseExcitationIncrement = 180;  % phase increment between averages of excitation pulse (90er of CPMG)
  Seq.SeqAverage.TXPhaseRefocusIncrement = 0;  % phase increment between averages of 180er of CPMG
  Seq.SeqAverage.AQPhaseIncrement = 180;  % phase increment between averages of acquisition window at Echoes
else
  Seq.SeqAverage.average = Seq.nAverages;
end


%% evaluation settings

Seq.SubtractEmptyMeasurement = false;  % Run the same measurement without sample to subtract it.
Seq.doT2iLaplace1D = true;  % inverse Laplace transform
Seq.FitT2 = true;  % single and double exponential fit
Seq.fitExpT2.CorrectFrequencyDrift = false;
Seq.fitExpT2.CorrectFrequencyOffset = false;
Seq.fitExpT2.CorrectPhaseOffset = false;
Seq.fitExpT2.EndOffset = false;
Seq.fitExpT2.CorrectPhaseOffsetAtEnd = false;
Seq.fitExpT2.SingleExp = true;
Seq.fitExpT2.DoubleExp = true;


%% display settings
Seq.plotSeq = 0;  % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)

Seq.Plot = 1;
if Seq.doT2iLaplace1D
  Seq.iLaplace1D.hParent = 201;
end

Seq.PlotT2 = false;  % plots are done manually in the measurement loop
Seq.doDisplay = false;


%% measurement loop
time_prev = 0;

for iLoop = 1:nLoop
  %% timing of loop repetitions
  % The timing of the loop repetitions is only approximate with some expected
  % delays and jitter.
  if iLoop > 1
    waitTime = (time_prev-now())*24*3600 + tRepetitionLoop;
    fprintf('Waiting %.1f seconds for repetition #%d of %d...', ...
      waitTime, iLoop, nLoop);
    pause(waitTime);
    fprintf('done\n');
  end
  time_prev = now();


  %% Frequency sweep
  HW.FindFrequencySweep.maxTime = 10;
  % HW.FindFrequencySweep.nMeasurements = 1;
  HW.FindFrequencyPause = 10;
  Seq.Find_Frequency_interval = true;  % let the sequence decide when to find the magnet frequency


  %% actual measurement
  [data, SeqOut, mySave] = sequence_RecoveryCPMG(HW, Seq, mySave);


  %% evaluation
  [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data);


  %% display results
  hf_results = clf(figure(1001));
  hax = axes('Parent', hf_results);
  T2 = data.T2;
  hl = plot(hax, T2(SeqOut.FitT2AtTau1(1)).timeCorrected, ...
    T2(SeqOut.FitT2AtTau1(1)).dataPhaseCorrectedReal, 'x');
  hold(hax, 'on');
  legendElem = hl;
  legendEntry = {'real'};
  if SeqOut.fitExpT2.DoubleExp
   legendElem(end+1) = plot(hax, T2(SeqOut.FitT2AtTau1(1)).functionTimeLog, ...
      T2(SeqOut.FitT2AtTau1(1)).functionAmpDoubleLog, '-.');
    legendEntry{end+1} = 'fit double';
  end
  if SeqOut.fitExpT2.SingleExp
    legendElem(end+1) = plot(hax, T2(SeqOut.FitT2AtTau1(1)).functionTimeLog, ...
      T2(SeqOut.FitT2AtTau1(1)).functionAmpSingleLog, '--');
    legendEntry{end+1} = 'fit single';
  end
  set(hax, 'XScale', 'log');
  set(hax, 'YLim', [0, 1.1*max(T2(SeqOut.FitT2AtTau1(1)).dataPhaseCorrectedReal(:))]);
  grid(hax, 'on');

  xlabel(hax, 'time in s');
  ylabel(hax, 'amplitude');
  title(hax, {[' ', ' ', ' ']});
  grid(hax, 'on');
  legend(legendElem, legendEntry);


  %% export results to disk
  if SeqOut.saveResults
    if ~exist(pathName, 'dir')
      mkdir(pathName);
    end

    fileName = sprintf('fit_%s_%03d.png', SeqOut.sampleName, iLoop);
    if verLessThan('MATLAB', '25.1')
      print(hf_results, '-dpng', fullfile(pathName, fileName));
    else
      if verLessThan('MATLAB', '26.1')
        % FIXME: It is unclear when this issue will be fixed in MATLAB. Adjust
        %        the above version check accordingly.
        % Docked figures that do not have the focus might be printed black in
        % MATLAB R2025a and MATLAB R2025b. To work around that, try to focus
        % the figure before exporting the .png image.
        focus(hf(iFig));
      end
      exportapp(hf_results, fullfile(pathName, fileName));
    end

    if SeqOut.doT2iLaplace1D
      fileName = sprintf('iLaplace_%s_%03d.png', SeqOut.sampleName, iLoop);
      if verLessThan('MATLAB', '25.1')
        print(SeqOut.iLaplace1D.hParent, '-dpng', fullfile(pathName, fileName));
      else
        if verLessThan('MATLAB', '26.1')
          % FIXME: It is unclear when this issue will be fixed in MATLAB. Adjust
          %        the above version check accordingly.
          % Docked figures that do not have the focus might be printed black in
          % MATLAB R2025a and MATLAB R2025b. To work around that, try to focus
          % the figure before exporting the .png image.
          focus(hf(iFig));
        end
        exportapp(SeqOut.iLaplace1D.hParent, fullfile(pathName, fileName));
      end
    end

    fileName = sprintf('data_%s_%03d.mat', SeqOut.sampleName, iLoop);
    data = rmfield(data, 'SeqAverage');
    data = rmfield(data, 'lift');
    data.T2 = rmfield(data.T2, 'hParent');
    data.T2 = rmfield(data.T2, 'hFigure');
    data.T2.fitExpSettings = rmfield(data.T2.fitExpSettings, 'hParent');
    save(fullfile(pathName, fileName), 'data');
  end

end
