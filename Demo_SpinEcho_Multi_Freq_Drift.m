%% Demo Sequence "Echo Multi Freq Drift"
% This demo sequence demonstrates how to create an echo train using
% "sequence_EchoStandard".
% The measurement is repeated in a loop. The system frequency drift caused by
% the magnet temperature drift is tracked and compensated every loop iteration.


%% Echo Multi Freq Drift
% Preparations
LoadSystem;                                                   % load system parameters
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, [], 1);    % find magnet frequency


% Define parameters required for the measurement
% Sequence parameter
Seq.tEcho       = 5e-3;                                       % echo time
Seq.nEchos      = 5;                                          % number of echoes
Seq.fast        = 0;                                          % fast mode for tEcho < 2 ms
Seq.p90         = get_PulseDuration(HW, 90, HW.TX(1).AmpDef); % duration of 1st TX pulse
Seq.p180        = get_PulseDuration(HW, 180, HW.TX(1).AmpDef);% duration of 2nd TX pulse
Seq.plotSeq     = 0;                                          % plot sequence without gradients
Seq.plot        = 1;                                          % plot data
Seq.AQEcho      = 0.8;                                        % ratio of the echo time where data will be acquired [0...1[, use 0 for acquiring one sample
Seq.AQFID       = 0.8;                                        % ratio of the echo time / 2 where data will be acquired ]0...1[

Seq.Shim        = [[0, 0, 0]*1e-3, ...                        % shim values [x,y,z] in T/m
                   00000/(HW.Gamma.H1/2/pi)];                 % B0 offset in Tesla (B0 coil not available by default)
Seq.fSample     = HW.RX(1).fSample/125/10;                    % sample frequency of the echo HW.RX.fSample = 125e6 Hz possible dividers 1,[4,5,...,8191]
Seq.fSampleFID  = HW.RX(1).fSample/125/10;                    % sample frequency of the FID  HW.RX.fSample = 125e6 Hz possible dividers 1,[4,5,...,8191]

Seq.average     = 1;                                          % averages
Seq.averageBreak= 1;                                          % waiting time in seconds between two measurements for averaging
Seq.loops       = 50;                                         % loops; measurements are repeated, no averaging
Seq.LoopBreak   = 0.5;                                        % waiting time in seconds between two loops
Seq.fOffsetTime = 0.5e-3;                                     % maximum time used to determine the frequency offset for every loop


% prepare figure for frequency drift
hf = 10;
if ~isgraphics(hf, 'figure')
  hf = figure(hf);
  % move to right to not cover data figure at its default position
  fPos = get(hf, 'Position');
  fPos(1) = fPos(1) + fPos(3);
  set(hf, 'Position', fPos);
else
  hf = clf(figure(hf));
end
hax = axes(hf);
box(hax, 'on');
title(hax, 'Thermal frequency drift over loops');
ylabel(hax, 'Larmor frequency in MHz');
xlabel(hax, 'time in s');
set(get(hax, 'YRuler'), 'TickLabelFormat', '%10.6f');
grid(hax, 'on');
hold(hax, 'on');

% measurement loop
for loop = 1:Seq.loops
  if Seq.loops > 1
    disp(['Loops left: ' num2str(Seq.loops-loop+1)]);
    if loop == 1
      Seq.StartSequenceTime = now*24*3600 + Seq.LoopBreak + 1;
    else
      Seq.StartSequenceTime = SeqOut.StartSequenceTime + Seq.LoopBreak;
      if Seq.LoopBreak <= 2
        % Don't repeat hardware initialization in each iteration, see manual (HW)
        HW.tRepInit = 0.05;
      end
    end
  end

  % create pulse program and start measurement
  [data, SeqOut] = sequence_EchoStandard(HW, Seq);

  if SeqOut.loops > 1
    % generate LoopData structure
    if loop == 1
      clear LoopData;
      LoopData.fOffset = zeros(Seq.loops, 1);
      LoopData.StartTime = zeros(Seq.loops, 1);
      Seq.plotSeq = [];
    end
    Seq.plotAllHandle = SeqOut.plotAllHandle;
    Seq.plotRaiseWindow = false;

    % store LoopData
    LoopData.StartTime(loop) = SeqOut.StartSequenceTime;
    % calculate mean frequency offset at FID
    tAQEndIdx = find(data(1).time_all(1:SeqOut.AQ(1).nSamples(1),1,1)>SeqOut.fOffsetTime, 1, 'first');
    LoopData.fOffset(loop) = mean(diff(unwrap(angle(data(1).data(1:min([SeqOut.AQ(1).nSamples(1),tAQEndIdx]),1,1))))) * ...
      SeqOut.AQ(1).fSample(1)/2/pi;
    % compensate frequency drift for next iteration
    HW.fLarmor = HW.fLarmor - LoopData.fOffset(loop);
    LoopData.fLarmor(loop) = HW.fLarmor;

    % plot frequency drift
    if loop == 1
      hl = plot(hax, LoopData.StartTime(1:loop)-LoopData.StartTime(1), ...
        LoopData.fLarmor(1:loop)/1e6);
    else
      set(hl, 'XData', LoopData.StartTime(1:loop)-LoopData.StartTime(1), ...
        'YData', LoopData.fLarmor(1:loop)/1e6);
    end
  end
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
