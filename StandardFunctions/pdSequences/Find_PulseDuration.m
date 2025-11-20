function [HW, mySave] = Find_PulseDuration(HW, mySave, minTime, doPlot, iterations, tPulse90, T1, Seq)
%% Find 90 degrees pulse duration and store corresponding calibration value in file
%
%   [HW, mySave] = Find_PulseDuration(HW, mySave, minTime, doPlot, iterations, tPulse90, T1, Seq)
%
% This function searches the correct pulse length for a 90 degrees pulse at the
% default amplitude and stores it in HW and mySave. Additionally, it appends the
% corresponding factor for output voltage on the TX channel to B1 field strength
% to the PaUout2AmplitudeCal.m file in the User folder.
%
% It determines the pulse length by iteratively increasing the flip angle and
% performing a 1D Spin Echo measurement along the y-axis of the sample. The
% amplitude at the inner half of the acquired 1D images is used for calculating
% the 90 degrees pulse length (i.e. the efficiency of the coil).
%
% Consistency checks with the apparent 90 degrees and 180 degrees  pulse length
% at the center and the inner half of the acquired 1D images are performed to
% evaluate the quality of the results.
%
%
% INPUT:
%
% All input parameters but HW are optional. If they are omitted or empty,
% default values are used.
%
%   HW
%       HW structure or object
%
%   mySave
%       mySave structure (necessary for minTime)
%
%   minTime
%       Minimum time in seconds since the last time Find_PulseDuration was
%       successfully executed before a new value for the coil efficiency is
%       determined again (default: 1000).
%
%   doPlot
%       Plot sequence and data (bool, default: 1)
%
%   iterations
%       Number of repetitions of the pulse length search (default: 1).
%       If set to 0, the value given by tPulse90 (see below) is used to set the
%       coil efficiency.
%
%   tPulse90
%       Estimated pulse length for the 90 degrees pulse in seconds (default:
%       HW.FindPulseDuration.tPulse90Estimated)
%
%   T1
%       approximate T1 time in seconds of the used sample (default: 0.2)
%
%   Seq
%       structure with additional settings that are passed to
%       sequence_PulseDuration (unless they are overridden by the other
%       parameters). (Default: none)
%       Fields that are used in this function include:
%
%     skipFrequencySweep
%         Do not search magnet frequency before starting the measurement
%         (default: false).
%
% Further settings can be changed by defining custom default values in the
% structure HW.FindPulseDuration (see sequence_PulseDuration).
%
%
% OUTPUT:
%
%   HW
%       HW structure or object
%
%   mySave
%       mySave structure
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 2,                        mySave     = [];      end
if nargin < 3 || isempty(minTime),    minTime    = 1000;    end
if nargin < 4,                        doPlot     = [];      end
if nargin < 5 || isempty(iterations), iterations = 1;       end
if nargin < 6 || isempty(tPulse90),   tPulse90   = HW.FindPulseDuration.tPulse90Estimated; end
if nargin < 7,                        T1         = HW.FindPulseDuration.T1Estimated; end
if nargin < 8,                        Seq        = struct(); end

if isemptyfield(Seq, 'skipFrequencySweep')
  Seq.skipFrequencySweep = false;
end
if isemptyfield(Seq, 'useGammaX')
  Seq.useGammaX = false;
end

%% initialization
% first time clear mySave
if isempty(mySave)
  mySave.lastTime_PulseDuration = 0;
  mySave.HW = HW;
end

if isemptyfield(mySave, 'lastTime_PulseDuration'), mySave.lastTime_PulseDuration = 0; end

% If time since last frequency search has not yet exceeded "minTime", do not
% search again.
if (now*24*3600-mySave.lastTime_PulseDuration <= minTime) && (iterations ~= 0)
  return;
end

if ~Seq.skipFrequencySweep
  oldFindFrequencyPause = HW.FindFrequencyPause;
  hwGuard = onCleanup(@() setfield(HW, 'FindFrequencyPause', oldFindFrequencyPause));  %#ok<SFLD>
  HW.FindFrequencyPause = T1*3;
  % find magnet frequency
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1000);
  % find frequency more accurately
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0, 0, 1, HW.tFlip90Def*HW.FindFrequencySweep.tPulseDivider, 1, 512);
  delete(hwGuard);
end

Seq.tPulse90Estimated = tPulse90;
for t = 1:iterations
  disp('Searching RF pulse length...')
  Seq.T1Estimated = T1;
  if t > 1
    Seq.tPulse90Estimated = SeqLoop.dataPulseDuration.tFlip90;
  end
  Seq.doPlot = doPlot;
  AQ = [];
  TX = [];
  Grad = [];

  [SeqLoop, mySave] = sequence_PulseDuration(HW, Seq, AQ, TX, Grad, mySave);

  if HW.FindFrequencyPause > 0
    % Pause to relax spin system
    if HW.FindFrequencyPause > 3
      disp(['Waiting ' num2str(HW.FindFrequencyPause) ' seconds after finding pulse duration']);
    end
    sleep(HW.FindFrequencyPause);
  end
end


%% Store calculated values or display them
if iterations > 0
  data = SeqLoop.dataPulseDuration;
else
  % Set the tPulse90 without searching.
  data.tFlip90 = tPulse90;
  data.savePulseFile = true;
  % calculate magnetic flux density of the transmit coil at HW.TX.AmpDef amplitude
  data.B1 = ((pi/2)/data.tFlip90) / HW.GammaDef;

  if ~isemptyfield(Seq, 'AQSlice') && ~isemptyfield(Seq.AQSlice(1), 'iDevice')
    SeqLoop.AQSlice.iDevice = Seq.AQSlice(1).iDevice;
  else
    SeqLoop.AQSlice.iDevice = 1;
  end
end

% Store calculated values or display them
if Seq.useGammaX
  PaUout2AmplitudeFieldname = 'PaUout2AmplitudeX';
  fLarmorPulseDuration = HW.fLarmorX;
else
  PaUout2AmplitudeFieldname = 'PaUout2Amplitude';
  fLarmorPulseDuration = HW.fLarmor;
end

newPaUout2Amplitude = HW.TX(SeqLoop.AQSlice.iDevice).(PaUout2AmplitudeFieldname);
TXVoltage = HW.TX(SeqLoop.AQSlice.iDevice).AmpDef / HW.TX(SeqLoop.AQSlice.iDevice).PaUout2Amplitude(HW.TX(SeqLoop.AQSlice.iDevice).ChannelDef);
newPaUout2Amplitude(HW.TX(SeqLoop.AQSlice.iDevice).ChannelDef) = data.B1 / TXVoltage;
comment = sprintf('%s (tFlip90 = %.3f us @ %.3f V @ %.6f MHz) from 1d Spin Echo by %s', ...
  datestr(now, 'yyyy-mm-ddTHH:MM:SS'), data.tFlip90*1e6, TXVoltage, fLarmorPulseDuration/1e6, mfilename());

newCalLine = sprintf('HW.TX(%d).%s = [%.6f, %.6f]*1e-6;', ...
  SeqLoop.AQSlice.iDevice, PaUout2AmplitudeFieldname, newPaUout2Amplitude*1e6);

if ~isempty(HW.TX(SeqLoop.AQSlice.iDevice).CoilName)
  newCalLine = sprintf('if strcmp(HW.TX(%d).CoilName, ''%s''),  %s  end', ...
    SeqLoop.AQSlice.iDevice, HW.TX(SeqLoop.AQSlice.iDevice).CoilName, newCalLine);
end

if SeqLoop.AQSlice.iDevice > 1
  newCalLine = sprintf('if numel(HW.TX) >= %d,  %s  end', ...
    SeqLoop.AQSlice.iDevice, newCalLine);
end

newCalLine = sprintf('%s  %% %s\n', newCalLine, comment);

if ~isempty(HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath) ...
    && (isemptyfield(mySave, 'DummySerial') ...
        || mySave.DummySerial(min(SeqLoop.AQSlice.iDevice, numel(mySave.DummySerial))) <= 0)
  addFirstLine = ~exist(HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath, 'file');
  if ~exist(fileparts(HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath), 'dir')
    mkdir(fileparts(HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath));
  end
  fid = fopen(HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath, 'a+');
  fid_protect = onCleanup(@() fclose(fid));
  if addFirstLine
    fwrite(fid, ['% factor from voltage amplitude at the coil input to B1+ field strength in T/V', sprintf('\n')]);
  end
  if data.savePulseFile
    fwrite(fid, newCalLine);
    fprintf('\nA new line was added to the following file:\n%s\n%s\n', ...
      HW.TX(SeqLoop.AQSlice.iDevice).PaUout2AmplitudePath, newCalLine);
  else
    fwrite(fid, ['% ', newCalLine]);  % add line as comment
  end
  delete(fid_protect);
  [~, name, ~] = fileparts(fopen(fid));
  clear(name);clear('name')

  % save the time of the last RF pulse duration search
  mySave.lastTime_PulseDuration = now*24*3600;
elseif data.savePulseFile
  fprintf('\n');
  fprintf('\nPlease add the following line to your LoadMySystem.m file:\n%s\n', ...
    newCalLine);
end
if ~data.savePulseFile
  fprintf('\n');
  warnStr = ['Determination of pulse length unsuccessful!\n', ...
    'If you want to use the uncertain best guess value anyway, ', ...
    'please manually append or un-comment the following line in your PaUout2AmplitudeCal.m file:\n%s\n'];
  warning('PD:sequence_PulseDuration', warnStr, newCalLine);
end

fprintf('90%s pulse duration: %.3f %ss @ %.3f V\n', char(176), data.tFlip90*1e6, char(181), TXVoltage);

end
