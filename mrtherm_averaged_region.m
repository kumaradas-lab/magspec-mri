clear all;
close all;

fName = 'exp1';

% --- Set up save directory in Google Drive ---
baseFolder = 'G:\My Drive\MRI\My runs';             % main directory on Google Drive, can change this to match the path in your Google Drive
todayFolder = datestr(now, 'yyyy-mm-dd');           % folder named with today's date
saveDir = fullfile(baseFolder, todayFolder);         % full path for today's run

if ~exist(saveDir, 'dir')
    mkdir(saveDir);                                 % create folder if it doesn't exist
end

fName = fullfile(saveDir, 'exp1');                  % full base name for saving files


% Initialize Osensa temperature sensor
osensa_dev = enable_osensa("COM3");                  % Change when setting up

% Initialize MRI system parameters
LoadSystem; % Load system parameters (reset to default: HW Seq AQ TX Grad)
HW.fLarmor = 23.42e6;                        % set correct Larmor frequency (0.55 T)
HW.FindFrequencySweep.fCenter = HW.fLarmor;     % center sweep near expected resonance
HW.FindFrequencySweep.fRange  = 500e3;       % widen sweep to ±250 kHz just in case, can increase the range depending on the medium
HW.FindFrequencySweep.fOffsetFIDsStdMaxValue = 5000;  % allow more noise tolerance

Seq.Loops = 1; % Number of loop averages

% Define parameters
Seq.T1 = 100e-3;
Seq.tEcho = 3e-3; % try for 3, 5, 20
Seq.tRep = 8e-3;
resolution = 32; % original 32x32
thickness = 0.002; % original 0.002
pausetime = 2;
position = resolution / 2;
measurement_time = 60; % Run time in seconds

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = resolution;
Seq.AQSlice(1).nPhase(2) = resolution;
Seq.AQSlice(1).HzPerPixMin = 0;
Seq.AQSlice(1).sizeRead = 0.010;
Seq.AQSlice(1).sizePhase(2) = 0.010;
Seq.AQSlice(1).thickness = thickness;
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orientation = 'zx';                                 % 'xy', 'yz', 'zx' for one of the cardinal planes (read-phase)
switch orientation
  case 'xy'
  Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
  Seq.AQSlice(1).phi  = 0.5*pi;                   % 2nd rotation around y axis in RAD
  Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
  case 'yz'
  Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
  Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
  Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
  case 'zx'
  Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
  Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
  Seq.AQSlice(1).theta= -0.5*pi;                  % 3rd rotation around z axis in RAD
  otherwise
  if ~ischar(orientation)
    orientation = num2str(orientation);
  end
  error('Unknown orientation "%s"\n', orientation);
end

%% Set up sequence visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;
Seq.LoopPlot = 1;
Seq.AQSlice(1).plotkSpace = 1;
Seq.AQSlice(1).plotImage = 1;
Seq.AQSlice(1).plotPhase = 1;
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;
Seq.AQSlice(1).ZeroFillFactor = 4;
Seq.AQSlice(1).ThicknessPos = [0 0 -0.01]; % position of the slice

Seq.CorrectSliceRephase = 0;                        % Correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % Correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;

% Initialize data storage
i = 0;
tStart = tic;
Referencephase = 0;
phase_diff_figure = figure('Name', 'Phase Difference vs. Time');

while true
  i = i + 1;
  disp("Acquisition of Image " + num2str(i));
  
  time = toc(tStart);
  Timedata(i) = time;
  
  % Read Osensa temperature sensor at each acquisition
  TemperatureData(i) = osensa_dev.read_channel_temp();
  
  % Run MRI acquisition sequence
  [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
  Acquisitiondata(i) = SeqLoop;
  
  % Extract phase data from acquired image
  Imagephase = unwrap(angle(SeqLoop.data.Image(position, 1, position)));
  Phasedata(i) = Imagephase;
  
  % Compute phase difference
  if i < 4
    deltaphase = 0;
  elseif i == 5
    Referencephase = Phasedata(5);
    deltaphase = 0;
  else
    deltaphase = Imagephase - Referencephase;
  end
  
  Deltaphase(i) = deltaphase;
  
  % Plot phase difference over time
  figure(phase_diff_figure);
  hold on;
  plot(gca, time, deltaphase, '-o');
  xlabel('Time (s)');
  ylabel('Phase difference (rad)');
  title('Phase difference Vs Time');
  
  pause(2); % Pausing for next acquisition
  
  % Stop acquisition if measurement time is exceeded
  if Timedata(i) > measurement_time
    break;
  end
end

% plot image and phase of the last image acquired
figure5 = figure('Name', 'Delta Phase over time');
figure(figure5)
subplot(1,2,1)

% Imagemangnitude = abs(reshape(SeqLoop.data.Image(:,1,:), [[],32]));
% % changed by Isabella
Imagemangnitude = squeeze(abs(SeqLoop.data.Image(:, 1, :)));
imagesc(Imagemangnitude)
colorbar
axis equal
title('Magnitude of the Image')

subplot(1,2,2)
% Imagephases = angle(reshape(SeqLoop.data.Image(:,1,:), [[],32]));
% % changed by Isabella
Imagephases = squeeze(angle(SeqLoop.data.Image(:, 1, :)));
imagesc(Imagephases)
colorbar 
axis equal
title('Image Phasemap')
% Save full figure after both subplots are drawn
saveas(figure5, fullfile(saveDir, 'Image_Magnitude_Phase.png'));

% Save all data
save(fullfile(saveDir, [fName '.mat']), 'Timedata', 'TemperatureData', 'Phasedata', 'Deltaphase', 'Acquisitiondata');
csvwrite([fName '.csv'], [Timedata' TemperatureData' Phasedata' Deltaphase']);
saveas(phase_diff_figure, fullfile(saveDir, [fName '.png']));

% Close Osensa sensor
osensa_dev.close();
disp("Osensa Transmitter OFF");

% Plot delta phase vs temperature
figure;
plot(TemperatureData, Deltaphase, 'o-','LineWidth',1.5);
xlabel('Temperature (°C)');
ylabel('Phase Difference (rad)');
title('Phase Difference vs Temperature');
grid on;
p = polyfit(TemperatureData, Deltaphase, 1);
yfit = polyval(p, TemperatureData);
hold on;
plot(TemperatureData, yfit, '--r');
legend('Data', sprintf('Fit: y = %.3fx + %.3f', p(1), p(2)));
saveas(gcf, fullfile(saveDir, 'Phase_vs_Temperature.png'));

save(fullfile(saveDir,'Trial5_Res32_TR300_RoomTemp_oil_Apr9.mat'));