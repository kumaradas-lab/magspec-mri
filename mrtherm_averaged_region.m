clear all;
close all;

%fName = 'exp1';
% --- Set up save directory in Google Drive ---
baseFolder = 'G:\My Drive\MR_Thermometry\My runs';             % main directory on Google Drive, can change this to match the path in your Google Drive
todayFolder = datestr(now, 'yyyy-mm-dd');           % folder named with today's date
saveDir = fullfile(baseFolder, todayFolder);         % full path for today's run

if ~exist(saveDir, 'dir')
   mkdir(saveDir);                                 % create folder if it doesn't exist
end

timestamp = datestr(now, 'HHMMSS');
fName = fullfile(saveDir, ['exp1_' timestamp]);


% Initialize Osensa temperature sensor
osensa_dev = enable_osensa("COM3");

% Initialize MRI system parameters
LoadSystem; % Load system parameters (reset to default: HW Seq AQ TX Grad)
% HW.fLarmor = 23.42e6;                        % set correct Larmor frequency (0.55 T)
% HW.FindFrequencySweep.fCenter = 23.42e6;     % center sweep near expected resonance
% HW.FindFrequencySweep.fRange  = 500e3;       % widen sweep to ±250 kHz just in case
% HW.FindFrequencySweep.fOffsetFIDsStdMaxValue = 5000;  % allow more noise tolerance

Seq.Loops = 1; % Number of loop averages

% Define parameters
Seq.T1 = 3000e-3;    %change for water 
Seq.tEcho = 3e-3; % try for 3, 5, 20
Seq.tRep = 200e-3;    % try higher to stabilize the phase
resolution = 32; % original 32x32
thickness = 0.002; % original 0.002
pausetime = 2;
position = resolution / 2;
measurement_time = 300; % Run time in seconds


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
%% Disable ALL internal plotting to prevent Pure Devices crashes
Seq.plotSeqAQ = 0;
Seq.LoopPlot = 0;
Seq.AQSlice(1).plotkSpace = 0;
Seq.AQSlice(1).plotImage = 0;
Seq.AQSlice(1).plotImageHandle = [];
Seq.AQSlice(1).plotPhase = 0;

% Extra safety flags (PD code checks these too)
Seq.plot = 0;
Seq.AQPlot = 0;
Seq.AQSlice(1).PlotkSpace = 0;
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;
Seq.AQSlice(1).ZeroFillFactor = 4;
Seq.AQSlice(1).ThicknessPos = [0 0 -0.01]; % position of the slice

Seq.CorrectSliceRephase = 0;                        % Correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % Correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;
Seq.CorrectPhaseDuration = 2.5e-3;

% Initialize data storage
i = 0;
tStart = tic;
roiSize = 3;
Timedata = [];
TemperatureData = [];
Phasedata = [];
Acquisitiondata = {};

while true
  i = i + 1;
  fprintf ("Acquisition of Imege %d\n", i)
  
  time = toc(tStart);
  Timedata(i) = time;
  
  % Read Osensa temperature sensor at each acquisition
  TemperatureData(i) = osensa_dev.read_channel_temp();
  
  % Run MRI acquisition sequence
  [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
  Acquisitiondata{i} = SeqLoop;
  
  % Trying using 3x3 ROI innstead of single pixel for phase
  
  x1 = position - floor(roiSize/2);
  x2 = position + floor(roiSize/2);
  roi = SeqLoop.data.Image (x1:x2, 1, x1:x2);   %3x3 ROI in central slice

  roi_mean_phase = angle(mean(roi(:)));
  Phasedata(i)=roi_mean_phase;
  
  % Stop acquisition if measurement time is exceeded
  if Timedata(i) > measurement_time
    break;
  end
end
%% ---- POST PROCESSING: unwrap and reference substraction
Phasedata_unwrapped = unwrap(Phasedata);
refIdx =5;
Referencephase = Phasedata_unwrapped(refIdx);

Deltaphase = Phasedata_unwrapped -Referencephase;

phase_diff_figure = figure ('Name', 'Last Image: Magnitude & Phase');
plot(Timedata, Deltaphase, '-o');
xlabel('Time(s)');
ylabel('Phase difference (rad)');
title('Phase Difference vs Time');
grid on;

% plot image and phase of the last image acquired
figure5 = figure('Name', 'Last Imgae: Magnitude & Phase');
subplot(1,2,1)
Imagemangnitude = squeeze(abs(SeqLoop.data.Image(:, 1, :)));  %changed by Isabella
imagesc(Imagemangnitude)
axis equal tight;
colorbar;
title('Magnitude of the Image')

subplot(1,2,2)
Imagephases = squeeze(angle(SeqLoop.data.Image(:, 1, :))); % changed by Isabella
imagesc(Imagephases)
axis equal tight;
colorbar;
title('Image Phasemap')
% Save full figure after both subplots are drawn
saveas(figure5, fullfile(saveDir, ['Image_Magnitude_Phase_' timestamp '.png']));


% Save all data
save([fName '.mat'], 'Timedata', 'TemperatureData', 'Phasedata', 'Phasedata_unwrapped', 'Deltaphase', 'Acquisitiondata');
csvwrite([fName '.csv'], [Timedata' TemperatureData' Phasedata' Deltaphase']);
saveas(phase_diff_figure, [fName '.png']);

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

saveas(gcf, fullfile(saveDir, ['Phase_vs_Temperature' timestamp '.png']));

save(fullfile(saveDir, ['Trial5_Res32_TR300_RoomTemp_oil_Apr9' timestamp '.mat']));