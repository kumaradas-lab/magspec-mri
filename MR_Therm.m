%% MRI Thermometry Acquisition Script (Clean Version)
clear; close all;

%% --- Set up save directory ---
baseFolder = 'G:\My Drive\MR_Thermometry\My runs';
todayFolder = datestr(now, 'yyyy-mm-dd');
saveDir = fullfile(baseFolder, todayFolder);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
timestamp = datestr(now, 'HHMMSS');
fName = fullfile(saveDir, ['exp1_' timestamp]);

osensa_dev = enable_osensa("COM3");

LoadSystem; 

%% --- Sequence parameters ---
Seq.Loops = 1;  
Seq.T1 = 3;          % T1 for water (s)
Seq.tEcho = 3e-3;     % TE (s)
Seq.tRep = 200e-3;     % TR (s)
Seq.CorrectPhaseDuration = 0.6e-3;  % s
resolution = 32;       
thickness = 0.002;     
%pausetime = 2;         
position = resolution/2; 
measurement_time = 420; % s

%% --- Acquisition parameters ---
Seq.AQSlice(1).nRead = resolution;
Seq.AQSlice(1).nPhase(2) = resolution;
Seq.AQSlice(1).sizeRead = 0.010;
Seq.AQSlice(1).sizePhase(2) = 0.010;
Seq.AQSlice(1).thickness = thickness;
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;
Seq.AQSlice(1).PhaseOS(2) = 2;  % phase oversampling
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;
Seq.AQSlice(1).ZeroFillFactor = 4;
Seq.AQSlice(1).ThicknessPos = [0 0 -0.01];

%% --- Orientation ---
orientation = 'zx';
switch orientation
    case 'xy'
        Seq.AQSlice(1).alfa = 0; Seq.AQSlice(1).phi = 0.5*pi; Seq.AQSlice(1).theta = 0;
    case 'yz'
        Seq.AQSlice(1).alfa = 0; Seq.AQSlice(1).phi = 0; Seq.AQSlice(1).theta = 0;
    case 'zx'
        Seq.AQSlice(1).alfa = 0; Seq.AQSlice(1).phi = 0; Seq.AQSlice(1).theta = -0.5*pi;
    otherwise
        error('Unknown orientation "%s"', orientation);
end

%% --- Visualization disabled for stability ---
Seq.plot = 0; Seq.AQPlot = 0;
Seq.AQSlice(1).plotkSpace = 0;
Seq.AQSlice(1).plotImage = 0;
Seq.AQSlice(1).plotPhase = 0;

%% --- Data storage ---
Timedata = [];
TemperatureData = [];
Phasedata = [];
Acquisitiondata = {};
roiSize = 3;
%Seq.CorrectPhase =0;
%% --- Start acquisition ---
tStart = tic;
i = 0;

while true
    i = i + 1;
    fprintf('Acquisition of Image %d\n', i);
    
    time = toc(tStart);
    Timedata(i) = time;
    
    % Read Osensa temperature
    TemperatureData(i) = osensa_dev.read_channel_temp();
    
    % Run MRI acquisition
    [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
    Acquisitiondata{i} = SeqLoop;
    
    % --- Extract ROI phase ---
    x1 = position - floor(roiSize/2);
    x2 = position + floor(roiSize/2);
    roi = SeqLoop.data.Image(x1:x2, 1, x1:x2);
    roi_mean_phase = angle(mean(roi(:)));
    Phasedata(i) = roi_mean_phase;
    
    % Stop acquisition if measurement time exceeded
    if Timedata(i) > measurement_time
        break;
    end
end

%% --- Post-processing: unwrap & reference subtraction ---
Phasedata_unwrapped = unwrap(Phasedata); % unwrap entire series
refIdx = 5; % baseline image (ignore first few transient images)
Referencephase = Phasedata_unwrapped(refIdx);
Deltaphase = Phasedata_unwrapped - Referencephase;
idx = refIdx:length(Timedata);

Timedata_p = Timedata(idx);
TemperatureData_p = TemperatureData(idx);
Deltaphase_p = Deltaphase(idx);


%% --- Plot phase difference over time ---
figure('Name','Phase Difference vs Time');
plot(Timedata_p, Deltaphase_p, '-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Phase Difference (rad)');
title('Phase Difference vs Time'); grid on;
saveas(gcf, fullfile(saveDir, ['Phase_vs_Time_' timestamp '.png']));

%% --- Plot phase difference over Osensa Temperature ---
figure('Name','Phase Difference vs Temperature');
plot(TemperatureData_p, Deltaphase_p, '-o','LineWidth',1.5);
xlabel('Temperature(°C)'); ylabel('Phase Difference (rad)');
title('Phase Difference vs Temperature'); grid on;
pT = polyfit(TemperatureData_p, Deltaphase_p, 1);
yfitT = polyval(pT, TemperatureData_p);
hold on;
plot(TemperatureData_p, yfitT, '--r');

legend('Data', sprintf('Fit: y = %.3fx + %.3f', pT(1), pT(2)));
saveas(gcf, fullfile(saveDir, ['Phase_vs_Temperature_' timestamp '.png']));


%% --- Plot last acquired image ---
figure('Name','Last Image Magnitude & Phase');
subplot(1,2,1);
lastAcquisition = Acquisitiondata{end};
Imagemangnitude = squeeze(abs(lastAcquisition.data.Image(:, 1, :)));
imagesc(Imagemangnitude); axis equal tight; colorbar;
title('Magnitude');
subplot(1,2,2);
Imagephase = squeeze(angle(lastAcquisition.data.Image(:,1,:)));
imagesc(Imagephase); axis equal tight; colorbar;
title('Phase');
saveas(gcf, fullfile(saveDir, ['LastImage_MagPhase_' timestamp '.png']));

%% --- Save data ---
save([fName '.mat'], 'Timedata','TemperatureData','Phasedata','Phasedata_unwrapped','Deltaphase','Acquisitiondata');
csvwrite([fName '.csv'], [Timedata' TemperatureData' Phasedata' Deltaphase']);

%% --- Close Osensa ---
osensa_dev.close();
disp('Osensa Transmitter OFF');

%% ---- ΔPhase vs ΔTemperature (PRF thermometry + alpha estimate) ----
ReferenceT = TemperatureData(refIdx);
DeltaTemp = TemperatureData - ReferenceT;
DeltaTemp_p = DeltaTemp(idx);
Deltaphase_p = Deltaphase(idx);

figure;
plot(DeltaTemp_p, Deltaphase_p, 'o','LineWidth',1.5);
xlabel('\Delta Temperature (°C)');
ylabel('\Delta Phase (rad)');
title('PRF Thermometry: \Delta\phi vs \DeltaT');
grid on;
hold on;

% Linear fit: slope = rad / °C
p = polyfit(DeltaTemp_p, Deltaphase_p, 1);
plot(DeltaTemp_p, polyval(p, DeltaTemp_p), '--r');

%% ---- Estimate alpha ----
gamma = 2*pi*42.58e6;   % rad/T/s
B0 = 0.55;              % Tesla
TE = Seq.tEcho;         % seconds

alpha_est = p(1) / (gamma * B0 * TE);   % fractional / °C
alpha_ppm = alpha_est * 1e6;            % ppm / °C

%% ---- Display results ----
fprintf('Estimated PRF coefficient alpha:\n');
fprintf('  alpha = %.3e /°C (%.3f ppm/°C)\n', alpha_est, alpha_ppm);

legend('Data', ...
       sprintf('\\Delta\\phi = %.3f\\DeltaT + %.3f', p(1), p(2)), ...
       'Location','best');

saveas(gcf, fullfile(saveDir, ['DeltaPhase_vs_DeltaTemp_' timestamp '.png']));


%% --- Save EVERYTHING (workspace snapshot, safe for big files) ---
save([fName '_ALL.mat'], '-v7.3');
