%% MRI Thermometry Acquisition Script 
clear; close all;

%% --- Set up save directory ---
baseFolder = 'G:\My Drive\MR_Thermometry\My runs';
todayFolder = datestr(now, 'yyyy-mm-dd');
saveDir = fullfile(baseFolder, todayFolder);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% --- Load System Parameters ---
osensa_dev = enable_osensa("COM3");
LoadSystem; 

%% --- Sequence parameters ---
Seq.Loops = 1;        % s
Seq.T1 = 3;           % T1 for water (s)
Seq.tEcho = 3e-3;     % TE (s)
Seq.tRep = 200e-3;    % TR (s)
Seq.CorrectPhaseDuration = 0.6e-3;  % s
resolution = 32;       
thickness = 0.002;            
position = resolution/2; 
measurement_time = 500; % s

%% --- Naming convention ---
dateStr = datestr (now, 'yyyymmdd');
TE_ms = round(Seq.tEcho*1e3);
TR_ms = round(Seq.tRep*1e3);
sample = 'water';                  %<-- change if needed
orientation ='zx';

% auto-increment run number
runNum = 1;
while true
    testName =sprintf('%s_%s_TE%dms_TR%dms_RES%d_%s_run%02d.mat',...
        dateStr, sample, TE_ms, TR_ms, resolution, orientation, runNum);
    if ~exist(fullfile(saveDir, testName), 'file')
        break;
    end
    runNum = runNum+1;
end 
runID = sprintf('run%02d', runNum);

baseName = sprintf('%s_%s_TE%dms_TR%dms_RES%d_%s_%s', ...
    dateStr, sample, TE_ms, TR_ms, resolution, orientation, runID);

fName = fullfile(saveDir, baseName);

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
%roiSize = 3;
%Seq.CorrectPhase =0;
%% --- Start acquisition ---
tStart = tic;
i = 0;

while true
    i = i + 1;
    fprintf('Acquisition of Image %d\n', i);
   
    Timedata(i) = toc(tStart);
   
    TemperatureData(i) = osensa_dev.read_channel_temp();  % Read Osensa temperature
    
    % Run MRI acquisition
    [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
    Acquisitiondata{i} = SeqLoop;
    
    % --- Extract ROI phase ---
    x1 = position - floor(roiSize/2);
    x2 = position + floor(roiSize/2);
    roi = SeqLoop.data.Image(x1:x2, 1, x1:x2);
    Phasedata(i) = angle(SeqLoop.data.Image(position, 1, position));
    
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
saveas(gcf, fullfile(saveDir, ['Phase_vs_Time_' baseName '.png']));

%% --- Plot phase difference over Osensa Temperature ---
figure('Name','Phase Difference vs Temperature');
plot(TemperatureData_p, Deltaphase_p, '-o','LineWidth',1.5);
xlabel('Temperature(°C)'); ylabel('Phase Difference (rad)');
title('Phase Difference vs Temperature'); grid on;

pT = polyfit(TemperatureData_p, Deltaphase_p, 1);
%yfitT = polyval(pT, TemperatureData_p);
hold on;
plot(TemperatureData_p, polyval(pT, TemperatureData_p), '--r');

legend('Data', sprintf('Fit: y = %.3fx + %.3f', pT(1), pT(2)));
saveas(gcf, fullfile(saveDir, ['Phase_vs_Temperature_' baseName '.png']));


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
saveas(gcf, fullfile(saveDir, ['LastImage_MagPhase_' baseName '.png']));

%% --- Save data ---
save([fName '.mat'], 'Timedata','TemperatureData','Phasedata','Phasedata_unwrapped','Deltaphase','Acquisitiondata');
writematrix([Timedata_p' TemperatureData_p' DeltaT_Osensa' DeltaT_MRI' Deltaphase_p'], ...
    [fName '_processed.csv']);

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

%% phase difference map
refIdx = 5;
lateIdx = length(Acquisitiondata);

refImage  = squeeze(Acquisitiondata{refIdx}.data.Image(:,1,:));
lateImage = squeeze(Acquisitiondata{lateIdx}.data.Image(:,1,:));

magRef  = abs(refImage);
phiRef  = angle(refImage);
phiLate = angle(lateImage);

DeltaPhi = angle(exp(1i*(phiLate - phiRef)));

figure('Name','Magnitude, Phase, and Phase Difference');

subplot(1,3,1);
imagesc(magRef); axis equal tight; colorbar;
title(sprintf('Magnitude (frame %d)', refIdx));

subplot(1,3,2);
imagesc(phiRef); axis equal tight; colorbar;
title(sprintf('Phase (frame %d)', refIdx));

subplot(1,3,3);
imagesc(DeltaPhi); axis equal tight; colorbar;
title('\Delta\phi = \phi_{late} - \phi_{ref}');

saveas(gcf, fullfile(saveDir, ['Mag_Phase_DeltaPhi_' baseName '.png']));



%% ---- Estimate alpha ----
gamma = 2*pi*42.58e6;   % rad/T/s
B0 = 0.55;              % Tesla
TE = Seq.tEcho;         % seconds

alpha_est = p(1) / (gamma * B0 * TE);   % fractional / °C
alpha_ppm = alpha_est * 1e6;            % / °C

% Compute MRI temperature change
DeltaT_MRI = Deltaphase_p ./ (gamma * alpha_used * B0 * TE);

% Reference Osensa the same way
ReferenceT = TemperatureData(refIdx);
DeltaT_Osensa = TemperatureData_p - ReferenceT;

figure;
plot(Timedata_p, DeltaT_Osensa, '-o','LineWidth',1.5); hold on;
plot(Timedata_p, DeltaT_MRI, '-s','LineWidth',1.5);

xlabel('Time (s)');
ylabel('\Delta Temperature (°C)');
title('MRI vs Osensa Temperature Change');
legend('Osensa','MRI');
grid on;

saveas(gcf, fullfile(saveDir, ['TempComparison_' baseName '.png']));

%% --- Correlation plot ---
figure;
plot(DeltaT_Osensa, DeltaT_MRI, 'o','LineWidth',1.5);
xlabel('Osensa \DeltaT (°C)');
ylabel('MRI \DeltaT (°C)');
title('MRI vs Osensa Temperature Correlation');
grid on;
hold on;

p_corr = polyfit(DeltaT_Osensa, DeltaT_MRI, 1);
plot(DeltaT_Osensa, polyval(p_corr, DeltaT_Osensa), '--r');

legend('Data', sprintf('y = %.2fx + %.2f', p_corr(1), p_corr(2)));

saveas(gcf, fullfile(saveDir, ['TempCorrelation_' baseName '.png']));

%% --- Final Delta T MRI computation ---
DeltaT_final_MRI = DeltaT_MRI(end);
DeltaT_final_Osensa = DeltaT_Osensa(end);

fprintf('\nFinal Temperature Change:\n');
fprintf('MRI: %.3f °C\n', DeltaT_final_MRI);
fprintf('Osensa: %.3f °C\n', DeltaT_final_Osensa);

%% ---- Display results ----
fprintf('Estimated PRF coefficient alpha:\n');
fprintf('  alpha = %.3e /°C (%.3f ppm/°C)\n', alpha_est, alpha_ppm);

legend('Data', ...
       sprintf('\\Delta\\phi = %.3f\\DeltaT + %.3f', p(1), p(2)), ...
       'Location','best');

saveas(gcf, fullfile(saveDir, ['DeltaPhase_vs_DeltaTemp_' timestamp '.png']));


%% --- Save EVERYTHING (workspace snapshot, safe for big files) ---
save([fName '_ALL.mat'], '-v7.3');
