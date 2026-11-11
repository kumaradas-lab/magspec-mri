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
Seq.tEcho = 10e-3;     % TE (s)
Seq.tRep = 150e-3;    % TR (s)
Seq.CorrectPhaseDuration = 7e-3;  % try 0.6 for TE 3ms and 13 for TE 15ms, 15 for TE 18ms

resolution = 32;       
thickness = 0.002;            
position = resolution/2; 
measurement_time = 200; % s

%% --- Naming convention ---
dateStr = datestr(now, 'yyyymmdd');
TE_ms = round(Seq.tEcho*1e3);
TR_ms = round(Seq.tRep*1e3);
sample = 'water';                  %<-- change if needed
orientation = 'zx';

% auto-increment run number
runNum = 1;
while true
    testName = sprintf('%s_%s_TE%dms_TR%dms_RES%d_%s_run%02d.mat',...
        dateStr, sample, TE_ms, TR_ms, resolution, orientation, runNum);
    if ~exist(fullfile(saveDir, testName), 'file')
        break;
    end
    runNum = runNum + 1;
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
        Seq.AQSlice(1).alfa = 0; 
        Seq.AQSlice(1).phi = 0.5*pi; 
        Seq.AQSlice(1).theta = 0;
        
    case 'yz'
        Seq.AQSlice(1).alfa = 0; 
        Seq.AQSlice(1).phi = 0; 
        Seq.AQSlice(1).theta = 0;
        
    case 'zx'
        Seq.AQSlice(1).alfa = 0; 
        Seq.AQSlice(1).phi = 0; 
        Seq.AQSlice(1).theta = -0.5*pi;
        
    otherwise
        error('Unknown orientation "%s"', orientation);
end

%% --- Visualization disabled for stability ---
Seq.plot = 0; 
Seq.AQPlot = 0;
Seq.AQSlice(1).plotkSpace = 0;
Seq.AQSlice(1).plotImage = 0;
Seq.AQSlice(1).plotPhase = 0;

%% --- Data storage ---
Timedata = [];
TemperatureData = [];
Phasedata = [];
Acquisitiondata = {};

roiSize = 9;

% Store frequency correction information
FrequencyCorrectionData = {};

%% --- Start acquisition ---
tStart = tic;
i = 0;

[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);
HW.FindFrequencySweep.maxTime = 600;

while true
    
    i = i + 1;
    
    fprintf('Acquisition of Image %d\n', i);
   
    Timedata(i) = toc(tStart);
   
    % Read Osensa temperature
    TemperatureData(i) = osensa_dev.read_channel_temp();
    
    % Run MRI acquisition
    [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
    if i == 1

    disp('--- SeqLoop fields ---');
    disp(fieldnames(SeqLoop));

    if isfield(SeqLoop,'data')
        disp('--- SeqLoop.data fields ---');
        disp(fieldnames(SeqLoop.data));
    end

end
    Acquisitiondata{i} = SeqLoop;
    
    %% --- Extract ROI phase ---
    
    x1 = position - floor(roiSize/2);
    x2 = position + floor(roiSize/2);
    
    roi = SeqLoop.data.Image(x1:x2, x1:x2);
    
    % Store Full Complex ROI over time
    if i == 1
        roi_complex_time = zeros([size(roi), 1]);
    else
        roi_complex_time(:,:,i) = roi;
    end
    
    if i == 1
        roi_complex_time(:,:,1) = roi;
    end

    % Also store single pixel (for comparison)
    Phasedata(i) = angle(SeqLoop.data.Image(position, position));
    
    %% --- Store frequency correction information ---
    %
    % IMPORTANT:
    % sequence_Flash returns SeqLoop and mySave here, not SeqOut.
    % Therefore we cannot directly use SeqOut.Correct_tfoffset
    % unless those fields are stored inside SeqLoop.
    
    if isfield(SeqLoop, 'Correct_tfoffset')
        FrequencyCorrectionData{i}.time = SeqLoop.Correct_tfoffset;
    elseif isfield(SeqLoop.data, 'Correct_tfoffset')
        FrequencyCorrectionData{i}.time = SeqLoop.data.Correct_tfoffset;
    else
        FrequencyCorrectionData{i}.time = [];
    end
    
    if isfield(SeqLoop, 'Correct_foffset')
        FrequencyCorrectionData{i}.offset = SeqLoop.Correct_foffset;
    elseif isfield(SeqLoop.data, 'Correct_foffset')
        FrequencyCorrectionData{i}.offset = SeqLoop.data.Correct_foffset;
    else
        FrequencyCorrectionData{i}.offset = [];
    end
    
    if isfield(SeqLoop, 'Correct_foffsetStd')
        FrequencyCorrectionData{i}.std = SeqLoop.Correct_foffsetStd;
    elseif isfield(SeqLoop.data, 'Correct_foffsetStd')
        FrequencyCorrectionData{i}.std = SeqLoop.data.Correct_foffsetStd;
    else
        FrequencyCorrectionData{i}.std = [];
    end

    %% --- Stop acquisition if measurement time exceeded ---
    if Timedata(i) > measurement_time
        break;
    end
end

%% --- Post-processing: unwrap & reference subtraction ---

Phasedata_unwrapped = unwrap(Phasedata);  % single pixel

roi_complex_mean = squeeze(mean(mean(roi_complex_time,1),2)); 

roi_phase_mean = angle(roi_complex_mean);
roi_phase_mean_unwrapped = unwrap(roi_phase_mean);

refIdx = 5; % baseline image

%--- single pixel ---
Referencephase = Phasedata_unwrapped(refIdx);

Deltaphase = Phasedata_unwrapped - Referencephase;

%--- ROI per pixel ---
roi_ref = roi_phase_mean_unwrapped(refIdx);

Deltaphase_roi_mean = roi_phase_mean_unwrapped - roi_ref;

disp(std(Deltaphase_roi_mean))

idx = refIdx:length(Timedata);

Timedata_p = Timedata(idx);
TemperatureData_p = TemperatureData(idx);
Phasedata_p = Phasedata_unwrapped(idx);   %new
Deltaphase_p = Deltaphase(idx);
Deltaphase_roi_mean_p = Deltaphase_roi_mean(idx);

%% --- Plot phase difference over time ---

figure('Name','Phase Difference vs Time');

plot(Timedata_p, Deltaphase_p, '-o','LineWidth',1.5);

xlabel('Time (s)');
ylabel('Phase Difference (rad)');
title('Phase Difference vs Time');
grid on;

saveas(gcf, fullfile(saveDir, ...
    ['Phase_vs_Time_' baseName '.png']));

%% --- Plot phase difference vs Osensa Temperature ---

figure('Name','Phase Difference vs Temperature');

plot(TemperatureData_p, Deltaphase_p, '-o','LineWidth',1.5);

xlabel('Temperature(°C)');
ylabel('Phase Difference (rad)');
title('Phase Difference vs Temperature');
grid on;

pT = polyfit(TemperatureData_p, Deltaphase_p, 1);

hold on;

plot(TemperatureData_p, ...
    polyval(pT, TemperatureData_p), '--r');

legend('Data', ...
    sprintf('Fit: y = %.3fx + %.3f', pT(1), pT(2)));

saveas(gcf, fullfile(saveDir, ...
    ['Phase_vs_Temperature_' baseName '.png']));

%% --- Plot last acquired image ---

figure('Name','Last Image Magnitude & Phase');

subplot(1,2,1);

lastAcquisition = Acquisitiondata{end};

Imagemangnitude = squeeze(abs(lastAcquisition.data.Image(:, 1, :)));

imagesc(Imagemangnitude);
axis equal tight;
colorbar;

title('Magnitude');

subplot(1,2,2);

Imagephase = squeeze(angle(lastAcquisition.data.Image(:,1,:)));

imagesc(Imagephase);
axis equal tight;
colorbar;

title('Phase');

saveas(gcf, fullfile(saveDir, ...
    ['LastImage_MagPhase_' baseName '.png']));

%% ---- Estimate alpha & constants ----

gamma_val = 2*pi*42.58e6;   % rad/T/s
B0 = 0.55;                  % Tesla
TE = Seq.tEcho;             % seconds

% Linear fit for PRF coefficient
DeltaTemp = TemperatureData_p - TemperatureData_p(1);

p = polyfit(DeltaTemp, Deltaphase_p, 1);

alpha_est = p(1) / (gamma_val * B0 * TE);   % fractional / °C

alpha_ppm = alpha_est * 1e6;

alpha_used = alpha_est;   % ensure consistency

%% --- Final Delta T MRI computation ---

refImage = squeeze(Acquisitiondata{refIdx}.data.Image(:,1,:));

lateImage = squeeze(Acquisitiondata{end}.data.Image(:,1,:));

phiRef = angle(refImage);
phiLate = angle(lateImage);

DeltaPhi_pixel = Deltaphase(end);

DeltaPhi_roi = Deltaphase_roi_mean(end);

DeltaT_final_MRI_pixel = ...
    DeltaPhi_pixel / ...
    (gamma_val * alpha_used * B0 * TE);

DeltaT_final_MRI_mean = ...
    DeltaPhi_roi / ...
    (gamma_val * alpha_used * B0 * TE);

% Osensa ΔT
DeltaT_final_Osensa = ...
    TemperatureData(end) - TemperatureData(refIdx);

fprintf('Echo time: %.5f s\n', Seq.tEcho);

fprintf('Repetition time: %.5f s\n', Seq.tRep);

fprintf('Correct Phase duration: %.5f s', ...
    Seq.CorrectPhaseDuration);

fprintf('\nFinal Temperature Change:\n');

fprintf('MRI (pixel): %.5f °C\n', ...
    DeltaT_final_MRI_pixel);

fprintf('MRI (mean ROI): %.5f °C\n', ...
    DeltaT_final_MRI_mean);

fprintf('Osensa: %.5f °C\n', ...
    DeltaT_final_Osensa);

%% --- Plot MRI vs Osensa ΔT ---

DeltaT_MRI_pixel = ...
    Deltaphase_p ./ ...
    (gamma_val * alpha_used * B0 * TE);

DeltaT_MRI_roi = ...
    Deltaphase_roi_mean(idx) ./ ...
    (gamma_val * alpha_used * B0 * TE);

DeltaT_Osensa = ...
    TemperatureData_p - TemperatureData_p(1);

fprintf('\nDEBUG:\n');

fprintf('DeltaPhi_pixel (last): %.4f rad\n', ...
    DeltaPhi_pixel);

fprintf('Code DeltaT: %.4f °C\n', ...
    DeltaT_final_MRI_pixel);

figure;

plot(Timedata_p, DeltaT_Osensa, '-o','LineWidth',1.5);
hold on;

plot(Timedata_p, DeltaT_MRI_pixel, '-s','LineWidth',1.5);

plot(Timedata_p, DeltaT_MRI_roi, '-d','LineWidth',1.5);

legend('Osensa', 'MRI pixel', 'MRI ROI');

xlabel('Time (s)');
ylabel('\Delta Temperature (°C)');

title('MRI vs Osensa Temperature Change');

legend('Osensa','MRI pixel', 'MRI ROI');

grid on;

saveas(gcf, fullfile(saveDir, ...
    ['TempComparison_' baseName '.png']));

%% --- Correlation plot ---

figure;

plot(DeltaT_Osensa, DeltaT_MRI_pixel, ...
    'o','LineWidth',1.5);

xlabel('Osensa \DeltaT (°C)');
ylabel('MRI \DeltaT (°C)');

title('MRI pixel vs Osensa Temperature Correlation');

grid on;
hold on;

p_corr = polyfit(DeltaT_Osensa, ...
    DeltaT_MRI_pixel, 1);

plot(DeltaT_Osensa, ...
    polyval(p_corr, DeltaT_Osensa), '--r');

legend('Data', ...
    sprintf('y = %.2fx + %.2f', ...
    p_corr(1), p_corr(2)));

saveas(gcf, fullfile(saveDir, ...
    ['TempCorrelation(pixel)_' baseName '.png']));

%% --- Display PRF coefficient ---

fprintf('Estimated PRF coefficient alpha:\n');

fprintf('  alpha = %.3e /°C (%.3f ppm/°C)\n', ...
    alpha_est, alpha_ppm);

%% --- Data quality checks ---

% Expected PRF coefficient for sanity check
alpha_expected = -0.010e-6;

% Expected phase slope at this TE
expected_phase_slope = ...
    gamma_val * B0 * Seq.tEcho * alpha_expected;

fprintf('\n--- DATA QUALITY CHECK ---\n');

fprintf('Expected PRF phase slope: %.6f rad/C\n', ...
    expected_phase_slope);

fprintf('Measured phase change: %.6f rad\n', ...
    Deltaphase(end));

fprintf('Osensa temperature change: %.6f C\n', ...
    DeltaT_final_Osensa);

expected_phase_change = ...
    expected_phase_slope * DeltaT_final_Osensa;

fprintf('Expected phase change from Osensa: %.6f rad\n', ...
    expected_phase_change);

if abs(expected_phase_change) > eps
    
    fprintf('Measured/expected phase ratio: %.3f\n', ...
        Deltaphase(end) / expected_phase_change);
        
else
    
    fprintf('Measured/expected phase ratio: undefined (Osensa ΔT ≈ 0)\n');
    
end

%% --- ROI Signal Quality Analysis ---

ROI_Magnitude = squeeze(mean(mean(abs(roi_complex_time),1),2));

ROI_Complex_Mean = squeeze(mean(mean(roi_complex_time,1),2));

ROI_Coherence = abs(ROI_Complex_Mean) ./ ...
    (ROI_Magnitude + eps);

figure;

subplot(2,1,1);

plot(Timedata, ROI_Magnitude, '-o');

ylabel('ROI Magnitude');
title('ROI Signal Stability');
grid on;

subplot(2,1,2);

plot(Timedata, ROI_Coherence, '-o');

xlabel('Time (s)');
ylabel('Phase Coherence');
ylim([0 1.05]);
grid on;

%% --- Frequency correction information ---

% Use the last acquisition for diagnostics

frequencyTime = FrequencyCorrectionData{end}.time;
frequencyOffset = FrequencyCorrectionData{end}.offset;
frequencyStd = FrequencyCorrectionData{end}.std;

if ~isempty(frequencyTime) && ...
        ~isempty(frequencyOffset) && ...
        ~isempty(frequencyStd)
    
    fprintf('\n--- FREQUENCY CORRECTION DATA ---\n');
    
    fprintf('Frequency tracking data successfully retrieved.\n');
    
    fprintf('Frequency offset mean = %.6f Hz\n', ...
        mean(frequencyOffset(:)));
    
    fprintf('Frequency offset std = %.6f Hz\n', ...
        std(frequencyOffset(:)));
    
    fprintf('Frequency tracking standard error mean = %.6f Hz\n', ...
        mean(frequencyStd(:)));
    
    fprintf('Frequency tracking standard error max = %.6f Hz\n', ...
        max(frequencyStd(:)));
    
    figure;
    
    errorbar(frequencyTime, ...
        frequencyOffset, ...
        frequencyStd);
    
    xlabel('Time (s)');
    ylabel('Frequency offset (Hz)');
    title('Frequency Tracking');
    grid on;
    
    saveas(gcf, fullfile(saveDir, ...
        ['FrequencyTracking_' baseName '.png']));
        
else
    
    fprintf('\n--- FREQUENCY CORRECTION DATA ---\n');
    fprintf(['Frequency correction fields were not available ', ...
        'inside SeqLoop.\n']);
    
    fprintf(['The warning printed by get_CorrectedPhase is still ', ...
        'valid and should be treated as a data-quality warning.\n']);
    
end

%% --- Close Osensa ---

osensa_dev.close();

disp('Osensa Transmitter OFF');



%% --- Frequency Drift Analysis ---

nImages = length(FrequencyCorrectionData);

FrequencyMean = nan(nImages,1);
FrequencySpread = nan(nImages,1);
FrequencyUncertainty = nan(nImages,1);

for k = 1:nImages

    f = FrequencyCorrectionData{k}.offset;
    s = FrequencyCorrectionData{k}.std;

    if ~isempty(f)
        FrequencyMean(k) = mean(f(:));
        FrequencySpread(k) = std(f(:));
    end

    if ~isempty(s)
        FrequencyUncertainty(k) = mean(s(:));
    end

end

%% --- Match all data lengths ---

Timedata_col = Timedata(:);
FrequencyMean = FrequencyMean(:);

n = min(numel(Timedata_col), numel(FrequencyMean));

t = Timedata_col(1:n);
f = FrequencyMean(1:n);

%% --- Remove invalid measurements ---

valid = isfinite(t) & isfinite(f);

t = t(valid);
f = f(valid);

%% --- Frequency drift relative to baseline ---

refIdx = 5;

% Find the corresponding reference measurement
refTime = Timedata_col(refIdx);

[~, refFrequencyIdx] = min(abs(t-refTime));

DeltaFrequency = f - f(refFrequencyIdx);

%% --- Plot frequency drift ---

figure;

plot(t, DeltaFrequency, '-o', ...
    'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('\Delta Frequency (Hz)');

title('Frequency Drift During MRI Acquisition');

grid on;

%% --- Equivalent phase change ---

DeltaPhase_frequency = ...
    2*pi*Seq.tEcho*DeltaFrequency;

figure;

plot(t, DeltaPhase_frequency, '-s', ...
    'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Equivalent phase change (rad)');

title('Frequency-Related Phase Evolution');

grid on;

fprintf('\nFrequency drift analysis:\n');

fprintf('Number of frequency measurements: %d\n', numel(f));

fprintf('Total frequency change: %.4f Hz\n', ...
    DeltaFrequency(end));

fprintf('Equivalent phase change: %.4f rad\n', ...
    DeltaPhase_frequency(end));


%% --- DIAGNOSTIC: FREQUENCY-ADJUSTED PHASE ---

% Frequency drift contribution
% This is a diagnostic estimate, not yet a validated correction.

FrequencyPhaseCorrection = ...
    2*pi*Seq.tEcho*DeltaFrequency;

% Match frequency data to phase data
nCommon = min([ ...
    numel(Deltaphase_p), ...
    numel(FrequencyPhaseCorrection)]);

PhaseOriginal = Deltaphase_p(1:nCommon);

PhaseFrequencyAdjusted = ...
    PhaseOriginal - FrequencyPhaseCorrection(1:nCommon);

% Convert to temperature
TemperatureMRI_original = ...
    PhaseOriginal ./ ...
    (gamma_val*alpha_used*B0*TE);

TemperatureMRI_adjusted = ...
    PhaseFrequencyAdjusted ./ ...
    (gamma_val*alpha_used*B0*TE);

% Plot comparison

figure;

plot(Timedata_p(1:nCommon), ...
    TemperatureMRI_original(1:nCommon), ...
    '-o','LineWidth',1.5);

hold on;

plot(Timedata_p(1:nCommon), ...
    TemperatureMRI_adjusted, ...
    '-s','LineWidth',1.5);

plot(Timedata_p(1:nCommon), ...
    DeltaT_Osensa(1:nCommon), ...
    '-d','LineWidth',1.5);

xlabel('Time (s)');
ylabel('Temperature change (^oC)');

legend('Original MRI',...
       'Frequency-adjusted MRI',...
       'Osensa');

title('Effect of Frequency Drift Compensation');

grid on;


%% --- FREQUENCY VS MRI PHASE CORRELATION ---

% Convert everything to column vectors
t_all = Timedata(:);
f_all = FrequencyMean(:);
phase_all = Phasedata_unwrapped(:);

% Check original lengths
fprintf('\n--- ARRAY LENGTH CHECK ---\n');

fprintf('Timedata: %d\n',numel(t_all));
fprintf('FrequencyMean: %d\n',numel(f_all));
fprintf('Phasedata_unwrapped: %d\n',numel(phase_all));

% Match array lengths
nCommon = min([ ...
    numel(t_all), ...
    numel(f_all), ...
    numel(phase_all)]);

t = t_all(1:nCommon);
f = f_all(1:nCommon);
phase = phase_all(1:nCommon);

% Remove invalid measurements
valid = isfinite(t) & ...
        isfinite(f) & ...
        isfinite(phase);

% Apply the SAME logical index to all arrays
t = t(valid);
f = f(valid);
phase = phase(valid);

% Reference measurement
refIdx = 5;

if numel(t) < refIdx
    error('Not enough valid measurements for reference index.');
end

DeltaF = f - f(refIdx);
DeltaPhi = phase - phase(refIdx);

%% Correlation

R = corrcoef(DeltaF,DeltaPhi);

fprintf('\n--- FREQUENCY / PHASE CORRELATION ---\n');

fprintf('Correlation coefficient: %.5f\n',R(1,2));

%% Plot

figure;

subplot(2,1,1);

plot(t,DeltaF,'-o');

ylabel('\Delta Frequency (Hz)');
title('Frequency Drift');

grid on;

subplot(2,1,2);

plot(t,DeltaPhi,'-s');

xlabel('Time (s)');
ylabel('\Delta Phase (rad)');

title('MRI Phase Evolution');

grid on;

%% Linear relationship

pFreqPhase = polyfit(DeltaF,DeltaPhi,1);

fprintf('Phase/frequency slope: %.8f rad/Hz\n',...
    pFreqPhase(1));

fprintf('Theoretical slope at TE: %.8f rad/Hz\n',...
    2*pi*Seq.tEcho);


%% --- FREQUENCY / PHASE REGRESSION ---

figure;

scatter(DeltaF,DeltaPhi,60,'filled');

hold on;

% Measured regression
xFit = linspace(min(DeltaF),max(DeltaF),100);

yFit = polyval(pFreqPhase,xFit);

plot(xFit,yFit,'r-','LineWidth',2);

% Theoretical frequency-phase relationship
theoreticalPhase = 2*pi*Seq.tEcho*xFit;

plot(xFit,theoreticalPhase,'k--','LineWidth',2);

xlabel('Frequency change (Hz)');
ylabel('MRI phase change (rad)');

title('Frequency Drift vs MRI Phase');

legend('Measured data',...
       'Measured regression',...
       'Theoretical frequency relationship',...
       'Location','best');

grid on;

saveas(gcf,fullfile(saveDir,...
    ['FrequencyPhaseRegression_' baseName '.png']));

%% --- Frequency Predicted Phase Correction ---

% Frequency change relative to reference image
DeltaFrequency = DeltaFrequency(:);

% Theoretical phase caused by frequency drift
DeltaPhasePredicted = ...
    2*pi*Seq.tEcho*DeltaFrequency;

% Compare with measured MRI phase
MeasuredPhase = Deltaphase(:);

% Match array lengths safely
n = min(length(DeltaPhasePredicted), ...
        length(MeasuredPhase));

DeltaPhasePredicted = DeltaPhasePredicted(1:n);
MeasuredPhase = MeasuredPhase(1:n);
t = Timedata(1:n);

% Remove reference offset
DeltaPhasePredicted = ...
    DeltaPhasePredicted - DeltaPhasePredicted(1);

MeasuredPhase = ...
    MeasuredPhase - MeasuredPhase(1);

figure;

plot(t, MeasuredPhase, '-o', ...
    'LineWidth', 1.5);

hold on;

plot(t, DeltaPhasePredicted, '-s', ...
    'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Phase (rad)');

legend('Measured MRI phase', ...
       'Predicted frequency phase');

title('Measured Phase vs Frequency-Predicted Phase');

grid on;
%% --- Save data ---

save([fName '.mat'], ...
    'Timedata', ...
    'TemperatureData', ...
    'Phasedata', ...
    'Phasedata_unwrapped', ...
    'Deltaphase', ...
    'Acquisitiondata', ...
    'FrequencyCorrectionData');


%% --- Export processed data to CSV ---

% Force all variables into column vectors
Timedata_p = Timedata_p(:);
TemperatureData_p = TemperatureData_p(:);
Phasedata_p = Phasedata_p(:);

DeltaT_Osensa = DeltaT_Osensa(:);
DeltaT_MRI_pixel = DeltaT_MRI_pixel(:);
DeltaT_MRI_roi = DeltaT_MRI_roi(:);

Deltaphase_p = Deltaphase_p(:);


% Find common minimum number of elements
minLen = min([ ...
    numel(Timedata_p), ...
    numel(TemperatureData_p), ...
    numel(Phasedata_p), ...
    numel(DeltaT_Osensa), ...
    numel(DeltaT_MRI_pixel), ...
    numel(DeltaT_MRI_roi), ...
    numel(Deltaphase_p)]);

% Construct output matrix
data_out = [ ...
    Timedata_p(1:minLen), ...
    TemperatureData_p(1:minLen), ...
    Phasedata_p(1:minLen), ...
    DeltaT_Osensa(1:minLen), ...
    DeltaT_MRI_pixel(1:minLen), ...
    DeltaT_MRI_roi(1:minLen), ...
    Deltaphase_p(1:minLen)];

% Column headers
headers = { ...
    'Time_s', ...
    'Temperature_C', ...
    'Phase_unwrapped_rad', ...
    'DeltaT_Osensa_C', ...
    'DeltaT_MRI_pixel_C', ...
    'DeltaT_MRI_ROI_C', ...
    'DeltaPhase_rad'};

% Write CSV with headers
csvFile = [fName '_processed.csv'];

T = array2table(data_out, 'VariableNames', headers);

writetable(T, csvFile);

