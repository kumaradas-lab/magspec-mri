clear; clc; close all;

% Repetition times to compare
TRs = [200 300 400]; %[200 400 600 800];

% Colors for each TR
colors = lines(length(TRs));

figure;

%% --------- Delta T Osensa ----------
subplot(3,1,1);
hold on; grid on;

for k = 1:length(TRs)

    filename = sprintf('*TR%dms*_processed.csv', TRs(k));

    f = dir(filename);

    if isempty(f)
        warning('No file found for TR = %d ms', TRs(k));
        continue;
    end

    data = readmatrix(f(1).name);

    time = data(:,1);
    deltaOsensa = data(:,3);

    plot(time, deltaOsensa, '-o', ...
        'Color', colors(k,:), ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('TR = %d ms',TRs(k)));
end

xlabel('Time (s)');
ylabel('\DeltaT Osensa (^\circC)');
title('Osensa Temperature');
legend('Location','best');


%% --------- Delta T MRI Pixel ----------
subplot(3,1,2);
hold on; grid on;

for k = 1:length(TRs)

    filename = sprintf('*TR%dms*_processed.csv', TRs(k));

    f = dir(filename);

    if isempty(f)
        continue;
    end

    data = readmatrix(f(1).name);

    time = data(:,1);
    deltaPixel = data(:,4);

    plot(time, deltaPixel, '-o', ...
        'Color', colors(k,:), ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('TR = %d ms',TRs(k)));
end

xlabel('Time (s)');
ylabel('\DeltaT MRI Pixel (^\circC)');
title('MRI Pixel');


%% --------- Delta T MRI ROI ----------
subplot(3,1,3);
hold on; grid on;

for k = 1:length(TRs)

    filename = sprintf('*TR%dms*_processed.csv', TRs(k));

    f = dir(filename);

    if isempty(f)
        continue;
    end

    data = readmatrix(f(1).name);

    time = data(:,1);
    deltaROI = data(:,5);

    plot(time, deltaROI, '-o', ...
        'Color', colors(k,:), ...
        'LineWidth',1.8, ...
        'DisplayName',sprintf('TR = %d ms',TRs(k)));
end

xlabel('Time (s)');
ylabel('\DeltaT MRI ROI (^\circC)');
title('MRI ROI');

sgtitle('Comparison of Different Repetition Times');


% clear; clc; close all;
% 
% TRs = [200 300 400]; %[200 400 600 800];
% colors = lines(length(TRs));
% 
% figure;
% hold on;
% grid on;
% 
% for k = 1:length(TRs)
% 
%     filename = sprintf('*TR%dms*_processed.csv', TRs(k));
%     f = dir(filename);
% 
%     if isempty(f)
%         warning('No file found for TR = %d ms', TRs(k));
%         continue;
%     end
% 
%     data = readmatrix(f(1).name);
% 
%     t = data(:,1);
% 
%     % Plot all three measurements
%     plot(t, data(:,3), '-', ...
%         'Color', colors(k,:), ...
%         'LineWidth',2, ...
%         'DisplayName', sprintf('Osensa  TR=%d ms',TRs(k)));
% 
%     plot(t, data(:,4), '--', ...
%         'Color', colors(k,:), ...
%         'LineWidth',2, ...
%         'DisplayName', sprintf('MRI Pixel  TR=%d ms',TRs(k)));
% 
%     plot(t, data(:,5), ':', ...
%         'Color', colors(k,:), ...
%         'LineWidth',2.5, ...
%         'DisplayName', sprintf('MRI ROI  TR=%d ms',TRs(k)));
% 
% end
% 
% xlabel('Time (s)');
% ylabel('\DeltaT (^oC)');
% title('Temperature Change for Different Repetition Times');
% 
% legend('Location','eastoutside');
