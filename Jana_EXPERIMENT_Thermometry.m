%% This code Measures T1 values using sequence_InversionSnapshot and Calculates the Sample-Temperature,
%% and compares it to the actual Osensa-Temperature
clear all; close all;
% --- Set up save directory in Google Drive ---
baseFolder = 'G:\My Drive\MR_Thermometry\My runs';             % main directory on Google Drive, can change this to match the path in your Google Drive
todayFolder = datestr(now, 'yyyy-mm-dd');           % folder named with today's date
saveDir = fullfile(baseFolder, todayFolder);         % full path for today's run

if ~exist(saveDir, 'dir')
    mkdir(saveDir);                                 % create folder if it doesn't exist
end

timestamp = datestr(now, 'HHMMSS');

% Preparations
LoadSystem;                                     % Load system parameters

%define concentration of the sample used
concentration = 3; %input concentration in mM

% Define parameters required for the measurement, e.g. oil T1 ~100e-3 sec
Seq.T1Estimated   = 100e-3;                     % estimated mean T1
% Seq.T1EstimatedMin = Seq.T1Estimated/3;         % minimum estimated T1,                         e.g. Seq.T1Estimated/3
% Seq.T1EstimatedMax = Seq.T1Estimated*3;         % maximum estimated T1,                         e.g. Seq.T1Estimated*3

%Global Variable to store - T1, Temp and Time respectively
X = zeros(60,14);
i=0;
total_time = 0;

figure1 = figure('Name', '1/T1 vs Time');
figure2 = figure('Name', '1/T1 vs Temperature');
figure3 = figure('Name', 'Osensa- vs MR-Temperature');
figure4 = figure('Name', 'Temperature-Error vs Time');


while true
    disp(' ');
    disp("##-- Measurement  " + num2str(i+1)+ " --##")
    disp("##-- Measuring T1 --##")

    %Start Timer
    tStart = tic ; 

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
    %t1 = t1*0.001; % convert t1 from ms to s
    
    %#####---- Adding Temperature code from here ----#
    % Import Osensapy Python Module
    osensa_port = 'COM3';   %change
    ospy = py.importlib.import_module('osensaMatlab');%change
    transmitter = ospy.Transmitter(osensa_port, uint16(247));
    disp(' ');
    disp("##-- Measuring Osensa-Temp --##")
    osensatemp = round(transmitter.read_channel_temp("A"),3);
    reading_time = round(toc(tStart),3);
   
    disp("The Osensa temperature is " + num2str(osensatemp));
    disp("The time for this reading is " +num2str(reading_time));
    
    total_time=total_time+reading_time;
    disp("Time elapsed: " + num2str(total_time));
    
    disp("##-- Doing calculation --##")

    r1 = 1/T1.tau; %/concentration; % make relaxivity concentration independant
    disp(['T1.tau = ', num2str(T1.tau), ' s']);


    disp(['Relaxivity r1 = ' num2str(r1) ' Hz']);
    disp("##-- Getting MRI-Temp --##")

    syms x;
    % Multiple equations determined by measurements of differnt
    % concentrations all leading to the determination of Temperature
    % based on Khalids experiment: https://drive.google.com/drive/u/1/folders/1c_4MUUdrv7qmMD6vtJd2EbW4qG1jLtsF
    eqn1 = r1 == -0.0013048*x^2 + 0.027533*x + 5.2391; %equation based on 1mM experiments
    eqn2 = r1 == -0.00111915*x^2 + 0.021059*x + 5.17005; %equation based on 2mM experiments
    eqn3 = r1 == -0.00158477*x^2 + 0.052*x + 4.631; %equation based on 3mM experiments
    eqn4 = r1 == -0.0018003*x^2 + 0.0675375*x + 4.4573; %equation based on 4mM experiments
    eqn5 = r1 == -0.00130026*x^2 + 0.031408*x + 5.07944; %equation based on 5mM experiments
    eqn6 = r1 == -0.00138682*x^2 + 0.03734667*x + 5.03225; %equation based on 6mM experiments

    S1 = solve(eqn1, x);
    S2 = solve(eqn2, x);
    S3 = solve(eqn3, x);
    S4 = solve(eqn4, x);
    S5 = solve(eqn5, x);
    S6 = solve(eqn6, x);

    Temperature1 = double(S1);
    Temperature2 = double(S2);
    Temperature3 = double(S3);
    Temperature4 = double(S4);
    Temperature5 = double(S5);
    Temperature6 = double(S6); % to show S numerically
    %Take absolute of temperature for Temperature > 15 deg. C!!
    Temperature1 = abs(Temperature1(2));
    Temperature2 = abs(Temperature2(2));
    Temperature3 = min(abs(Temperature3(2)));
    Temperature4 = abs(Temperature4(2));
    Temperature5 = abs(Temperature5(2));
    Temperature6 = abs(Temperature6(2));

    disp(['The Temperature is T1 = ' ,num2str(Temperature1), ' deg. C']);
    disp(['The Temperature is T2 = ' ,num2str(Temperature2), ' deg. C']);
    disp(['The Temperature is T3 = ' ,num2str(Temperature3), ' deg. C']);
    disp(['The Temperature is T4 = ' ,num2str(Temperature4), ' deg. C']);
    disp(['The Temperature is T5 = ' ,num2str(Temperature5), ' deg. C']);
    disp(['The Temperature is T6 = ' ,num2str(Temperature6), ' deg. C']);

    Temperature = (Temperature1+Temperature2+Temperature3+Temperature4+Temperature5+Temperature6)/6;

    disp(['The average Temperature is T = ' ,num2str(Temperature), ' deg. C']);


    
    error = osensatemp - Temperature;
    
    %Save T1,temp and time values
    i=i+1;
    X(i,1)= T1.tau;
    X(i,2)= round(1/(T1.tau),4);
    X(i,3)= osensatemp;
    X(i,4) = total_time;
    X(i,5) = reading_time;
    X(i,6) = Temperature1;
    X(i,7) = Temperature2;
    X(i,8) = Temperature3;
    X(i,9) = Temperature4;
    X(i,10) = Temperature5;
    X(i,11) = Temperature6;
    X(i,12) = Temperature;
    X(i,13) = concentration;
    X(i,14) = error;
    
    figure(figure1)
    axes1 = gca; hold on;
    plot(gca, X(1:i,4), X(1:i,2), 'go');
    plot(gca, X(1:i,4), X(1:i,2), 'b-');
    xlabel('Time in s');
    ylabel('1/T1 in Hz');
    %hold (axes1, 'all');
    
    figure(figure2)
    axes2 = gca; hold on;
    plot(gca, X(1:i,3), X(1:i,2), 'mo');
    plot(gca, X(1:i,3), X(1:i,2), 'r-');
    xlabel('Osensa-Temperature in deg C');
    ylabel('1/T1 in Hz');
    
    figure(figure3)
    axes3 = gca; hold on;
    plot(gca, X(1:i,3), X(1:i,12), 'mo');
    plot(gca, X(1:i,3), X(1:i,12), 'r-');
    xlabel('Osensa-Temperature in deg C');
    ylabel('MRI-Temperature in deg C');
    
    figure(figure4)
    axes4 = gca; hold on;
    plot(gca, X(1:i,4), X(1:i,14), 'mo');
    plot(gca, X(1:i,4), X(1:i,14), 'r-');
    xlabel('Time in s');
    ylabel('Temperature-Error in deg C');
    
    if (osensatemp <= 26) %(i==1)%SET END OF CODE HERE
       break
    end
    
end

%Saving values to CSV
T = array2table(X);
T.Properties.VariableNames(1:14) = {'T1 (s)', '1/T1 (Hz)','Temp(C)','Time Elapsed(s)', 'Reading Time(s)', 'Temperature1 (C)', 'Temperature2 (C)', 'Temperature3 (C)', 'Temperature4 (C)', 'Temperature5 (C)', 'Temperature6 (C)', 'Average Temperature (C)', 'concentration (mM=mmole/L)','error (C)'};
writetable(T,fullfile(saveDir, ['MR-Thermometry_concentr'  num2str(concentration) '3.csv']))

%Save your figure as JPG
saveas(figure1, fullfile(saveDir, ['MR-Thermometry_R1vstime_concentr'  num2str(concentration) '3.jpg']))
saveas(figure2, fullfile(saveDir, ['MR-Thermometry_R1vsTemp_concentr'  num2str(concentration) '3.jpg']))
saveas(figure3, fullfile(saveDir, ['MR-Thermometry_OsensaTvsMRIT_concentr' num2str(concentration) '3.jpg']))
saveas(figure4, fullfile(saveDir, ['MR-Thermometry_error_concentr'  num2str(concentration) '3.jpg']))

%Close osens port 
transmitter.close()