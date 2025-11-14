%%GRE Sequence Phase Drive analysis (modified Exaple_Flash_2D code)
   
   %% Gradient Echo 2D (Flash 2D)
% This example acquires a 2D gradient echo image with Ernst angle excitation.

%%
% --- Set up save directory in Google Drive ---
baseFolder = 'G:\My Drive\MR_Thermometry\My runs';             % main directory on Google Drive, can change this to match the path in your Google Drive
todayFolder = datestr(now, 'yyyy-mm-dd');           % folder named with today's date
saveDir = fullfile(baseFolder, todayFolder);         % full path for today's run

if ~exist(saveDir, 'dir')
    mkdir(saveDir);                                 % create folder if it doesn't exist
end

timestamp = datestr(now, 'HHMMSS');

%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                      % number of loop averages 1...

%%PARAMETERS TO DEFINE
Seq.T1 = 100e-3;                                    % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 3e-3;                                   % echo time in seconds e.g. 4e-3
Seq.tRep = 8e-3;                                   % repetition time in seconds (default is Seq.tEcho*2)
resolution = 16;
thickness = 0.002;
pausetime = 2;
position = resolution/2;

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = resolution;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(2) = resolution;                      % number of pixels in phase direction
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.010;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.010;                % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = thickness;                   % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect; %Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

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
% Seq.AQSlice(1).alfa = 0.5*pi;                     % un-comment to exchange read and phase direction

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor

Seq.CorrectSliceRephase = 0;                        % Correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % Correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;

i = 0;
tStart = tic;
% figure1 = figure('Name', 'Amplitude at one position over time');
% figure2 = figure('Name', 'Phase at one position over time');
figure3 = figure('Name', 'Phase vs. Temperature');
figure4 = figure('Name', 'Phasedifference vs. Temperature');


Timedata = zeros(100,1);
Temperatures = zeros(100,1);
Phasedata = zeros(100,1);
Deltaphase = zeros(100,1);
Deltatemp = zeros(100,1);

while true
    i = i+1;
    repeat = strcat('v',num2str(i));
    osensa_port = 'COM3';       %change it
    ospy = py.importlib.import_module('osensaMatlab'); % change osensAssi
    transmitter = ospy.Transmitter(osensa_port, uint16(247));

    disp("Acquisition of Image " + num2str(i));
    
    time.(repeat) = toc(tStart);
    %% Here the actual sequence is run:
    [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

    %measure the current temperature with Osensa
    disp(' ');
    disp("##-- Measuring Temp --##")
    temp = round(transmitter.read_channel_temp("A"),3);
   
    disp("The measured temperature is " + num2str(temp));
    
    %Save image data for each repetition
    Acquisitiondata.(repeat) = SeqLoop;
    Timedata(i,1) = time.(repeat);
    Temperatures(i, 1) = temp;
    
    Imagephase = unwrap(angle(SeqLoop.data.Image(position, 1, position)));
    %Average the phase over an area x=20, y = 15:25
%     phase = [];
%     for y = (15:25)
%         phase(end+1) =(unwrap(angle(SeqLoop.data.Image(y, 1, 20)))); %first entry defines row position=y
%     end
%     Imagephase = mean(phase);
    
    if Imagephase < 0
        Imagephase = Imagephase + 2*pi;
    end
    
    Phasedata(i, 1) = Imagephase;
    
    Referencephase = Phasedata(5, 1);
    Referencetemp = Temperatures(5, 1);
    
    if (i <= 5)
        deltaphase = 0;
        deltatemp = 0;
    end
    
    if (i > 5)
        deltaphase = Imagephase - Referencephase;
        deltatemp = temp - Referencetemp;
    end
        
    Deltaphase(i, 1) = deltaphase;
    Deltatemp(i,1) = deltatemp;
    %plot the amplitude over time
%     figure(figure1)
%     axes1 = gca; hold on;
%     plot(gca, time.(repeat), abs(SeqLoop.data.Image(20, 1, 20)), 'o')
%     xlabel('Time in s')
%     ylabel('Image Amplitude')
    %plot the phase over time
%     figure(figure2)
%     axes2 = gca; hold on;
%     plot(gca, time.(repeat), unwrap(angle(SeqLoop.data.Image(20, 1, 20))), 'o')
%     xlabel('Time in s')
%     ylabel('Image Phase')
    %plot the phase vs temperature
    figure(figure3)
    axes3 = gca; hold on;
    plot(gca, temp, Imagephase, '-o')
    xlabel('Temperature in C')
    ylabel('Image Phase in rad')
    
    figure(figure4)
    axes4 = gca; hold on;
    plot(gca, deltatemp, deltaphase, '-o')
    xlabel('Temperaturedifference in C')
    ylabel('Phasedifference in rad')
    
    pause(pausetime)
    
    if (temp < 27)
        break
    end
end

%plot amplitude over time to see when the sample reaches steadz state (this
%is when we want to start with phase analysis, the images before will be ignored)

%plot image and phase of the last image acquired
figure5 = figure('Name', 'Delta Phase over time');
figure(figure5)
% subplot(2,2,1)
% Imagemangnitude = abs(reshape(SeqLoop.data.Image(:,1,:), [32,32]));
% imagesc(Imagemangnitude)
% axis equal
% title('Magnitude of the Image')
% 
% subplot(2,2,2)
% Imagephases = angle(reshape(SeqLoop.data.Image(:,1,:), [32,32]));
% imagesc(Imagephases)
% axis equal
% title('Image Phasemap')
% 
% subplot(2,2,3)
plot(Deltatemp(1:i,1), Deltaphase(1:i,1), '-o')
xlabel('Temperaturedifference in C')
ylabel('Delta Phase in rad')

% save all the data in a csv file
%save temperature, phasevalue at chosen location and time of each acquisition
Alldata = zeros(length(Timedata), 5);
Alldata(:, 1) = Timedata;
Alldata(:, 2) = Temperatures;
Alldata(:, 3) = Phasedata;
Alldata(:, 4) = Deltaphase;
Alldata(:, 5) = Deltatemp;

Table = array2table(Alldata);
Table.Properties.VariableNames(1:5)={'Time of measurement in s', 'Temperature in C', 'Absolute Phase at chosen location in rad', 'Delta Phase at location in rad','Delta Temperature in C'};
writetable(Table, fullfile(saveDir, ['Phasedata_resolution16_TE3ms_TR8ms_trial5' timestamp '.csv']))

saveas(figure5, fullfile(saveDir, ['DeltaphasevsT_resolution16_TE3ms_TR8ms_trial5' timestamp '.jpg']));

%% added by Jana for unwraped phase analysis
% hf = figure(10); clf(hf);
% hax = axes(hf);
% plot(hax, SeqLoop.data.time_of_tRep(:,2,15), (angle(SeqLoop.data.data(:,2,20))));
% title(hax, 'Phase of acquired signal');
% ylabel(hax, 'phase in rad');
% xlabel(hax, 'Time in s');
% grid(hax, 'on');
% % 
% hf = figure(11); clf(hf);
% hax = axes(hf);
% plot(hax, SeqLoop.data.time_of_tRep(:,2,15), unwrap(angle(SeqLoop.data.data(:,2,20))));
% title(hax, 'Unwrapped Phase of acquired signal');
% ylabel(hax, 'phase in rad');
% xlabel(hax, 'Time in s');
% grid(hax, 'on');

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
