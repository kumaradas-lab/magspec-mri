%% Spin Echo 2D Difference
% Use this sequence to compare two measurements.
%
% First set up a Spin Echo 2D image with proper values.
% Choose if you want to measure the same sample or two different ones.
% Choose sequence parameters for second measurement if using the same sample.
% Hint: It is possible to save and load data.

LoadSystem;                                          % load system parameters (reset to default: HW Seq AQ TX Grad)

LoadDataFromFile= 1;                                 % [] for normal measurement or [pathname filename] and 1 to load data.
LoadDataFromFile= [];                                % uncomment for normal measurement

if LoadDataFromFile == 1
  [filename, pathname] = uigetfile('*.mat', 'File Selector');
  LoadDataFromFile=[pathname filename];
end
if isempty(LoadDataFromFile)

  % Use same sample?
  Seq.SameSample = 0;                                         % use same sample (1) or two samples (0)

  % Save data
  Seq.SaveSeqName = ['Test' '_' datestr(now, 'dd.mm.yyyy_HH.MM.SS')];     % save sequence (file name) or not save empty [];
  % Seq.SaveSeqName = [];                                                % path CurrentFolder\output\'Seq.SaveSeqName'

  %% Spin Echo 2D / Turbo Spin Echo 2D parameters

  Seq.Loops = 1;                                              % number of loop averages 1...

  Seq.tEcho = 5e-3;                                           % echo time in seconds e.g. 5e-3
  Seq.RepetitionTime = 100e-3;                                % repetition time in seconds

  % % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice(1).nRead = 32;                                  % number of Pixels in read, if nRead>1 nPhase(1)=1
  Seq.AQSlice(1).nPhase(2) = 32;                              % number of Pixels in phase(2)
  Seq.AQSlice(1).HzPerPixMin = 500;                           % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
  Seq.AQSlice(1).sizeRead = 0.01;                             % size in read direction in meter
  Seq.AQSlice(1).sizePhase(2) = 0.01;                         % size in phase(2) direction in meter
  Seq.AQSlice(1).thickness = 0.005;                           % slice thickness in meter
  Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;          % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
  Seq.AQSlice(1).inversionPulse = @Pulse_Rect;                % inversion pulse function

  % % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice(1).PhaseOS(2) = 2;                              % oversampling phase(2)  1...

  % % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice(1).TurboFactor = 1;                             % number of image k-lines per excitation
  % Seq.AQSlice(1).TurboBreak = Seq.tRep;                       % break between last echo and next excitation

  % % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

  % % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Seq.plotSeq = 1:3;                                          % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
  Seq.LoopPlot = 1;                                           % plot every loop
  Seq.AQSlice(1).plotkSpace = 1;                              % plot k-space
  Seq.AQSlice(1).plotImage = 1;                               % plot image
  Seq.AQSlice(1).plotPhase = 1;                               % plot phase of k-space or image


  %% run sequence_Spin_Echo with first sample
  if ~Seq.SameSample
    disp('Put first sample into the magnet and press any key')
    commandwindow();
    pause();
  end
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0); % find magnet frequency rough
  [HW, mySave] = Find_Frequency_FID(HW, mySave, 0); % find magnet frequency fine
  disp(['fLarmor = ' num2str(HW.fLarmor, '%9.0f') ' Hz']);
  [SeqDiff.Seq1, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

  %% Sequence parameter for second measurment if same sample.
  if Seq.SameSample
    % Seq.RepetitionTime = 50e-3;                                 % repetition time in seconds
    % Seq.tEcho = 25e-3;                                          % echo time in seconds e.g. 5e-3
  end

  %% run sequence_Spin_Echo with second sample
  if ~Seq.SameSample
    disp('Put second sample into the magnet and press any key')
    commandwindow();
    pause();
  end
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0); % find magnet frequency rough
  [HW, mySave] = Find_Frequency_FID(HW, mySave, 0); % find magnet frequency fine
  disp(['fLarmor = ' num2str(HW.fLarmor, '%9.0f') ' Hz']);
  [SeqDiff.Seq2, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

  fLarmorDiff = SeqDiff.Seq1.HW.fLarmor - SeqDiff.Seq2.HW.fLarmor;
  disp(['fLarmorDiff = ' num2str(fLarmorDiff, '%9.0f') ' Hz']);
else
  load(LoadDataFromFile);
end

%% image 1 and 2
SeqDiff.image1Complex = permute(squeeze(SeqDiff.Seq1.data.Image), [2,1]) / SeqDiff.Seq1.AQSlice(1).VoxelVolume * SeqDiff.Seq1.AQSlice(1).AreaCoil / SeqDiff.Seq1.AQSlice(1).AmplitudeUnitScale;
SeqDiff.image2Complex = permute(squeeze(SeqDiff.Seq2.data.Image), [2,1]) / SeqDiff.Seq2.AQSlice(1).VoxelVolume * SeqDiff.Seq2.AQSlice(1).AreaCoil / SeqDiff.Seq2.AQSlice(1).AmplitudeUnitScale;

hf = figure(5);
clf(hf);

% Image 1
hax(1) = subplot(2,2,1, 'Parent', hf);
imagesc(SeqDiff.Seq1.data.Ticks(1).Read, ...
        SeqDiff.Seq1.data.Ticks(2).Phase, ...
        abs(zeroFill_image(abs(SeqDiff.image1Complex), size(SeqDiff.image1Complex).*2, 1.4)), ...
        'Parent', hax(1));
set(hax(1), 'YDir', 'normal');
colormap(hax(1), gray);
title(hax(1), ['Image 1 amplitude in ' SeqDiff.Seq1.AQSlice(1).AmplitudeUnit ]);
xlabel(hax(1), 'Read in m');
ylabel(hax(1), 'Phase in m');
set(hax(1), 'DataAspectRatio', [1 1 1]);


hax(2) = subplot(2,2,2, 'Parent', hf);
imagesc(angle(SeqDiff.image1Complex), 'Parent', hax(2));
set(hax(2), 'CLim', [-pi,pi]);
set(hax(2), 'YDir', 'normal');
% colorbar
colormap(hax(2), gray);
set(hax(2), 'DataAspectRatio', [SeqDiff.Seq1.AQSlice(1).nRead/SeqDiff.Seq1.AQSlice(1).nPhase(2) ...
                                SeqDiff.Seq1.AQSlice(1).sizeRead/SeqDiff.Seq1.AQSlice(1).sizePhase(2) 1]);
title(hax(2), 'Image 1 phase in rad');


% Image 2
hax(3) = subplot(2,2,3, 'Parent', hf);
imagesc(SeqDiff.Seq2.data.Ticks(1).Read, ...
        SeqDiff.Seq2.data.Ticks(2).Phase, ...
        abs(zeroFill_image(abs(SeqDiff.image2Complex),size(SeqDiff.image2Complex).*2,1.4)), ...
        'Parent', hax(3));
set(hax(3), 'YDir', 'normal');
% colorbar
colormap(hax(3), gray);
title(hax(3), ['Image 2 amplitude in ' SeqDiff.Seq2.AQSlice(1).AmplitudeUnit ]);
xlabel(hax(3), 'Read in m');
ylabel(hax(3), 'Phase in m');
set(hax(3), 'DataAspectRatio', [1 1 1]);

hax(4) = subplot(2,2,4, 'Parent', hf);
imagesc(angle(SeqDiff.image2Complex), 'Parent', hax(4));
set(hax(4), 'CLim', [-pi,pi]);
set(hax(4), 'YDir', 'normal');
% colorbar
colormap(hax(4), gray)
set(hax(4), 'DataAspectRatio', [SeqDiff.Seq1.AQSlice(1).nRead/SeqDiff.Seq1.AQSlice(1).nPhase(2) ...
                                SeqDiff.Seq1.AQSlice(1).sizeRead/SeqDiff.Seq1.AQSlice(1).sizePhase(2) 1]);
title(hax(4), 'Image 2 phase in rad');

%% subtract both images absolute and angle

imageDiff = abs(SeqDiff.image1Complex) - abs(SeqDiff.image2Complex);
imageDiffAngle = mod(angle(SeqDiff.image1Complex ./ SeqDiff.image2Complex)+pi, pi*2) - pi; % angle -pi to pi; figure(1);plot((-10:0.1:10),mod((-10:0.1:10)+pi,pi*2)-pi);
% imageDiffAngle = unwrap2Dmiddle(imageDiffAngle); % use unwrap?  comment line if not needed.

hf = figure(6);
clf(hf);
hax(7) = subplot(1,2,1, 'Parent', hf);
imagesc(SeqDiff.Seq1.data.Ticks(1).Read, ...
        SeqDiff.Seq1.data.Ticks(2).Phase, ...
        real(zeroFill_image((imageDiff), size(imageDiff).*2, 1.4)), ...
        'Parent', hax(7));
set(hax(7), 'YDir', 'normal');
% colorbar
colormap(hax(7), gray);
title(hax(7), ['Difference (1-2) of image amplitude in ' SeqDiff.Seq1.AQSlice(1).AmplitudeUnit]);
xlabel(hax(7), 'Read in m');
ylabel(hax(7), 'Phase in m');
set(hax(7), 'DataAspectRatio', [1 1 1]);


hax(8) = subplot(1,2,2, 'Parent', hf);
imagesc(imageDiffAngle, 'Parent', hax(8));
set(hax(8), 'CLim', [-pi,pi]);
set(hax(8), 'YDir', 'normal');
% colorbar
colormap(hax(8), gray);
set(hax(8), 'DataAspectRatio', [SeqDiff.Seq1.AQSlice(1).nRead/SeqDiff.Seq1.AQSlice(1).nPhase(2) ...
                                SeqDiff.Seq1.AQSlice(1).sizeRead/SeqDiff.Seq1.AQSlice(1).sizePhase(2) 1]);
title(hax(8), 'Difference (1-2) of image phase in rad');


%% Image diff complex

imageDiffComplex = SeqDiff.image1Complex-SeqDiff.image2Complex;

hf = figure(7);
clf(hf);
hax(5) = subplot(1,2,1, 'Parent', hf);
imagesc(SeqDiff.Seq2.data.Ticks(1).Read, ...
        SeqDiff.Seq2.data.Ticks(2).Phase, ...
        abs(zeroFill_image(abs(imageDiffComplex),size(imageDiffComplex).*2,1.4)),'Parent',hax(5));
set(hax(5), 'YDir', 'normal');
% colorbar
colormap(hax(5), gray);
title(hax(5), ['Difference complex 1-2 amplitude in ' SeqDiff.Seq2.AQSlice(1).AmplitudeUnit]);
xlabel(hax(5), 'Read in m');
ylabel(hax(5), 'Phase in m');
set(hax(5), 'DataAspectRatio', [1 1 1]);

hax(6) = subplot(1,2,2, 'Parent', hf);
imagesc(angle(imageDiffComplex), 'Parent', hax(6));
set(hax(6), 'CLim', [-pi,pi]);
set(hax(6), 'YDir', 'normal');
% colorbar
colormap(hax(6), gray);
set(hax(6), 'DataAspectRatio', [SeqDiff.Seq1.AQSlice(1).nRead/SeqDiff.Seq1.AQSlice(1).nPhase(2) ...
                                SeqDiff.Seq1.AQSlice(1).sizeRead/SeqDiff.Seq1.AQSlice(1).sizePhase(2) 1]);

title(hax(6), 'Image complex 1-2 phase in rad');

%% Save data
if ~isempty(SeqDiff.Seq1.SaveSeqName) && isempty(LoadDataFromFile)
  if ~exist(fullfile(HW.RootPath, 'output'), 'dir')
    mkdir(fullfile(HW.RootPath, 'output'));
  end
  save(fullfile(HW.RootPath, 'output', [SeqDiff.Seq1.SaveSeqName '.mat']), 'SeqDiff');
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
