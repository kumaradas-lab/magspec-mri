%% Spin Echo V.1       10.10.2013
%                      Toni, Driessle - Pure Devices GmbH
% Spin Echo 1 2 and 3 dimension
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

warning('PD:sequence_Spin_Echo_Multi:deprecated', ...
  ['This script is deprecated and might be removed in a future version ' ...
   'of OpenMatlab. Consider using the function "sequence_Spin_Echo" with ' ...
   'Seq.AQSlice(1).TurboFactor>1 instead.']);

if exist('HW', 'var') && isa(HW, 'PD.HW')
    ResetStructs;  % reset variables: Seq AQ TX Grad
else
    % mySave.DummySerial = 0;  % for Dummy measurements without an MRI device (0 or [] for real measurements or 1 for dummy measurements)
    LoadSystem;  % load system parameters (reset all variables HW Seq AQ TX Grad)
end
HW.FindFrequencyPause=2;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 20);  % Find magnet frequency

iDevice = 1;  % FIXME: Support multiple MMRT devices


% HW.Grad.TimeDelay=52e-6-15e-6;
HW.Grad(iDevice).tRamp=18e-6;
HW.Grad(iDevice).tEC=18e-6;

HW.Grad(iDevice).SliceTimeDelay=[0,0,0]*1e-6;
HW.Grad(iDevice).ReadTimeDelay=[0,0,0]*1e-6+0e-6;
HW.Grad(iDevice).PhaseTimeDelay=[0,0,0]*1e-6;

Seq.average=1;                                          % Number of averages    1...
Seq.averageBreak=0.2;                                     % Pause between two averages in seconds

% Seq.tRep=50e-3;                                        % Repetition time in seconds
Seq.tEcho=5e-3;                                        % Echo time in seconds

Seq.plotSeqTR=[1:3];                                      % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
Seq.plotSeq=1:3;                                        % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)

Seq.SteadyState_PreShots=2;
Seq.SteadyState_PostShots=2;

Seq.CLTime=10e-6;

% Readout at echo
% rotate all image coordinates around... rotate
Seq.AQSlice(1).plotkSpace=1;                               % Plot k-Space
Seq.AQSlice(1).plotImage=1;                                % Plot Image
Seq.AQSlice(1).plotPhase=1;                                % Plot phase of k-Space and Image
Seq.AQSlice(1).alfa=0.5*pi;                                % first rotation around x axis in RAD
Seq.AQSlice(1).phi=0.0*pi;                                 % second rotation around y axis in RAD
Seq.AQSlice(1).theta=0.5*pi;                               % third rotation around z axis in RAD
Seq.AQSlice(1).HzPerPixMin=500;                            % Bandwith per pixel in Hz
Seq.AQSlice(1).nRead=64;                                   % Number of Pixels in read,  nPhase(1)=1 if nRead>1
Seq.AQSlice(1).nPhase(1)=1;                                % Number of Pixels in CSI, nRead=1 if nPhase(1)>1
Seq.AQSlice(1).nPhase(2)=128;                              % Number of Pixels in phase
Seq.AQSlice(1).nPhase(3)=1;                                % Number of Pixels in phase3D
Seq.AQSlice(1).sizeRead=0.0064+0e12+(Seq.AQSlice(1).nRead==1)*1e12;                            % Image size in read in meter
Seq.AQSlice(1).sizePhase(1)=0.0064+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).sizePhase(2)=0.0064*2+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).sizePhase(3)=0.0064+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).thickness=0.002+0e12;                       % Image thickness in slice direction  used for 2D and 3D! in meter
% Seq.AQSlice(1).ReadOSUsedForImage=17;                      % Number of samples used in CSI
Seq.AQSlice(1).ReadOS=4;                                   % Oversampling read         (recommended >=2) 1...
Seq.AQSlice(1).PhaseOS(1)=1;                               % Oversampling phase        1...
Seq.AQSlice(1).PhaseOS(2)=1;                               % Oversampling phase CSI        1...
Seq.AQSlice(1).PhaseOS(3)=1;                               % Oversampling phase3D      1...
Seq.AQSlice(1).UsetRep=1+Seq.SteadyState_PreShots+(1:prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)).';
Seq.AQSlice(1).UseAQWindow=1;
Seq.AQSlice(1).ReadCoordinate=1;                           % direction of Read:  x = 1,  y = 2, z = 3
Seq.AQSlice(1).PhaseCoordinate(1)=1;                       % direction of Phase:   x = 1,  y = 2, z = 3
Seq.AQSlice(1).PhaseCoordinate(2)=2;                       % direction of Phase:   x = 1,  y = 2, z = 3
Seq.AQSlice(1).PhaseCoordinate(3)=3;                       % direction of Phase:   x = 1,  y = 2, z = 3
Seq.AQSlice(1).SliceCoordinate=3;                          % direction of Slice:   x = 1,  y = 2, z = 3

% Generate tRep times
Seq.tRep=Seq.tEcho*ones(1,1+prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)+Seq.SteadyState_PreShots+Seq.SteadyState_PostShots);

% 90 degrees slice excitation
Seq.Slice(1).Pulse.Function=@Pulse_Sinc_2;
Seq.Slice(1).Pulse.FlipAngle=90;
Seq.Slice(1).Thickness=Seq.AQSlice(1).thickness;
Seq.Slice(1).Pulse.Phase=0;
Seq.Slice(1).CenterOfPulse=Seq.tEcho/2;
Seq.Slice(1).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
Seq.Slice(1).GradTimeDelay=HW.Grad(iDevice).SliceTimeDelay;
Seq.Slice(1).UseAtRepetitionTime=1;

% 180 degrees (slice) inversion
Seq.Slice(2).Pulse.Function=@Pulse_RaisedCos;
Seq.Slice(2).Pulse.FlipAngle=180;
Seq.Slice(2).Thickness=Seq.AQSlice(1).thickness+1e12;  % use inversion slice?
Seq.Slice(2).CenterOfPulse=0;
% Seq.Slice(2).tRamp=100e-6;
Seq.Slice(2).Pulse.Phase=0;
Seq.Slice(2).GradDephaseSign=0;

% Seq.Slice(2).Pulse.PhaseIncrement=HW.Constant.GoldenAngle138;
Seq.Slice(2).Pulse.PhaseIncrement=180;
Seq.Slice(2).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
% Seq.Slice(2).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(2);
% Seq.Slice(2).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
Seq.Slice(2).UseAtRepetitionTime=2:length(Seq.tRep);
Seq.Slice(2).UseAtRepetitionTimeDephase=Seq.Slice(2).UseAtRepetitionTime-1;

% Readout at tEcho
Seq.Read(1).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
Seq.Read(1).CenterOfReadout=Seq.tEcho/2;
Seq.Read(1).Phase=0;
Seq.Read(1).PhaseIncrement=Seq.Slice(2).Pulse.PhaseIncrement;
Seq.Read(1).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
Seq.Read(1).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Read(1).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
% Seq.Read(1).Overdrive=0;
% Seq.Read(1).tRamp=HW.Grad.tRamp*1;
% Seq.Read(1).tEC=0e-3;
% generate readout timings
Seq=get_ReadParameter(Seq, HW);
% Seq.Read(1).CenterOfDephase=Seq.Read(1).CenterOfDephase- 0e-3;
Seq.Slice(2).CenterOfRephase=Seq.Read(1).CenterOfDephase;
Seq.Slice(2).GradRephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Slice(2).CenterOfDephase=Seq.Read(1).CenterOfRephase;
Seq.Slice(2).GradDephaseLength=Seq.Read(1).GradRephaseLength;


% phase encoding aligned with read dephase and rephase
Seq.Phase(1).sizePhase=Seq.AQSlice(1).sizePhase(1);
Seq.Phase(1).nPhase=Seq.AQSlice(1).nPhase(1);
Seq.Phase(1).PhaseOS=Seq.AQSlice(1).PhaseOS(1);
Seq.Phase(1).StepIncrement=1;
Seq.Phase(1).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(1).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(1).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(1).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(1).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(1);
Seq.Phase(1).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(1).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;

% phase3D encoding aligned with read dephase and rephase
Seq.Phase(2).sizePhase=Seq.AQSlice(1).sizePhase(2);
Seq.Phase(2).nPhase=Seq.AQSlice(1).nPhase(2);
Seq.Phase(2).PhaseOS=Seq.AQSlice(1).PhaseOS(2);
Seq.Phase(2).StepIncrement=Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1);
Seq.Phase(2).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(2).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(2).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(2).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(2).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(2);
Seq.Phase(2).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;

% phase3D encoding aligned with read dephase and rephase
Seq.Phase(3).sizePhase=Seq.AQSlice(1).sizePhase(3);
Seq.Phase(3).nPhase=Seq.AQSlice(1).nPhase(3);
Seq.Phase(3).PhaseOS=Seq.AQSlice(1).PhaseOS(3);
Seq.Phase(3).StepIncrement=Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)*Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2);
Seq.Phase(3).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(3).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(3).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(3).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(3).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(3);
Seq.Phase(3).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;

% generate slice parameters
Seq=get_SliceParameter(Seq, HW);
% generate phase parameters
Seq=get_PhaseParameter(Seq, HW);

% % Readout at FID
% Seq.AQSlice(2).plotkSpace=0;                               % Plot k-Space
% Seq.AQSlice(2).plotImage=1;                                % Plot Image
% Seq.AQSlice(2).plotPhase=2;                                % Plot phase of k-Space and Image
% Seq.AQSlice(2).alfa=1.5*pi;                                % first rotation around x axis in RAD
% Seq.AQSlice(2).phi=0.0*pi;                                 % second rotation around y axis in RAD
% Seq.AQSlice(2).theta=-0.0*pi;                               % third rotation around z axis in RAD
% Seq.AQSlice(2).nRead=96;                                   % Number of Pixels in read
% Seq.AQSlice(2).nPhase(1)=1;                                  % Number of Pixels in phase
% Seq.AQSlice(2).nPhase(2)=numel(Seq.tRep);                                  % Number of Pixels in phase
% Seq.AQSlice(2).nPhase(3)=1;                                  % Number of Pixels in phase
% Seq.AQSlice(2).nPhase3D=1;                                 % Number of Pixels in phase3D
% Seq.AQSlice(2).sizeRead=0.0064+1e12;                            % Image size in read in meter
% Seq.AQSlice(2).sizePhase(1)=0.0064+1e12;                           % Image size in phase in meter
% Seq.AQSlice(2).sizePhase(2)=0.0064+1e12;                           % Image size in phase in meter
% Seq.AQSlice(2).sizePhase(3)=0.0064+1e12;                           % Image size in phase in meter
% Seq.AQSlice(2).sizePhase3D=0.0064+1e12;                         % Image size in phase3D in meter
% Seq.AQSlice(2).thickness=0.0005+1e12;                       % Image thickness in slice direction  used for 2D and 3D! in meter
% Seq.AQSlice(2).ReadOS=2;                                   % Oversampling read         (recommended >=2) 1...
% Seq.AQSlice(2).PhaseOS(1)=1;                                  % Oversampling phase        1...
% Seq.AQSlice(2).PhaseOS(2)=1;                                  % Oversampling phase        1...
% Seq.AQSlice(2).PhaseOS(3)=1;                                  % Oversampling phase        1...
% Seq.AQSlice(2).UsetRep=1:prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS);
% Seq.AQSlice(2).UseAQWindow=1;
% Seq.AQSlice(2).ReadCoordinate=1;                                                   % direction of Read:  x = 1,  y = 2, z = 3
% Seq.AQSlice(2).PhaseCoordinate=2;                                                  % direction of Phase:   x = 1,  y = 2, z = 3
% Seq.AQSlice(2).Phase3DCoordinate=3;                                                % direction of Phase3D:   x = 1,  y = 2, z = 3
% Seq.AQSlice(2).SliceCoordinate=3;                                                  % direction of Slice:   x = 1,  y = 2, z = 3
%
% Seq.Read(2).HzPerPixelMin=1000;
% Seq.Read(2).CenterOfReadout=0.5/Seq.Read(2).HzPerPixelMin+Seq.Slice(1).CenterOfRephase+Seq.Slice(1).GradRephaseLength/2+HW.Grad.tEC;

% generate readout timings
Seq=get_ReadParameter(Seq, HW);

% AQ
% AQ=add_AQ(Seq.Read(1).AQ, Seq.Read(2).AQ);
AQ=Seq.Read(1).AQ;
% AQ.Frequency=AQ.Frequency(4);
% AQ.Phase=AQ.Phase(4);
% AQ.fSample=AQ.fSample(4);
% AQ.nSamples=AQ.nSamples(4);
% AQ.ResetPhases=[];
% TX
TX=add_TX(Seq.Slice(1).TX, Seq.Slice(2).TX);

% Gradients of image encoding
Grad=add_Grad(Grad, Seq.Slice(1).Grad);
Grad=add_Grad(Grad, Seq.Slice(1).GradRephase);
Grad=add_Grad(Grad, Seq.Slice(2).GradDephase);
Grad=add_Grad(Grad, Seq.Slice(2).Grad);
Grad=add_Grad(Grad, Seq.Slice(2).GradRephase);
Grad=add_Grad(Grad, Seq.Phase(1).GradDephase);
Grad=add_Grad(Grad, Seq.Phase(2).GradDephase);
Grad=add_Grad(Grad, Seq.Phase(3).GradDephase);
Grad=add_Grad(Grad, Seq.Read(1).GradDephase);
Grad=add_Grad(Grad, Seq.Read(1).Grad);
Grad=add_Grad(Grad, Seq.Read(1).GradRephase);
Grad=add_Grad(Grad, Seq.Phase(1).GradRephase);
Grad=add_Grad(Grad, Seq.Phase(2).GradRephase);
Grad=add_Grad(Grad, Seq.Phase(3).GradRephase);
Seq.Loops=1;
% Seq.average=1;
for Loop=1:Seq.Loops
  if Loop==1; Seq.Reinitialize=1;end
  if Seq.Loops>1; Seq.TimeToNextSequence=0.5; else Seq.TimeToNextSequence=[]; end
  if Loop==2
    Seq.Reinitialize=0;
    Seq.TimeFromLastSequence=Seq.TimeToNextSequence;
  end
  if Loop==Seq.Loops
    Seq.TimeFromLastSequence=Seq.TimeToNextSequence;
    Seq.TimeToNextSequence=[];
  end

  % Run sequence
  [~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);

  Channel = 1;  % FIXME: Add support for multiple acquisition channels?
  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

  % discard data from all but the used channels
  data = data(iAQ);

  % % get data of FID
  % dataFID = get_kSpaceAndImage(data, SeqOut.AQSlice(2));
  % dataFID = plot_kSpaceAndImage(dataFID, SeqOut.AQSlice(2));

  % get data of echo
  data = get_kSpaceAndImage(data, SeqOut.AQSlice(1));
  % data.Image(max(abs(data.Image(:)))/4>abs(data.Image(:))) = NaN;
  data = plot_kSpaceAndImage(data, SeqOut.AQSlice(1));
  data.SeqOut = SeqOut;
end
%%
if isfield(data,'ImageCsi')
%%
index=[4,6,4];
figure(2)
plot(data.ImageCsiFrequency(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsi(:,index(1),index(2),index(3))))%%
xlim([-20 20])
%%
index=[5,4,4];
figure(3)
plot(data.ImageCsiFrequencyZero(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsiRawZero(:,index(1),index(2),index(3))))
xlim([-20 20])
%%
end
