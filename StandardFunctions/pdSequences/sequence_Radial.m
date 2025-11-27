%% Radial V.2           12.03.2014
%                      Toni, Driessle - Pure Devices GmbH
% Radial 1 2 and 3 dimension

function [SeqLoop, mySave] = sequence_Radial(HW, Seq, AQ, TX, Grad, mySave)
% LoadSystem
% mySave.DummySerial=0;                                   % for Dummy measurements without an MRI device (0 or [] for real measurements or 1 for dummy measurements)
% Seq=[];
% Seq.Loops=1;
% Seq.LoopsBreak=[];
% Seq.LoopsBreakExactly=0;
% Seq.LoopPlot=0;
%
% Seq.plotSeqTR=[1:3];                                      % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeq=1:3;                                        % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
%
% Seq.average=1;                                          % Number of averages    1...
% Seq.averageBreak=1;                                     % Pause between two averages in seconds
%
%
% Seq.T1=100e-3;                                             % T1 of probe excitation is acos(exp(-Seq.tRep/Seq.T1))/pi*180
% Seq.tEcho=7.5e-3;                                        % Echo time in seconds
%
% % Seq.tRep=10e-3;                                        % Repetition time in seconds
% Seq.tRep=Seq.tEcho*2 ;                               % Repetition time in seconds
%
% % get to steady state
% Seq.SteadyState_PreShots=10;
% Seq.SteadyState_PostShots=10;
% Seq.CorrectPhase=1;
% Seq.CorrectPlotFrequency=0;
%
% % Readout at echo
% % rotate all image coordinates around... rotate
% Seq.AQSlice(1).plotkSpace=1;                               % Plot k-Space
% Seq.AQSlice(1).plotImage=1;                                % Plot Image
% Seq.AQSlice(1).plotPhase=1;                                % Plot phase of k-Space and Image
% Seq.AQSlice(1).alfa=0.5*pi;                                % first rotation around x axis in RAD
% Seq.AQSlice(1).phi=0.0*pi;                                 % second rotation around y axis in RAD
% Seq.AQSlice(1).theta=0.0*pi;                               % third rotation around z axis in RAD
% Seq.AQSlice(1).HzPerPixMin=200;                            % Bandwith per pixel in Hz
% Seq.AQSlice(1).nRead=32;                                   % Number of Pixels in read,  nPhase(1)=1 if nRead>1
% Seq.AQSlice(1).nPhase(1)=1;                                % Number of Pixels in CSI, nRead=1 if nPhase(1)>1
% Seq.AQSlice(1).nPhase(2)=32;                               % Number of Pixels in phase
% Seq.AQSlice(1).nPhase(3)=1;                                % Number of Pixels in phase3D
% Seq.AQSlice(1).sizeRead=0.0064+(Seq.AQSlice(1).nRead==1)*1e12;                            % Image size in read in meter
% Seq.AQSlice(1).sizePhase(1)=0.0064+0e12;                           % Image size in phase in meter
% Seq.AQSlice(1).sizePhase(2)=0.0064+0e12;                           % Image size in phase in meter
% Seq.AQSlice(1).sizePhase(3)=0.0064+0e12;                           % Image size in phase in meter
% Seq.AQSlice(1).thickness=0.001+0e12;                       % Image thickness in slice direction  used for 2D and 3D! in meter
% Seq.AQSlice(1).ReadOSUsedForImage=17;                     % Number of samples used in CSI
% Seq.AQSlice(1).ReadOS=16*Seq.AQSlice(1).nPhase(1);                                   % Oversampling read         (recommended >=2) 1...
% Seq.AQSlice(1).PhaseOS(1)=1;                                  % Oversampling phase        1...
% Seq.AQSlice(1).PhaseOS(2)=32;                                  % Oversampling phase CSI        1...
% Seq.AQSlice(1).PhaseOS(3)=1;                                % Oversampling phase3D      1...
% Seq.AQSlice(1).UsetRep=Seq.SteadyState_PreShots+(1:prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS));
% Seq.AQSlice(1).UseAQWindow=1+(Seq.CorrectPhase>0);
% Seq.AQSlice(1).ReadCoordinate=1;                                                   % direction of Read:  x = 1,  y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(1)=1;                                                  % direction of Phase:   x = 1,  y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(2)=2;                                                  % direction of Phase:   x = 1,  y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(3)=3;                                                  % direction of Phase:   x = 1,  y = 2, z = 3
% Seq.AQSlice(1).SliceCoordinate=3;                                                  % direction of Slice:   x = 1,  y = 2, z = 3
% Seq.AQSlice(1).plotImageHandle=111+Loop;
%
% HW.Grad.SliceTimeDelay=[0,0,0]*1e-6;
% HW.Grad.ReadTimeDelay=[0,0,0]*1e-6;
% HW.Grad.PhaseTimeDelay=[0,0,0]*1e-6;
if ~isfield(Seq, 'AQSlice')
  Seq.AQSlice = struct();
end
if isemptyfield(Seq.AQSlice, 'HzPerPixMin')
  Seq.AQSlice(1).HzPerPixMin = 500;
end
if isemptyfield(Seq, 'Loops')
  Seq.Loops = 1;
end
if ~isfield(Seq, 'LoopsBreak')
  Seq.LoopsBreak = [];
end
if isemptyfield(Seq, 'LoopsBreakExactly')
  Seq.LoopsBreakExactly = 0;
end
if isemptyfield(Seq, 'LoopPlot')
  Seq.LoopPlot = 1;
end
if isemptyfield(Seq, 'tRep')
  Seq.tRep = 1/Seq.AQSlice(1).HzPerPixMin*2.5 + 1e-3;
end
if isemptyfield(Seq, 'average')
  Seq.average = 1;
end
if isemptyfield(Seq, 'T1')
  Seq.T1 = 100e-3;
end
if isemptyfield(Seq, 'SteadyState_PreShots')
  Seq.SteadyState_PreShots = 10;
end
if isemptyfield(Seq,'SteadyState_PostShots')
  Seq.SteadyState_PostShots = Seq.SteadyState_PreShots;
end
if isemptyfield(Seq, 'CorrectPhase')
  Seq.CorrectPhase = 0;
end
if isemptyfield(Seq, 'CorrectPlotFrequency')
  Seq.CorrectPlotFrequency = Seq.CorrectPhase;
end
if isemptyfield(Seq, 'StartWithKSpaceCenterNonlinear')
  Seq.StartWithKSpaceCenterNonlinear = 0;
end
if isemptyfield(Seq, 'StartWithGradientBeforeExcitationPulse')
  Seq.StartWithGradientBeforeExcitationPulse = 0;
end
if ~isfield(Seq,'AQSlice')
  Seq.AQSlice = struct();
end
if isemptyfield(Seq.AQSlice, 'ReadRadial')
  Seq.AQSlice(1).ReadRadial = 1;
end
if isemptyfield(Seq.AQSlice, 'FixedReadGrad')
  Seq.AQSlice(1).FixedReadGrad = 1;
end
if isemptyfield(Seq.AQSlice, 'ReadCoordinate')
  Seq.AQSlice(1).ReadCoordinate = 1;
end
if isemptyfield(Seq.AQSlice, 'PhaseCoordinate')
  Seq.AQSlice(1).PhaseCoordinate = [1 2 3];
end
if isemptyfield(Seq.AQSlice, 'SliceCoordinate')
  Seq.AQSlice(1).SliceCoordinate = 3;
end
if isemptyfield(Seq.AQSlice, 'plotImageHandle')
  Seq.AQSlice(1).plotImageHandle = 110;
end
if isemptyfield(Seq.AQSlice, 'UseAQWindow')
  Seq.AQSlice(1).UseAQWindow = 1+(Seq.CorrectPhase>0);
end
if isemptyfield(Seq.AQSlice, 'PhaseOS')
  Seq.AQSlice(1).PhaseOS = [1 1 1];
end
if isemptyfield(Seq.AQSlice, 'sizePhase')
  Seq.AQSlice(1).sizePhase = HW.Grad.ImageVol(2:2:6) - HW.Grad.ImageVol(1:2:6);
end
if isemptyfield(Seq.AQSlice, 'sizeRead')
  Seq.AQSlice(1).sizeRead = HW.Grad.ImageVol(2) - HW.Grad.ImageVol(1);
end
if isemptyfield(Seq.AQSlice, 'thickness')
  Seq.AQSlice(1).thickness = 1e12;
end
if isemptyfield(Seq.AQSlice, 'plotPhase')
  Seq.AQSlice(1).plotPhase = 1;
end
if isemptyfield(Seq.AQSlice, 'plotkSpace')
  Seq.AQSlice(1).plotkSpace = 1;
end
if isemptyfield(Seq.AQSlice, 'plotImage')
  Seq.AQSlice(1).plotImage = 1;
end
if isemptyfield(Seq.AQSlice, 'alfa')
  Seq.AQSlice(1).alfa = 0;
end
if isemptyfield(Seq.AQSlice, 'phi')
  Seq.AQSlice(1).phi = 0;
end
if isemptyfield(Seq.AQSlice, 'theta')
  Seq.AQSlice(1).theta = 0;
end
if isemptyfield(Seq.AQSlice, 'nRead')
  Seq.AQSlice(1).nRead = 32;
end
if isemptyfield(Seq.AQSlice, 'nPhase')
  Seq.AQSlice(1).nPhase = [1 32 1];
end
if Seq.AQSlice(1).nPhase(1) == 0
  Seq.AQSlice(1).nPhase(1) = 1;
end
if isscalar(Seq.AQSlice(1).nPhase)
  Seq.AQSlice(1).nPhase(2) = 1;
end
if Seq.AQSlice(1).nPhase(2) == 0
  Seq.AQSlice(1).nPhase(2) = 1;
end
if numel(Seq.AQSlice(1).nPhase) == 2
  Seq.AQSlice(1).nPhase(3) = 1;
end
if Seq.AQSlice(1).nPhase(3) == 0
  Seq.AQSlice(1).nPhase(3) = 1;
end
if Seq.AQSlice(1).PhaseOS(1) == 0
  Seq.AQSlice(1).PhaseOS(1) = 1;
end
if isscalar(Seq.AQSlice(1).PhaseOS)
  Seq.AQSlice(1).PhaseOS(2) = 1;
end
if Seq.AQSlice(1).PhaseOS(2) == 0
  Seq.AQSlice(1).PhaseOS(2) = 1;
end
if numel(Seq.AQSlice(1).PhaseOS) == 2
  Seq.AQSlice(1).PhaseOS(3) = 1;
end
if Seq.AQSlice(1).PhaseOS(3) == 0
  Seq.AQSlice(1).PhaseOS(3) = 1;
end
if Seq.AQSlice(1).sizePhase(1) == 0
  Seq.AQSlice(1).sizePhase(1) = HW.Grad.ImageVol(2) - HW.Grad.ImageVol(1);
end
if isscalar(Seq.AQSlice(1).sizePhase)
  Seq.AQSlice(1).sizePhase(2) = HW.Grad.ImageVol(4) - HW.Grad.ImageVol(3);
end
if Seq.AQSlice(1).sizePhase(2) == 0
  Seq.AQSlice(1).sizePhase(2) = HW.Grad.ImageVol(4) - HW.Grad.ImageVol(3);
end
if numel(Seq.AQSlice(1).sizePhase) == 2
  Seq.AQSlice(1).sizePhase(3) = HW.Grad.ImageVol(6) - HW.Grad.ImageVol(5);
end
if Seq.AQSlice(1).sizePhase(3) == 0
  Seq.AQSlice(1).sizePhase(3) = HW.Grad.ImageVol(6) - HW.Grad.ImageVol(5);
end
if isemptyfield(Seq, 'plotSeqEnd')
  Seq.plotSeqEnd = [];
end
if isemptyfield(Seq.AQSlice, 'ReadOS')
  Seq.AQSlice(1).ReadOS = max(2, ceil(16000/(Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead)/2)*2);
end
if isemptyfield(Seq.AQSlice, 'ReadOSUsedForImage')
  Seq.AQSlice(1).ReadOSUsedForImage = min(Seq.AQSlice(1).ReadOS, 17);
end
if isemptyfield(Seq.AQSlice, 'excitationPulse')
  Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_1;
end


for Loop = 1:Seq.Loops

  if Loop==1
    init.AQ=AQ;
    init.TX=TX;
    init.Seq=Seq;
    init.Grad=Grad;
  else
    AQ=init.AQ;
    TX=init.TX;
    Grad=init.Grad;
    Seq=init.Seq;
  end

  if Loop==1 || Seq.LoopsBreakExactly==0
    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10);  % Find magnet frequency
  end
  % HW.Grad.TimeDelay=52e-6-15e-6;
  % HW.Grad.tRamp=18e-6;
  % HW.Grad.tEC=18e-6;


  if Loop==1
    % Seq.plotSeqTR= 1:3;  % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
    % Seq.plotSeq = 1:3;  % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
  else
    Seq.plotSeqTR = [];  % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
    Seq.plotSeq = [];  % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
  end

  % generate grid with rotations
  [Seq.kgrid.all1, Seq.kgrid.all2, Seq.kgrid.all3] = ndgrid( ...
        linspace(-0.5, 0.5-1/(Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)), Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1))*(Seq.AQSlice(1).nPhase(1)-1)/Seq.AQSlice(1).sizePhase(1),...
        linspace(-0.5, 0.5-1/(Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2)), Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2))*(Seq.AQSlice(1).nPhase(2)-1)/Seq.AQSlice(1).sizePhase(2),...
        linspace(-0.5, 0.5-1/(Seq.AQSlice(1).nPhase(3)*Seq.AQSlice(1).PhaseOS(3)), Seq.AQSlice(1).nPhase(3)*Seq.AQSlice(1).PhaseOS(3))*(Seq.AQSlice(1).nPhase(3)-1)/Seq.AQSlice(1).sizePhase(3));
  Seq.kgrid.selectWall = true(size(Seq.kgrid.all1));
  if Seq.StartWithGradientBeforeExcitationPulse == 0
    Seq.kgrid.selectWall(2:end-1,2:end-1,2:end-1)=false;
    if size(Seq.kgrid.all1,1) == 1
      Seq.kgrid.selectWall(1,2:end-1,2:end-1) = false;
    end
    if size(Seq.kgrid.all1,2) == 1
      Seq.kgrid.selectWall(2:end-1,1,2:end-1) = false;
    end
    if size(Seq.kgrid.all1,3) == 1
      Seq.kgrid.selectWall(2:end-1,2:end-1,1) = false;
    end
  end
  Seq.kgrid.wall1 = Seq.kgrid.all1(Seq.kgrid.selectWall);
  Seq.kgrid.wall2 = Seq.kgrid.all2(Seq.kgrid.selectWall);
  Seq.kgrid.wall3 = Seq.kgrid.all3(Seq.kgrid.selectWall);
  % convert rotations to points on the unit sphere
  [Seq.kgrid.wallphi,Seq.kgrid.walltheta,Seq.kgrid.wallr]=cart2sph(Seq.kgrid.wall1(:),Seq.kgrid.wall2(:),Seq.kgrid.wall3(:));
  if Seq.StartWithGradientBeforeExcitationPulse == 0
    [~,Seq.kgrid.Index]=sortrows([Seq.kgrid.wallphi,Seq.kgrid.walltheta],[1,2]);
    Seq.kgrid.wall1(:)=Seq.kgrid.wall1(Seq.kgrid(1).Index);
    Seq.kgrid.wall2(:)=Seq.kgrid.wall2(Seq.kgrid(1).Index);
    Seq.kgrid.wall3(:)=Seq.kgrid.wall3(Seq.kgrid(1).Index);
    Seq.kgrid.wallphi=Seq.kgrid.wallphi(Seq.kgrid.Index);
    Seq.kgrid.walltheta=Seq.kgrid.walltheta(Seq.kgrid.Index);
    Seq.kgrid.wallr=Seq.kgrid.wallr(Seq.kgrid.Index);
  else
    Seq.kgrid.Index=1:numel(Seq.kgrid.wallphi);
  end
  % Seq.kgrid.VectorOfNormalisedWavenumber=linspace(0,1,Seq.AQSlice(1).nRead*Seq.AQSlice(1).ReadOS).';
  Seq.Read(1).alfa=linspace(0,2*pi*(1-1/numel(Seq.kgrid.wallr)),numel(Seq.kgrid.wallr))*0+Seq.AQSlice(1).alfa;
  Seq.Read(1).phi=Seq.kgrid.wallphi.'*1+Seq.AQSlice(1).phi;
  Seq.Read(1).theta=Seq.kgrid.walltheta.'*1+Seq.AQSlice(1).theta;
  if Seq.StartWithGradientBeforeExcitationPulse == 0
    Seq.Read(1).Wavenumber=Seq.kgrid.wallr.'*(~Seq.AQSlice(1).FixedReadGrad)+max(Seq.kgrid.wallr(:))*(Seq.AQSlice(1).FixedReadGrad)+1e-6;
  else
    Seq.Read(1).Wavenumber=Seq.kgrid.wallr.';
  end
  % Seq.Read(1).Resolution=1/(   1/(Seq.AQSlice(1).sizePhase(1)/(Seq.AQSlice(1).nPhase(1)-1)).^2 ...
  %                             +1/(Seq.AQSlice(1).sizePhase(2)/(Seq.AQSlice(1).nPhase(2)-1)).^2 ...
  %                             +1/(Seq.AQSlice(1).sizePhase(3)/(Seq.AQSlice(1).nPhase(3)-1)).^2).^0.5+zeros(size(Seq.Read(1).Wavenumber));
  Seq.Read(1).Resolution=1./Seq.Read(1).Wavenumber;

  % Generate tRep times
  Seq.tRep=Seq.tRep*ones(1,Seq.SteadyState_PostShots+Seq.SteadyState_PostShots+numel(Seq.kgrid.wallr));
  if isemptyfield(Seq.AQSlice, 'UsetRep')
    Seq.AQSlice(1).UsetRep = Seq.SteadyState_PreShots + (1:numel(Seq.kgrid.wallr)).';
  end

  % Ernst angle slice excitation
  Seq.Slice(1).Pulse.Function=Seq.AQSlice(1).excitationPulse;
  Seq.Slice(1).Pulse.FlipAngle=acos(exp(-Seq.tRep(1)/Seq.T1))/pi*180;
  Seq.Slice(1).Thickness=Seq.AQSlice(1).thickness;
  Seq.Slice(1).Pulse.Phase=0;
  Seq.Slice(1).Pulse.PhaseIncrement=HW.Constant.GoldenAngle138;
  Seq.Slice(1).CenterOfPulse=0;
  Seq.Slice(1).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
  Seq.Slice(1).GradTimeDelay=HW.Grad.SliceTimeDelay;
  Seq.Slice(1).UseAtRepetitionTime=1:numel(Seq.tRep);

  % generate slice parameters
  Seq=get_SliceParameter(Seq,HW);

  % Readout at tEcho

  % Seq.Read(1).alfa=linspace(0,2*pi*(1-1/numel(Seq.AQSlice(1).UsetRep)),numel(Seq.AQSlice(1).UsetRep))*0+Seq.AQSlice(1).alfa;
  % Seq.Read(1).phi=linspace(0,2*pi*(1-1/numel(Seq.AQSlice(1).UsetRep)),numel(Seq.AQSlice(1).UsetRep))*1+Seq.AQSlice(1).phi;
  % Seq.Read(1).theta=linspace(0,2*pi*(1-1/numel(Seq.AQSlice(1).UsetRep)),numel(Seq.AQSlice(1).UsetRep))*0+Seq.AQSlice(1).theta;
  % Seq.Read(1).sizeRead=sec(rem(linspace(0,2*pi*(1-1/numel(Seq.AQSlice(1).UsetRep)),numel(Seq.AQSlice(1).UsetRep))*1,pi/2))*Seq.AQSlice(1).sizeRead;
  % Seq.Read(1).sizeRead(
  Seq.Read(1).StartWithKSpaceCenter=1;
  Seq.Read(1).StartWithKSpaceCenterNonlinear=Seq.StartWithKSpaceCenterNonlinear;
  Seq.Read(1).StartWithGradientBeforeExcitationPulse=Seq.StartWithGradientBeforeExcitationPulse;
  Seq.Read(1).StartWithGradientTime=Seq.Slice(1).CenterOfPulse-Seq.Slice(1).GradLength/2;
  Seq.Read(1).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
  Seq.Read(1).CenterOfReadout=0;
  Seq.Read(1).PhaseIncrement=Seq.Slice(1).Pulse.PhaseIncrement;
  Seq.Read(1).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
  Seq.Read(1).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
  Seq.Read(1).GradTimeDelay=HW.Grad.ReadTimeDelay;
  Seq.Read(1).GradRephaseSign=1;
  Seq.Read(1).GradDephaseSign=-1;

  % generate readout timings
  Seq=get_ReadParameter(Seq,HW);

  Seq.Read(1).CenterOfDephase=[];
  Seq.Read(1).CenterOfRephase=[];
  if Seq.AQSlice(1).thickness<1000
    Seq.Read(1).CenterOfReadout=Seq.Slice(1).CenterOfPulse+Seq.Slice(1).GradLength/2+max(Seq.Slice(1).GradRephaseLength,Seq.Read(1).GradDephaseLength)-Seq.Read(1).StartOfReadoutGrad;
  else
    Seq.Read(1).CenterOfReadout=Seq.Slice(1).CenterOfPulse+Seq.Slice(1).Pulse.MaxLength/2-Seq.Read(1).StartOfReadoutGrad+Seq.Read(1).GradDephaseLength;
  end
  if Seq.StartWithKSpaceCenterNonlinear
    Seq.Read(1).CenterOfReadout=Seq.Read(1).CenterOfReadout-Seq.Read(1).GradDephaseLength+get_DeadTimeTX2RX(HW,Seq.Read(1).fSample);
  elseif Seq.StartWithGradientBeforeExcitationPulse
    Seq.Read(1).CenterOfReadout=Seq.Read(1).CenterOfReadout-Seq.Read(1).GradDephaseLength+get_DeadTimeTX2RX(HW,Seq.Read(1).fSample);
  end
  if isemptyfield(Seq, 'tEcho')
    Seq.tEcho = Seq.Read(1).CenterOfReadout;
  end

  Seq.Read(2).StartWithKSpaceCenter=1;
  Seq.Read(2).StartWithKSpaceCenterNonlinear=Seq.StartWithKSpaceCenterNonlinear;
  Seq.Read(2).StartWithGradientBeforeExcitationPulse=Seq.StartWithGradientBeforeExcitationPulse;
  Seq.Read(2).StartWithGradientTime=Seq.Slice(1).CenterOfPulse-Seq.Slice(1).GradLength/2;
  Seq.Read(2).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
  Seq.Read(2).CenterOfReadout=0;
  Seq.Read(2).PhaseIncrement=Seq.Slice(1).Pulse.PhaseIncrement;
  Seq.Read(2).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
  Seq.Read(2).GradTimeDelay=HW.Grad.ReadTimeDelay;
  Seq.Read(2).GradRephaseSign=1;
  Seq.Read(2).GradDephaseSign=-1;
  Seq.Read(2).CenterOfDephase=[];
  Seq.Read(2).CenterOfRephase=[];
  if Seq.AQSlice(1).thickness<1000
    Seq.Read(2).CenterOfReadout=Seq.Slice(1).CenterOfPulse+Seq.Slice(1).GradLength/2+max(Seq.Slice(1).GradRephaseLength,Seq.Read(1).GradDephaseLength)-Seq.Read(1).StartOfReadoutGrad;
  else
    Seq.Read(2).CenterOfReadout=Seq.Slice(1).CenterOfPulse+Seq.Slice(1).Pulse.MaxLength/2-Seq.Read(1).StartOfReadoutGrad+Seq.Read(1).GradDephaseLength;
  end
  if Seq.StartWithKSpaceCenterNonlinear
    Seq.Read(2).CenterOfReadout=Seq.Read(2).CenterOfReadout-Seq.Read(1).GradDephaseLength+get_DeadTimeTX2RX(HW,Seq.Read(1).fSample);
    temp=Seq.Read(1).nRead;
    Seq.Read(1).nRead=Seq.AQSlice(1).nRead; % reset nRead
    Seq.Read(2).nRead=Seq.AQSlice(1).nRead; % reset nRead
    Seq.AQSlice(1).nRead=temp;
  elseif Seq.StartWithGradientBeforeExcitationPulse
    Seq.Read(2).CenterOfReadout=Seq.Read(2).CenterOfReadout-Seq.Read(1).GradDephaseLength+get_DeadTimeTX2RX(HW,Seq.Read(1).fSample);
    temp=Seq.Read(1).nRead;
    Seq.Read(1).nRead=Seq.AQSlice(1).nRead; % reset nRead
    Seq.Read(2).nRead=Seq.AQSlice(1).nRead; % reset nRead
    Seq.AQSlice(1).nRead=temp;

  end

  Seq.Read(2).UseAtRepetitionTime=Seq.Read(1).UseAtRepetitionTime(1)-((Seq.SteadyState_PreShots:-1:1));
  Seq.Read(2).alfa=Seq.Read(1).alfa((Seq.SteadyState_PreShots:-1:1))*0+Seq.AQSlice(1).alfa;
  Seq.Read(2).phi=Seq.kgrid.wallphi((Seq.SteadyState_PreShots:-1:1)).'*1+Seq.AQSlice(1).phi;
  Seq.Read(2).theta=Seq.kgrid.walltheta((Seq.SteadyState_PreShots:-1:1)).'*1+Seq.AQSlice(1).theta;
  Seq.Read(2).sizeRead=Seq.kgrid.wallr((Seq.SteadyState_PreShots:-1:1)).'*(~Seq.AQSlice(1).FixedReadGrad)+max(Seq.kgrid.wallr(:))*(Seq.AQSlice(1).FixedReadGrad);
  Seq.Read(2).Resolution=Seq.Read(1).Resolution((Seq.SteadyState_PreShots:-1:1));

  % generate readout timings final
  Seq=get_ReadParameter(Seq,HW);



  % spoiler encoding aligned with read rephase
  Seq.Phase(1).sizePhase=mean(Seq.AQSlice(1).sizeRead/Seq.AQSlice(1).nRead);
  Seq.Phase(1).nPhase=1;
  Seq.Phase(1).PhaseOS=2;
  Seq.Phase(1).StepOrder=[1,1];
  Seq.Phase(1).CenterOfDephase=Seq.Read(1).CenterOfDephase;
  Seq.Phase(1).CenterOfRephase=Seq.Read(1).CenterOfRephase;
  Seq.Phase(1).GradDephaseLength=Seq.Read(1).GradDephaseLength;
  Seq.Phase(1).GradRephaseLength=Seq.Read(1).GradRephaseLength;
  Seq.Phase(1).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(1).GradTimeDelay=HW.Grad.PhaseTimeDelay;
  Seq.Phase(1).UseAtRepetitionTime=1:numel(Seq.tRep);


  % generate phase parameters
  Seq=get_PhaseParameter(Seq,HW);
  % Gradients of image encoding

  % % Readout at FID
  % Seq.AQSlice(2).plotkSpace=1;                               % Plot k-Space
  % Seq.AQSlice(2).plotImage=1;                                % Plot Image
  % Seq.AQSlice(2).plotPhase=2;                                % Plot phase of k-Space and Image
  % Seq.AQSlice(2).alfa=1.5*pi;                                % first rotation around x axis in RAD
  % Seq.AQSlice(2).phi=0.0*pi;                                 % second rotation around y axis in RAD
  % Seq.AQSlice(2).theta=-0.0*pi;                               % third rotation around z axis in RAD
  % Seq.AQSlice(2).HzPerPixMin=2000;                           % Bandwith per pixel in Hz
  % Seq.AQSlice(2).nRead=16;                                   % Number of Pixels in read
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
  % Seq.AQSlice(2).UsetRep=1:numel(Seq.tRep);
  % Seq.AQSlice(2).UseAQWindow=1;
  % Seq.AQSlice(2).ReadCoordinate=1;                                                   % direction of Read:  x = 1,  y = 2, z = 3
  % Seq.AQSlice(2).PhaseCoordinate=2;                                                  % direction of Phase:   x = 1,  y = 2, z = 3
  % Seq.AQSlice(2).Phase3DCoordinate=3;                                                % direction of Phase3D:   x = 1,  y = 2, z = 3
  % Seq.AQSlice(2).SliceCoordinate=3;                                                  % direction of Slice:   x = 1,  y = 2, z = 3
  % Seq.AQSlice(2).plotImageHandle=109;
  %
  % Seq.Read(2).useAQSlice=2;
  % Seq.Read(2).HzPerPixelMin=Seq.AQSlice(2).HzPerPixMin;
  % Seq.Read(2).CenterOfReadout=0.5/Seq.Read(2).HzPerPixelMin+Seq.Slice(1).CenterOfRephase+Seq.Slice(1).GradRephaseLength/2+HW.Grad.tEC;
  % Seq.Read(2).PhaseIncrement=Seq.Read(1).PhaseIncrement;
  % Seq.Read(2).UseAtRepetitionTime=Seq.AQSlice(2).UsetRep;
  % % generate readout timings
  % Seq=get_ReadParameter(Seq,HW);

  % AQ
  % AQ=add_AQ(Seq.Read(1).AQ, Seq.Read(2).AQ);
  % if Seq.CorrectPhase>0
  %     AQ=add_AQ(Seq.Read(2).AQ,Seq.Read(1).AQ);
  %     if Seq.Read(2).CenterOfReadout+Seq.Read(2).GradLength/2>Seq.Read(1).CenterOfDephase-Seq.Read(2).GradDephaseLength/2
  %         error('tEcho too short for selected AQ length (1/HzPerPixelMin)')
  %     end
  % else
       AQ=add_AQ(AQ, Seq.Read(1).AQ);
  %     if Seq.Slice(1).CenterOfRephase+Seq.Slice(1).GradRephaseLength/2>Seq.Read(1).CenterOfDephase-Seq.Read(2).GradDephaseLength/2
  %         error('tEcho too short for selected AQ length (1/HzPerPixelMin)')
  %     end
  %
  % end

  % TX
  TX=add_TX(TX, Seq.Slice(1).TX);

  if Seq.AQSlice(1).thickness<1000
    Grad=add_Grad(Grad, Seq.Slice(1).Grad);
    Grad=add_Grad(Grad, Seq.Slice(1).GradRephase);
  end
  if and(~Seq.StartWithKSpaceCenterNonlinear, ~Seq.StartWithGradientBeforeExcitationPulse)
    Grad=add_Grad(Grad, Seq.Read(1).GradDephase);
    Grad=add_Grad(Grad, Seq.Read(2).GradDephase);
  end
  Grad=add_Grad(Grad, Seq.Read(1).Grad);
  Grad=add_Grad(Grad, Seq.Read(2).Grad);
  if ~Seq.StartWithGradientBeforeExcitationPulse
    Grad=add_Grad(Grad, Seq.Read(1).GradRephase);
    Grad=add_Grad(Grad, Seq.Read(2).GradRephase);
    % Grad=add_Grad(Grad, Seq.Phase(1).GradRephase);
  end
  if Seq.LoopsBreakExactly
    if Loop==1
      Seq.Reinitialize=1;
      if Seq.Loops>1
        Seq.TimeToNextSequence=Seq.LoopsBreak;
        Seq.TimeFromLastSequence=[];
      else
        Seq.TimeToNextSequence=[];
        Seq.TimeFromLastSequence=[];
      end
    elseif Loop>=2 && Loop~=Seq.Loops
      Seq.Reinitialize=0;
      Seq.TimeFromLastSequence=Seq.LoopsBreak;
      Seq.TimeToNextSequence=Seq.LoopsBreak;
    elseif Loop==Seq.Loops
      Seq.Reinitialize=0;
      Seq.TimeFromLastSequence=Seq.LoopsBreak;
      Seq.TimeToNextSequence=[];
    end
  end

  if isempty(Seq.LoopsBreak)
    tb = 1;
  else
    tb = Seq.LoopsBreak;
  end
  if isempty(Seq.averageBreak)
    ta = 0;
  else
    ta = Seq.averageBreak;
  end
  if Seq.Loops>1
    disp(['Time to run remaining loops ' num2str((sum(Seq.tRep)*Seq.average+ta*(Seq.average-1))*(Seq.Loops-Loop+1)+tb*(Seq.Loops-Loop+1)) ' sec.']);
  end
  clear ta tb

  % Run sequence
  [ ~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);
  %     talker.mySequency.exctractArrayToFile(talker.mySequency.getCommandArray,'test.txt');
  if ~SeqOut.LoopsBreakExactly && ~isempty(SeqOut.LoopsBreak)
    init.Seq.StartSequenceTime = SeqOut.StartSequenceTime+SeqOut.SequenceTime+Seq.LoopsBreak;
  end
  % get data of FID
  % if Seq.CorrectPhase>0
  %    [dataFID]=get_kSpaceAndImage(data,SeqOut.AQSlice(2));
  %    if Seq.LoopPlot
  %      [dataFID]=plot_kSpaceAndImage(dataFID,SeqOut.AQSlice(2));
  %    end
  %    [data,SeqOut]=get_CorrectedPhase(data,SeqOut);
  % end
  % dataCorr=data.data;
  % data.data=data.dataWithoutCorrectPhase;
  % data.data=dataCorr;

  % FIXME: No image reconstruction currently
  % get data of echo
  % [data]=get_kSpaceAndImage(data,SeqOut.AQSlice(1));
  %  data.Image(max(abs(data.Image(:)))/4>abs(data.Image(:)))=nan;
  % if Seq.LoopPlot
  %   [data]=plot_kSpaceAndImage(data,SeqOut.AQSlice(1));
  % end
  if Loop == 1
    SeqLoop=SeqOut;
    SeqLoop.data=data;
    SeqLoop.data.data=zeros(size(data.data));
  end
  SeqLoop.dataLoop(Loop).data=data.data;
  if Seq.CorrectPhase>0
    SeqLoop.dataLoop(Loop).dataWithoutCorrectPhase = data.dataWithoutCorrectPhase;
  end
  SeqLoop.data.data=SeqLoop.data.data+data.data;
end
SeqLoop.AQSlice(1).plotImageHandle=2000;
SeqLoop.data.data=SeqLoop.data.data/SeqLoop.Loops;
% FIXME: No image reconstruction currently
% [SeqLoop.data]=get_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(1));
% [SeqLoop.data]=plot_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(1));


%%
% if 0;%isfield(data,'ImageCsi')
% %%
% index=[4,6,4];
% figure(2)
% plot(data.ImageCsiFrequency(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsi(:,index(1),index(2),index(3))))%%
% xlim([-20 20])
% %%
% index=[5,4,4];
% figure(3)
% plot(data.ImageCsiFrequencyZero(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsiRawZero(:,index(1),index(2),index(3))))
% xlim([-20 20])
% %%
% end
