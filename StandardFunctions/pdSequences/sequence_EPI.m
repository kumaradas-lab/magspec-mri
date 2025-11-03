%% EPI experiment
%
% Disclaimer: This function is not actively supported at the moment.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%%
warning('PD:sequence_EPI:deprecated', ...
  ['This script is deprecated and might be removed in a future version ', ...
   'of OpenMatlab. Consider using the function "sequence_Flash" with ', ...
   'Seq.AQSlice(1).EPIFactor>1 instead.']);

if exist('HW', 'var') && isa(HW, 'PD.HW')
  ResetStructs;  % reset variables: Seq AQ TX Grad
else
  % mySave.DummySerial = 0;   % for Dummy measurements without an MRI device (0 or [] for real measurements or 1 for dummy measurements)
  LoadSystem;  % load system parameters (reset all variables HW Seq AQ TX Grad)
end

iDevice = 1;  % FIXME: Support multiple MMRT devices
Channel = 1;  % FIXME: Support multiple acquisition channels?

HW.FindFrequencyPause = 0.5;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 100);   % Find magnet frequency
% HW.Grad(iDevice).TimeDelay=52e-6-15e-6;
HW.Grad(iDevice).SliceTimeDelay=[0,0,0]*1e-6+-0e-6;
HW.Grad(iDevice).ReadTimeDelay=[0,0,0]*1e-6;
HW.Grad(iDevice).PhaseTimeDelay=[0,0,0]*1e-6;
HW.Grad(iDevice).tRamp=18e-6;
HW.Grad(iDevice).tEC=18e-6;

Seq.average=1;                                          % Number of averages    1...
Seq.averageBreak=0.2;                                   % Pause between two averages in seconds

% Seq.tRep=10e-3;                                         % Repetition time in seconds
% Seq.T1=100e-3;                                          % T1 of probe excitation is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho=1.3e-3;                                       % Echo time in seconds
Seq.CLTime=0;                                           % command load time in seconds

Seq.plotSeqTR=[1:3];                                      % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
Seq.plotSeq=1:3;                                        % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)

% get to steady state
Seq.SteadyState_PreShots=1;
Seq.SteadyState_PostShots=1;


% Readout at echo
% rotate all image coordinates around... rotate
Seq.AQSlice(1).plotkSpace=1;                               % Plot k-Space 
Seq.AQSlice(1).plotImage=1;                                % Plot Image 
Seq.AQSlice(1).plotPhase=1;                                % Plot phase of k-Space and Image 
Seq.AQSlice(1).alfa=0.0*pi;                                % first rotation around x axis in RAD
Seq.AQSlice(1).phi=0.0*pi;                                 % second rotation around y axis in RAD
Seq.AQSlice(1).theta=0.0*pi;                               % third rotation around z axis in RAD
Seq.AQSlice(1).HzPerPixMin=2000;                           % Bandwith per pixel in Hz
Seq.AQSlice(1).nRead=16;                                   % Number of Pixels in read,  nPhase(1)=3 if nRead>1
Seq.AQSlice(1).nPhase(1)=8;                                % Number of Pixels in phase2D
Seq.AQSlice(1).nPhase(2)=4;                                % Number of Pixels in phase3D
Seq.AQSlice(1).nPhase(3)=1;                                % Number of Pixels in CSI, nRead=1 if nPhase(1)>1
Seq.AQSlice(1).sizeRead=0.010+0e12+(Seq.AQSlice(1).nRead==1)*1e12;                            % Image size in read in meter
Seq.AQSlice(1).sizePhase(1)=0.010+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).sizePhase(2)=0.010+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).sizePhase(3)=0.010+0e12;                           % Image size in phase in meter
Seq.AQSlice(1).thickness=0.01+1e16;                       % Image thickness in slice direction  used for 2D and 3D! in meter
Seq.AQSlice(1).ReadOSUsedForImage=17;                     % Number of samples used in CSI
Seq.AQSlice(1).ReadOS=4;                                   % Oversampling read         (recommended >=2) 1...
Seq.AQSlice(1).PhaseOS(1)=1;                                  % Oversampling phase        1...
Seq.AQSlice(1).PhaseOS(2)=1;                                  % Oversampling phase CSI        1...
Seq.AQSlice(1).PhaseOS(3)=1;                                % Oversampling phase3D      1...

tReps= cumsum(     ones(1+Seq.SteadyState_PreShots+Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)+Seq.SteadyState_PreShots,Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2)),1) ...
     +(cumsum(     ones(1+Seq.SteadyState_PreShots+Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)+Seq.SteadyState_PreShots,Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2)) ...
                   .*(1+Seq.SteadyState_PreShots+Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)+Seq.SteadyState_PreShots),2)...
                   - (1+Seq.SteadyState_PreShots+Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)+Seq.SteadyState_PreShots));

Seq.AQSlice(1).UsetRep = reshape(tReps(1+Seq.SteadyState_PreShots+(1:prod(Seq.AQSlice(1).nPhase(1))*prod(Seq.AQSlice(1).PhaseOS(1))),:), [], 1);
Seq.AQSlice(1).UseAQWindow=1;
Seq.AQSlice(1).AQWindowDirection=[1];
Seq.AQSlice(1).ReadCoordinate=3;                                                   % direction of Read:  x = 1,  y = 2, z = 3  
Seq.AQSlice(1).PhaseCoordinate(1)=1;                                                  % direction of Phase:   x = 1,  y = 2, z = 3  
Seq.AQSlice(1).PhaseCoordinate(2)=2;                                                  % direction of Phase:   x = 1,  y = 2, z = 3  
Seq.AQSlice(1).PhaseCoordinate(3)=3;                                                  % direction of Phase:   x = 1,  y = 2, z = 3  
Seq.AQSlice(1).SliceCoordinate=1;                                                  % direction of Slice:   x = 1,  y = 2, z = 3  

% Generate tRep times
Seq.tRep=Seq.tEcho+0*tReps;
Seq.tRep(end,:)=Seq.tEcho+HW.FindFrequencyPause;
Seq.tRep(end,end)=Seq.tEcho;
Seq.tRep=reshape(Seq.tRep,1,[]);

% 90 degrees slice excitation
% firstLength=[];
for t=1:2
     if isfield(Seq, 'Slice') 
        Seq=rmfield( Seq, 'Slice'); 
     end;
    Seq.Slice(1).Pulse.Function=@Pulse_Rect;
    Seq.Slice(1).Pulse.FlipAngle=90;
    Seq.Slice(1).Thickness=Seq.AQSlice(1).thickness;
    Seq.Slice(1).Pulse.Phase=0;
    Seq.Slice(1).Pulse.PhaseIncrement=0;
    Seq.Slice(1).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
    Seq.Slice(1).GradTimeDelay=HW.Grad(iDevice).SliceTimeDelay;
    Seq.Slice(1).UseAtRepetitionTime=reshape(tReps(1,:),1,[]);
    Seq.Slice(1).CenterOfPulse=Seq.tRep(1)/2;
%     Seq.Slice(1).CenterOfRephase=Seq.tRep(1)*3/4;
%     Seq.Slice(1).GradRephaseLength=Seq.tRep(1)*1/4-200e-6;
    % generate slice parameters
    Seq=get_SliceParameter(Seq,HW);

    Seq.tRep(1)=(Seq.Slice(1).GradRephaseLength)*2+Seq.Slice(1).GradLength+Seq.Slice(1).tRamp*2+100e-6;
end


% Readout at tRep/2
Seq.Read(1).useAQSlice=1;
Seq.Read(1).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
Seq.Read(1).CenterOfReadout=Seq.tEcho/2;
Seq.Read(1).PhaseIncrement=0;
Seq.Read(1).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
Seq.Read(1).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep(1:2:end);
Seq.Read(1).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Read(1).GradDephaseSign=-1;
Seq.Read(1).GradSign=1;
Seq.Read(1).GradRephaseSign=-1;

Seq.Read(2).useAQSlice=1;
Seq.Read(2).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
Seq.Read(2).CenterOfReadout=Seq.tEcho/2;
Seq.Read(2).PhaseIncrement=0;
Seq.Read(2).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
Seq.Read(2).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep(2:2:end);
Seq.Read(2).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Read(2).GradDephaseSign=-1;
Seq.Read(2).GradSign=1;
Seq.Read(2).GradRephaseSign=-1;



% generate readout timings
Seq=get_ReadParameter(Seq,HW);

% First Readout Phase prepare
Seq.Read(3).useAQSlice=1;
Seq.Read(3).HzPerPixelMin=Seq.AQSlice(1).HzPerPixMin;
Seq.Read(3).CenterOfReadout=0;
Seq.Read(3).CenterOfDephase=Seq.Read(1).CenterOfRephase;
Seq.Read(3).GradDephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Read(3).PhaseIncrement=0;
Seq.Read(3).UseCoordinate=Seq.AQSlice(1).ReadCoordinate;
Seq.Read(3).UseAtRepetitionTime=reshape(tReps(1+Seq.SteadyState_PreShots,:),1,[]);
Seq.Read(3).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Read(3).GradRephaseSign=1;
Seq.Read(3).GradDephaseSign=-1;
% generate readout timings
Seq=get_ReadParameter(Seq,HW);



% phase encoding aligned with read dephase and rephase
Seq.Phase(1).sizePhase=Seq.AQSlice(1).sizePhase(2)/2*Seq.AQSlice(1).PhaseOS(2);
Seq.Phase(1).nPhase=1;
Seq.Phase(1).PhaseOS=2;
Seq.Phase(1).StepOrder=[1,1];
Seq.Phase(1).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(1).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(1).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(1).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(1).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(2);
Seq.Phase(1).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(1).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
Seq.Phase(1).GradRephaseSign=1;
Seq.Phase(1).GradDephaseSign=-1;

Seq.Phase(2).sizePhase=Seq.AQSlice(1).sizePhase(2);
Seq.Phase(2).nPhase=Seq.AQSlice(1).nPhase(2);
Seq.Phase(2).PhaseOS=Seq.AQSlice(1).PhaseOS(2);
Seq.Phase(2).StepIncrement=1;
Seq.Phase(2).CenterOfDephase=Seq.Read(1).CenterOfReadout;
Seq.Phase(2).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(2).GradDephaseLength=Seq.Read(1).GradDephaseLength+Seq.Read(1).GradRephaseLength+Seq.Read(1).GradLength;
Seq.Phase(2).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(2).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(2);
Seq.Phase(2).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(2).UseAtRepetitionTime=Seq.Read(3).UseAtRepetitionTime;
Seq.Phase(2).GradDephaseSign=1;

% phase encoding aligned with read dephase and rephase
Seq.Phase(3).sizePhase=Seq.AQSlice(1).sizePhase(1)/2*Seq.AQSlice(1).PhaseOS(1);
Seq.Phase(3).nPhase=1;
Seq.Phase(3).PhaseOS=2;
Seq.Phase(3).StepOrder=[1,1];
Seq.Phase(3).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(3).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(3).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(3).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(3).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(1);
Seq.Phase(3).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(3).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
Seq.Phase(3).GradRephaseSign=1;
Seq.Phase(3).GradRephaseSign=-1;

Seq.Phase(4).sizePhase=Seq.AQSlice(1).sizePhase(1)/Seq.AQSlice(1).nPhase(1);
Seq.Phase(4).nPhase=1;
Seq.Phase(4).PhaseOS=2;
Seq.Phase(4).StepOrder=[1,1];
Seq.Phase(4).CenterOfDephase=Seq.Read(1).CenterOfReadout;
Seq.Phase(4).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(4).GradDephaseLength=Seq.Read(1).GradDephaseLength+Seq.Read(1).GradRephaseLength+Seq.Read(1).GradLength;
Seq.Phase(4).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(4).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(1);
Seq.Phase(4).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(4).UseAtRepetitionTime=Seq.Read(3).UseAtRepetitionTime;
Seq.Phase(4).GradDephaseSign=+1;

% phase encoding aligned with read dephase and rephase
Seq.Phase(5).sizePhase=Seq.AQSlice(1).sizePhase(3)/2*Seq.AQSlice(1).PhaseOS(3);
Seq.Phase(5).nPhase=1;
Seq.Phase(5).PhaseOS=2;
Seq.Phase(5).StepOrder=[1,1];
Seq.Phase(5).CenterOfDephase=Seq.Read(1).CenterOfDephase;
Seq.Phase(5).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(5).GradDephaseLength=Seq.Read(1).GradDephaseLength;
Seq.Phase(5).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(5).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(3);
Seq.Phase(5).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(5).UseAtRepetitionTime=Seq.AQSlice(1).UsetRep;
Seq.Phase(5).GradRephaseSign=1;
Seq.Phase(5).GradRephaseSign=-1;

Seq.Phase(6).sizePhase=Seq.AQSlice(1).sizePhase(3)./Seq.AQSlice(1).nPhase(3);
Seq.Phase(6).nPhase=1;
Seq.Phase(6).PhaseOS=2;
Seq.Phase(6).StepOrder=[1,1];
Seq.Phase(6).CenterOfDephase=Seq.Read(1).CenterOfReadout;
Seq.Phase(6).CenterOfRephase=Seq.Read(1).CenterOfRephase;
Seq.Phase(6).GradDephaseLength=Seq.Read(1).GradDephaseLength+Seq.Read(1).GradRephaseLength+Seq.Read(1).GradLength;
Seq.Phase(6).GradRephaseLength=Seq.Read(1).GradRephaseLength;
Seq.Phase(6).UseCoordinate=Seq.AQSlice(1).PhaseCoordinate(3);
Seq.Phase(6).GradTimeDelay=HW.Grad(iDevice).ReadTimeDelay;
Seq.Phase(6).UseAtRepetitionTime=Seq.Read(3).UseAtRepetitionTime;
Seq.Phase(6).GradDephaseSign=+1;

% generate slice parameters
% Seq=get_SliceParameter(Seq,HW);
% generate phase parameters
Seq=get_PhaseParameter(Seq,HW);

% Seq.Phase(1).GradDephase(1).Amp
% Seq.Phase(2).GradDephase(1).Amp
% Seq.Phase(3).GradDephase(1).Amp
% Seq.Phase(4).GradDephase(1).Amp
% Seq.Phase(5).GradDephase(1).Amp
% Seq.Phase(6).GradDephase(1).Amp


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
% Seq.AQSlice(2).UsetRep=1:length(Seq.tRep);
% Seq.AQSlice(2).UseAQWindow=1;
% Seq.AQSlice(2).ReadCoordinate=1;                                                   % direction of Read:  x = 1,  y = 2, z = 3  
% Seq.AQSlice(2).PhaseCoordinate=2;                                                  % direction of Phase:   x = 1,  y = 2, z = 3  
% Seq.AQSlice(2).Phase3DCoordinate=3;                                                % direction of Phase3D:   x = 1,  y = 2, z = 3  
% Seq.AQSlice(2).SliceCoordinate=3;                                                  % direction of Slice:   x = 1,  y = 2, z = 3  
% 
% Seq.Read(3).HzPerPixelMin=Seq.AQSlice(2).HzPerPixMin;
% Seq.Read(3).CenterOfReadout=0.5/Seq.Read(3).HzPerPixelMin+Seq.Slice(1).CenterOfRephase+Seq.Slice(1).GradRephaseLength/2+HW.Grad.tEC;
% Seq.Read(3).PhaseIncrement=Seq.Read(1).PhaseIncrement;
% % generate readout timings
% Seq=get_ReadParameter(Seq,HW);

% AQ
% AQ=add_AQ(Seq.Read(1).AQ, Seq.Read(2).AQ);
 AQ=add_AQ(AQ,Seq.Read(1).AQ);
 AQ=add_AQ(AQ,Seq.Read(2).AQ);
% AQ=Seq.Read(1).AQ;
% TX
TX=add_TX(TX, Seq.Slice(1).TX);

% Gradients of image encoding
Grad=add_Grad(Grad, Seq.Slice(1).Grad);
Grad=add_Grad(Grad, Seq.Slice(1).GradRephase);

% % Grad=add_Grad(Grad, Seq.Read(3).GradDephase); 


% Grad=add_Grad(Grad, Seq.Phase(1).GradDephase);
% Grad=add_Grad(Grad, Seq.Phase(2).GradDephase);
% Grad=add_Grad(Grad, Seq.Phase(3).GradDephase);
 Grad=add_Grad(Grad, Seq.Read(1).GradDephase);
 Grad=add_Grad(Grad, Seq.Read(1).Grad);
 Grad=add_Grad(Grad, Seq.Read(1).GradRephase);
 Grad=add_Grad(Grad, Seq.Read(2).GradDephase);
 Grad=add_Grad(Grad, Seq.Read(2).Grad);
 Grad=add_Grad(Grad, Seq.Read(2).GradRephase);
% if double(Seq.AQSlice(1).nPhase(1)>1); Grad=add_Grad(Grad, Seq.Phase(1).GradDephase);end


% if double(Seq.AQSlice(1).nPhase(1)>1); Grad=add_Grad(Grad, Seq.Phase(1).GradRephase);end
if double(Seq.AQSlice(1).nPhase(1)>1); Grad=add_Grad(Grad, Seq.Phase(2).GradDephase);end
% if double(Seq.AQSlice(1).nPhase(2)>1); Grad=add_Grad(Grad, Seq.Phase(3).GradDephase);end
if double(Seq.AQSlice(1).nPhase(2)>1); Grad=add_Grad(Grad, Seq.Phase(3).GradRephase);end
if double(Seq.AQSlice(1).nPhase(2)>1); Grad=add_Grad(Grad, Seq.Phase(4).GradDephase);end
% if double(Seq.AQSlice(1).nPhase(3)>1); Grad=add_Grad(Grad, Seq.Phase(5).GradDephase);end
% if double(Seq.AQSlice(1).nPhase(3)>1); Grad=add_Grad(Grad, Seq.Phase(5).GradRephase);end
% if double(Seq.AQSlice(1).nPhase(3)>1); Grad=add_Grad(Grad, Seq.Phase(6).GradDephase);end
% Grad(1).Shim=0.7e-3;
% Grad=add_Grad(Grad, Seq.Phase(4).GradRephase);
% Seq.Loops=1;
% Seq.average=1;
% for Loop=1:Seq.Loops
% if Loop==1; Seq.Reinitialize=1;end
% if Seq.Loops>1; Seq.TimeToNextSequence=0.5; else Seq.TimeToNextSequence=[]; end
% if Loop==2;
%     Seq.Reinitialize=0; 
%     Seq.TimeFromLastSequence=Seq.TimeToNextSequence;
% end
% if Loop==Seq.Loops;
%     Seq.TimeFromLastSequence=Seq.TimeToNextSequence;
%     Seq.TimeToNextSequence=[];
% end

% Run sequence
[ ~, SeqOut, data, ~] = set_sequence(HW, Seq, AQ, TX, Grad);

% % get data of FID
%  [dataFID]=get_kSpaceAndImage(data,SeqOut.AQSlice(2));
%  [dataFID]=plot_kSpaceAndImage(dataFID,SeqOut.AQSlice(2));
% %  SeqOut.AQSlice(1).ReadOSUsedForImage=17;
% % get data of echo
% SeqOut.AQSlice(1).PhaseOS(2)=1;                                  % Oversampling phase CSI        1...
% SeqOut.AQSlice(1).UsetRep=SeqOut.AQSlice(1).UsetRep(1:2:end);

iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

 [data]=get_kSpaceAndImage(data(iAQ),SeqOut.AQSlice(1));
%  data.Image(max(abs(data.Image(:)))/4>abs(data.Image(:)))=nan;
 [data]=plot_kSpaceAndImage(data,SeqOut.AQSlice(1));
 data.SeqOut=SeqOut;
% % end
%%
if 0;%isfield(data,'ImageCsi')
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
