function [TXShape] = Pulse_LoadOwnPulseFromMat(HW, Center, BW, FlipAngle, maxpulseCount,  maxLength, FrequencyAdd, PhaseAdd)
%     Load own pulse from mat file
%
%     ComplexAmplitude                          % complex Amplitude in Hz of B1 fLarmor,        if Amplitude and Phase is not used
%     Amplitude                                 % abs of pulse Amplitude in Hz of B1 fLarmor,   if ComplexAmplitude is not used
%     Phase                                     % phase of pulse in deg,                        if ComplexAmplitude is not used
%     nPulseSegments;                           % Number of pulse segments,                     Set to [] if not used
%     DurationOfPulse;                          % Duration of shaped pulse,                     if not set it is sum of DurationOfPulseSegments
%     DurationOfPulseSegments;                  % Duration of pulse segments,                   if not set Segments have equale duration
%     CenterOfPulses                            % Center of NMR excitation                                      zero if not used
%     TimeOffsetToCenterOfPulse                 % Time from duration center to the center of NMR excitation     zero if not used
%     Bandwidth;                                % pulse bandwidth, needed for slice selectation                 zero if not used 

load DemoPulseData.mat

if ~exist('Amplitude','var');Amplitude=[];end;if isempty(Amplitude); Amplitude=abs(ComplexAmplitude); Phase=angle(ComplexAmplitude)/pi*180;end
if ~exist('nPulseSegments','var');nPulseSegments=[];end;if isempty(nPulseSegments); nPulseSegments=numel(Amplitude);end
if ~exist('Phase','var');Phase=[];end;if isempty(Phase); Phase=zeros(size(Amplitude));end
if ~exist('FrequencyOffset','var');FrequencyOffset=[];end;if isempty(FrequencyOffset); FrequencyOffset=zeros(size(Amplitude));end
if ~exist('DurationOfPulse','var');DurationOfPulse=[];end;if isempty(DurationOfPulse); DurationOfPulse=sum(DurationOfPulseSegments);end
if ~exist('DurationOfPulseSegments','var');DurationOfPulseSegments=[];end;if isempty(DurationOfPulseSegments); DurationOfPulseSegments=DurationOfPulse/nPulseSegments;end
if ~exist('ComplexAmplitude','var');ComplexAmplitude=[];end;if isempty(ComplexAmplitude); ComplexAmplitude=Amplitude.*exp(1i*Phase/180*pi);end
if ~exist('CenterOfPulses','var');CenterOfPulses=[];end;if isempty(CenterOfPulses); CenterOfPulses=0;end
if ~exist('TimeOffsetToCenterOfPulse','var');TimeOffsetToCenterOfPulse=[];end;if isempty(TimeOffsetToCenterOfPulse); TimeOffsetToCenterOfPulse=0;end
if ~exist('Bandwidth','var');Bandwidth=[];end;if isempty(Bandwidth); Bandwidth=0;end

if numel(DurationOfPulseSegments)==1; DurationOfPulseSegments=ones(size(Amplitude))*DurationOfPulseSegments;end
if numel(Phase)==1; Phase=ones(size(Amplitude))*Phase;end
if numel(FrequencyOffset)==1; FrequencyOffset=ones(size(Amplitude))*FrequencyOffset;end

% Set Bandwidth and Length of your pulse for the Slice calculation
if nargin==2
    if strcmp(Center,'Amp')                         % for 90 degrees rect pulse
        TXShape=nan;                                % PulseLength/HW.tFlip90Def  % what is the factor in amp of the pulse to get a 90 degrees flip in the same time as a 90 degrees rect pulse
    elseif strcmp(Center,'Time')
        TXShape=nan;                                % PulseLength * BW           % what is the factor in time of the pulse to get the same BW as an rect pulse  (BW of rect pulse = 1/length)
    elseif strcmp(Center,'Bandwidth')
        TXShape=Bandwidth;                          % BW, use to fix BW of pulse (set 'Time' to nan)
    elseif strcmp(Center,'Length')
        TXShape=DurationOfPulse;                    % Length, use to fix length of pulse  (set 'Amp' and 'Time' to nan)                    
    else
        TXShape=nan;
    end
else
    
    % put Shape values into TXShape
    TXShape(1).Start=CenterOfPulses-DurationOfPulse/2-TimeOffsetToCenterOfPulse+cumsum(DurationOfPulseSegments(:))-DurationOfPulseSegments(1)+Center;
    TXShape(1).Center=CenterOfPulses+Center;
    TXShape(1).CenterOffset=TimeOffsetToCenterOfPulse;
    TXShape(1).Duration=DurationOfPulseSegments(:);
    TXShape(1).Frequency=FrequencyOffset(:)+FrequencyAdd;
    TXShape(1).Phase=Phase(:)+PhaseAdd;
    TXShape(1).Amplitude=Amplitude(:)*2*pi/HW.Gamma.H1;
end

% generate a mat file
if 0
    nPulseSegments=100;                                                 % Number of pulse segments
    DurationOfPulse=200e-6;                                             % Duration of shaped pulse
    DurationOfPulseSegments=DurationOfPulse/nPulseSegments;             % Duration of pulse segment
    ComplexAmplitude=sinc(linspace(-2,2,nPulseSegments))*6000;                     % complex Amplitude in Hz of B1 fLarmor
    Amplitude=abs(ComplexAmplitude);                                    % abs of pulse Amplitude in Hz of B1 fLarmor
    Phase=angle(ComplexAmplitude)/pi*180;                               % phase of pulse in deg
    FrequencyOffset=0;                                                  % Frequency offset of pulse in deg
    CenterOfPulses=0;                                                   % Center of NMR excitation
    TimeOffsetToCenterOfPulse=0;                                        % Time from duration center to the center of NMR excitation
    Bandwidth=10e3;                                                     % pulse bandwidth. Set to zero if not used 
    save('DemoPulseData.mat',...
         'nPulseSegments',...
         'DurationOfPulse',...
         'DurationOfPulseSegments',...
         'ComplexAmplitude',...
         'Amplitude',...
         'Phase',...
         'FrequencyOffset',...
         'CenterOfPulses',...
         'TimeOffsetToCenterOfPulse',...
         'Bandwidth')
end

% translate a PulseData.mat file
if 0
    load PulseData.mat                                      % load pulse data with amp (Amplitude) and pha (Phase)
    Amplitude=amp;                                          % Amplitude in Hz of B1 fLarmor
    Phase=pha+90;                                              % Phase in deg

    DurationOfPulse=0.5e-3;                                 % Set duration of the shaped pulse
    
    CenterOfPulses=0;                                       % Center of NMR excitation
    TimeOffsetToCenterOfPulse=240e-6;                       % Time from duration center to the center of NMR excitation

    Bandwidth=20e3;                                         % pulse bandwidth. Set to zero if not used 
    
    save('PulseData2.mat',...
         'DurationOfPulse',...
         'Amplitude',...
         'Bandwidth',...
         'CenterOfPulses',...
         'TimeOffsetToCenterOfPulse',...
         'Phase')
end
