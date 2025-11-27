% settings for Stepper motor '14HS13-0804S' 

HW.Piezo.PiezoSN = 6;

HW.Piezo.Channel=               4;                         % Grad output channel

HW.Piezo.Frequency=             100;                      % Piezo Frequency in Hz
HW.Piezo.Start=                 0;                         % Start Time in s
HW.Piezo.DisplacementPhase=     0;                         % Displacement phase at Start Time in s
HW.Piezo.Periods=               [];                        % Piezo Periods in n
HW.Piezo.Duration=              10e-3;                     % Piezo Duration in s
HW.Piezo.Resistance=            100;                       % Parallel resistance to piezo to get 20 V offset voltage e.g. 100 Ohm
HW.Piezo.Inductance=            [];
HW.Piezo.OffsetCurrent=         0;                         % offset current  e.g. 0.1 A
HW.Piezo.DisplacementPerAmpere= 0;                         % Displacement in RAD/A
HW.Piezo.DisplacementPerVolt=   0;                         % Displacement in m/V
HW.Piezo.OffsetVoltage = [];
%         HW.Piezo.DisplacementPerAmpere= 0.032/0.5;                 % Displacement in RAD/A
HW.Piezo.DisplacementPerAmpere= 0.008*pi/360*1.8/2*1/0.7;  % Displacement in m/A
HW.Piezo.VoltageAmplitudeMax=   55;                        % Max Voltage at Piezo in V
HW.Piezo.CurrentAmplitudeMax=   0.8*2;                       % Max Current in Piezo in A
HW.Piezo.Capacity=              0;                         % Capacity of piezo
HW.Piezo.OffsetVoltage=         0;                         % offset voltage  e.g. 20 V
HW.Piezo.OffsetCurrent=         0;                         % offset voltage  e.g. 20 V
HW.Piezo.Resistance=            6.8;                        % resistance of Coil
HW.Piezo.Inductance=            10e-3+7.5e-3;                     % Inductance of Coil
HW.Grad.Inductance(4)=    HW.Piezo.Inductance;                     % Inductance of Coil
HW.Grad.PaUoutMax(HW.Piezo.Channel)=HW.Piezo.VoltageAmplitudeMax;
HW.Piezo.OffsetCurrent=         HW.Piezo.OffsetVoltage/HW.Piezo.Resistance;% offset current to get 20 V offset voltage e.g. 0.2 A

HW.Piezo.OffsetOverdrivePercent=        -20;
HW.Piezo.OffsetRampTime=                4.4e-3;
HW.Piezo.OffsetSetTime=                 0.6e-3;
HW.Piezo.OffsetBeforeAfterSeq=          1;
HW.Piezo.PreEmphasisPeriods=            10;
HW.Piezo.PreEmphasisOverdrivePercent=   00;
HW.Piezo.PreEmphasisRelativeAmplitude=  [ 0;  0;    (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;  (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;    1;1];
HW.Piezo.PreEmphasisRelativeTime=       [-1;  0;    1/3;                                              2/3;                                                1;2];
HW.Piezo.PostEmphasisPeriods=           2;
HW.Piezo.PostEmphasisOverdrivePercent=  -50;
HW.Piezo.PostEmphasisRelativeAmplitude= [ 1;    1;    (HW.Piezo.PostEmphasisOverdrivePercent+100)/100;    0;0];
HW.Piezo.PostEmphasisRelativeTime=      [-1;    0;    1/2;                                                 1;2];
HW.Piezo.UseAtRepetitionTime=           1;
HW.Piezo.Samples=                       900;

if isscalar(HW.Grad.SystemTimeDelay)
  HW.Grad.SystemTimeDelay = ones(HW.Grad.n, 1) * HW.Grad.SystemTimeDelay;
end
HW.Grad.SystemTimeDelay(HW.Piezo.Channel)=18.0e-06;              % Time delay of grad amp at piezo channel
HW.Grad.CoilMaxDcCurrent(HW.Piezo.Channel)=0.75*0.9*1;          % Fuse DC current in A
HW.Grad.CoilCurrentSquareTime(HW.Piezo.Channel)=0.9*0.9;      % Fuse delay in A^2*sec

HW.Grad.LoadIin2Amp(4)=HW.Piezo.DisplacementPerAmpere;
HW.Grad.LoadRin(4)=HW.Piezo.Resistance;

switch 'A' % Seq.plotSeq in ampere or meter
  case 'A'
    HW.Grad.AmpUnit{4}='A';
    HW.Grad.AmpUnitScale(4)=HW.Grad.LoadIin2Amp(4);
    HW.Grad.Name{4}='coil current';
  case 'm'
    HW.Grad.AmpUnit{4}=[char(181) 'm'];
    HW.Grad.AmpUnitScale(4)=1e-6;
    HW.Grad.Name{4}='displacement';
end

