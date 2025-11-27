% settings for piezo element SN #4
HW.Piezo.PiezoSN = 4;

HW.Piezo.Channel=               4;                         % grad output channel
HW.Piezo.Frequency=             1000;                      % piezo frequency in Hz
HW.Piezo.Start=                 0;                         % start time in s
HW.Piezo.DisplacementPhase=     0;                         % displacement phase at start time in degrees
HW.Piezo.Periods=               [];                        % piezo periods in n
HW.Piezo.Duration=              10e-3;                     % piezo duration in s
HW.Piezo.DisplacementPerVolt=   63e-6/150;                 % displacement in m/V
HW.Piezo.VoltageAmplitudeMax=   30;                        % max voltage at piezo in V
HW.Piezo.CurrentAmplitudeMax=   2;                         % max current in piezo in A
HW.Piezo.Capacity=              20.5e-6;                   % capacity of piezo in F
HW.Piezo.Resistance=            100;                       % parallel resistance to piezo to get 20 V offset voltage e.g. 100 Ohm
HW.Piezo.OffsetCurrent=         0.2;                       % offset current to get 20 V offset voltage e.g. 0.2 A
HW.Piezo.DisplacementPerAmpere = [];
HW.Piezo.Inductance = [];
HW.Piezo.OffsetVoltage = [];

HW.Piezo.OffsetOverdrivePercent = 0;
HW.Piezo.OffsetRampTime = 4.4e-3;
HW.Piezo.OffsetSetTime = 0.6e-3;
HW.Piezo.PreEmphasisPeriods=            10;
HW.Piezo.PreEmphasisOverdrivePercent=   20;
HW.Piezo.PreEmphasisRelativeAmplitude=  [ 0;  0;    (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;  (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;    1;1];
HW.Piezo.PreEmphasisRelativeTime=       [-1;  0;    1/3;                                              2/3;                                              1;2];
HW.Piezo.PostEmphasisPeriods=           2;
HW.Piezo.PostEmphasisOverdrivePercent=  -80;
HW.Piezo.PostEmphasisRelativeAmplitude= [ 1;    1;    (HW.Piezo.PostEmphasisOverdrivePercent+100)/100;    0;0];
HW.Piezo.PostEmphasisRelativeTime=      [-1;    0;    1/2;                                                 1;2];
HW.Piezo.UseAtRepetitionTime=           1;
HW.Piezo.Samples=                       1000;

HW.Grad.LoadIin2Amp(HW.Piezo.Channel) = 1;
HW.Grad.LoadRin(HW.Piezo.Channel) = 1;
if isscalar(HW.Grad.SystemTimeDelay)
  HW.Grad.SystemTimeDelay = ones(HW.Grad.n, 1) * HW.Grad.SystemTimeDelay;
end
HW.Grad.SystemTimeDelay(HW.Piezo.Channel) = 1e-06;         % time delay of grad amp at piezo channel
HW.Grad.CoilMaxDcCurrent(HW.Piezo.Channel) = 1*0.9;        % fuse DC current in A
HW.Grad.CoilCurrentSquareTime(HW.Piezo.Channel) = 0.9*0.9;  % fuse delay in A^2*sec


HW.Grad.AmpUnit{HW.Piezo.Channel} = 'A';
HW.Grad.AmpUnitScale(HW.Piezo.Channel) = 1/HW.Grad.LoadIin2Amp(HW.Piezo.Channel);
HW.Grad.Name{HW.Piezo.Channel} = 'piezo current';
