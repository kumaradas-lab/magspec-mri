
HW.Piezo.Channel=               4;                         % Grad output channel
if HW.Piezo.Channel==3
  HW.Grad.x=2;
  HW.Grad.y=4;
  HW.Grad.z=1;
  HW.Grad.B=3;
  if isstruct(HW);HW.Grad.xyzB=[HW.Grad.x,HW.Grad.y,HW.Grad.z,HW.Grad.B];    end           % index Gadienten Channal Zuordnung.
  if isstruct(HW);[~, HW.Grad.Channel2xyzB] =sort(HW.Grad.xyzB);  end  % Dac Channal zu Gradienten Richtung
else
  HW.Grad.x=2;
  HW.Grad.y=3;
  HW.Grad.z=1;
  HW.Grad.B=4;
  if isstruct(HW); HW.Grad.xyzB=[HW.Grad.x,HW.Grad.y,HW.Grad.z,HW.Grad.B];   end           % index Gadienten Channal Zuordnung.
  if isstruct(HW);[~, HW.Grad.Channel2xyzB] =sort(HW.Grad.xyzB);  end  % Dac Channal zu Gradienten Richtung
end

HW.Piezo.Frequency=             1000;                      % Piezo Frequency in Hz
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
switch 'PK2FVP2'
    case 'PK44LA2P2' % rund mit Loch
        HW.Piezo.DisplacementPerVolt=   9e-6/150;                  % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   10;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   2;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              2e-6;                      % Capacity of piezo
        HW.Piezo.OffsetVoltage=         HW.Piezo.VoltageAmplitudeMax/2+5;  % offset voltage  e.g. 20 V
    case 'PK2FVP2'  % 75V 5x5x40mm
        HW.Piezo.DisplacementPerVolt=   44.8e-6/75;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   32;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   2;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              16.5e-6;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         HW.Piezo.VoltageAmplitudeMax-3;  % offset voltage  e.g. 20 V
    case 'Pahl60/20'
        HW.Piezo.DisplacementPerVolt=   63e-6/150;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   35;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   4;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              21.5e-6;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         20;                        % offset voltage  e.g. 20 V
    case 'ccc'  
        HW.Piezo.DisplacementPerVolt=   44.8e-6/75;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   32;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   2;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              16.5e-6;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         0;  % offset voltage  e.g. 20 V
    case 'JK35HS280504' % Stepper motor
        HW.Piezo.DisplacementPerAmpere= 0.032/0.5;                 % Displacement in RAD/A
        HW.Piezo.VoltageAmplitudeMax=   60;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   0.6;                       % Max Current in Piezo in A
        HW.Piezo.Capacity=              0;                         % Capacity of piezo
        HW.Piezo.OffsetVoltage=         0;                         % offset voltage  e.g. 20 V
        HW.Piezo.OffsetCurrent=         0;                         % offset voltage  e.g. 20 V
        HW.Piezo.Resistance=            20;                        % resistance of Coil
        HW.Piezo.Inductance=            20e-3;                     % Inductance of Coil
    case '14HS13-0804S' % Stepper motor
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
    case 'PKP225D15B2' % Stepper motor
%         HW.Piezo.DisplacementPerAmpere= 0.032/0.5;                 % Displacement in RAD/A
        HW.Piezo.DisplacementPerAmpere= 0.008*pi/360*1.8/2*1/0.7;%?  % Displacement in m/A
        HW.Piezo.VoltageAmplitudeMax=   55;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   0.8*2;                       % Max Current in Piezo in A
        HW.Piezo.Capacity=              0;                         % Capacity of piezo
        HW.Piezo.OffsetVoltage=         0;                         % offset voltage  e.g. 20 V
        HW.Piezo.OffsetCurrent=         0;                         % offset voltage  e.g. 20 V
        HW.Piezo.Resistance=            2*557/403*557/500;%?                        % resistance of Coil
        HW.Piezo.Inductance=            3.5e-3;%?                     % Inductance of Coil
        HW.Grad.Inductance(4)=    HW.Piezo.Inductance;                     % Inductance of Coil
        HW.Grad.PaUoutMax(HW.Piezo.Channel)=HW.Piezo.VoltageAmplitudeMax;
  case 'TestCs'
        HW.Piezo.DisplacementPerVolt=   63e-6/150;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   35;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   2;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              9.5e-6;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         20;                        % offset voltage  e.g. 20 V
    case '7.6 Ohm'
        HW.Piezo.DisplacementPerVolt=   1;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   62;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   7.4;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              0;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         0;                        % offset voltage  e.g. 20 V
        HW.Piezo.Resistance=            13.6;                     % Parallel resistance to piezo to get 20 V offset voltage e.g. 100 Ohm
    otherwise
        HW.Piezo.DisplacementPerVolt=   15e-6/150;                 % Displacement in m/V
        HW.Piezo.VoltageAmplitudeMax=   1;                        % Max Voltage at Piezo in V
        HW.Piezo.CurrentAmplitudeMax=   1;                         % Max Current in Piezo in A
        HW.Piezo.Capacity=              1e-6;                   % Capacity of piezo
        HW.Piezo.OffsetVoltage=         20;                        % offset voltage  e.g. 20 V
end
HW.Piezo.OffsetCurrent=         HW.Piezo.OffsetVoltage/HW.Piezo.Resistance;% offset current to get 20 V offset voltage e.g. 0.2 A

HW.Piezo.OffsetOverdrivePercent=        -20;
HW.Piezo.OffsetRampTime=                4.4e-3;
HW.Piezo.OffsetSetTime=                 0.6e-3;
HW.Piezo.OffsetBeforeAfterSeq=          1;
HW.Piezo.PreEmphasisPeriods=            20;
HW.Piezo.PreEmphasisOverdrivePercent=   00;
HW.Piezo.PreEmphasisRelativeAmplitude=  [ 0;  0;    (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;  (HW.Piezo.PreEmphasisOverdrivePercent+100)/100;    1;1];
HW.Piezo.PreEmphasisRelativeTime=       [-1;  0;    1/3;                                              2/3;                                                1;2];
HW.Piezo.PostEmphasisPeriods=           2;
HW.Piezo.PostEmphasisOverdrivePercent=  -50;
HW.Piezo.PostEmphasisRelativeAmplitude= [ 1;    1;    (HW.Piezo.PostEmphasisOverdrivePercent+100)/100;    0;0];
HW.Piezo.PostEmphasisRelativeTime=      [-1;    0;    1/2;                                                 1;2];
HW.Piezo.UseAtRepetitionTime=           1;
HW.Piezo.Samples=                       900;

if numel(HW.Grad.SystemTimeDelay)==1; HW.Grad.SystemTimeDelay=ones(HW.Grad.n,1)*HW.Grad.SystemTimeDelay; end
HW.Grad.SystemTimeDelay(HW.Piezo.Channel)=18.0e-06;              % Time delay of grad amp at piezo channel
HW.Grad.CoilMaxDcCurrent(HW.Piezo.Channel)=0.75*0.9*1;          % Fuse DC current in A
HW.Grad.CoilCurrentSquareTime(HW.Piezo.Channel)=0.9*0.9;      % Fuse delay in A^2*sec

if HW.Piezo.DisplacementPerVolt
  if HW.Grad.PaCurrentControlled(HW.Piezo.Channel)
      HW.Grad.LoadIin2Amp(4)=1;
      HW.Grad.LoadRin(4)=HW.Piezo.Resistance;
      HW.Grad.AmpUnit{4}='A';
      HW.Grad.AmpUnitScale(4)=1/HW.Grad.LoadIin2Amp(4);
      HW.Grad.Name{4}='piezo current';
  else
      HW.Grad.LoadRin(4)=HW.Piezo.Resistance;
      HW.Grad.LoadIin2Amp(4)=HW.Grad.LoadRin(4)*HW.Piezo.DisplacementPerVolt;
      HW.Grad.AmpUnit{4} = [char(181) 'm'];
      HW.Grad.AmpUnitScale(4)=1e-6;
      HW.Grad.Name{4}='piezo displacement';
  end
else
  if HW.Grad.PaCurrentControlled(HW.Piezo.Channel)
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
%       case 'rad'
%         HW.Grad.AmpUnit{4}='RAD';
%         HW.Grad.AmpUnitScale(4)=1/HW.Grad.LoadIin2Amp(4);
%         HW.Grad.Name{4}='coil displacement';
      end
  else
      HW.Grad.AmpCurrentDependent(HW.Piezo.Channel)=0;
      HW.Grad.LoadRin(4)=HW.Piezo.Resistance;
      HW.Grad.LoadUin2Amp(4)=1;
      HW.Grad.AmpUnit{4} = 'V';
      HW.Grad.AmpUnitScale(4)=1;
      HW.Grad.Name{4}='coil voltage';
  end
end
