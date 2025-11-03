HW.Grad.ExtGradSN = 10;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.0002897053,0.0005113268,-0.001826162,-0.001214733]; % Offset Current 07-Nov-2016 15:07:18

HW.Grad.PaUin2PaIout=([0.3339845,0.3331593,0.3335963,0.3334057]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 07-Nov-2016 15:08:35

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = 19.5e-6;                 % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600

% load MRE Setting
HW.Grad.LoadRin(4)=100;  %Resistance parallel to Piezo
% set correct sequence plot names
HW.Grad.Name(4)={'Piezo Current'};
HW.Grad.AmpUnit(4)={'A'};
HW.Grad.AmpUnitScale(4)=1;
% set correct gradient efficiency
HW.Grad.LoadIin2Amp(4)=1;
