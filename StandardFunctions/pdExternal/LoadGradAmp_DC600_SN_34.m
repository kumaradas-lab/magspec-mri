HW.Grad.ExtGradSN = 34;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.002199488,0.002812103,0.001010186,0.001241282]; % Offset Current 15-Mar-2021 15:42:13
 
HW.Grad.PaUin2PaIout=([0.3345515,0.3352553,0.333306,0.3337563]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 15-Mar-2021 15:42:13

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;

switch HW.UserName
  case 'magnet_01_probe_1H'
    % 10mm Probenkopf H1
    HW.Grad.SystemTimeDelay(1:3) = [2.21384e-05  2.53016e-05  2.69616e-05]; % Time delay of grad amp
  case 'magnet_02_probe_1H'
    % 15mm Probenkopf H1
    HW.Grad.SystemTimeDelay(1:3) = [3.43072e-05  3.87304e-05  4.58992e-05]; % Time delay of grad amp
  otherwise
    % 10mm Probenkopf H1
    HW.Grad.SystemTimeDelay(1:3) = [2.21384e-05  2.53016e-05  2.69616e-05]; % Time delay of grad amp
end
    
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
