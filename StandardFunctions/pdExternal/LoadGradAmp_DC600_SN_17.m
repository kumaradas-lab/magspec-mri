HW.Grad.ExtGradSN = 17;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.002585665,0.001740887,-0.000628206,0.002483215]; % Offset Current 17-May-2018 10:48:37
 
HW.Grad.PaUin2PaIout=([0.330494,0.3329757,0.3323808,0.3348418]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 17-May-2018 10:48:38
  
HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = 19.5e-6;                 % Time delay of grad amp in s

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600

