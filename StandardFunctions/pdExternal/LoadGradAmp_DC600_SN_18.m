HW.Grad.ExtGradSN = 18;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.002307707,-0.001242077,0.003051957,-0.001274058]; % Offset Current 17-May-2018 12:26:58
 
HW.Grad.PaUin2PaIout=([0.3341128,0.3316297,0.3358662,0.3317328]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 17-May-2018 12:26:58
   
HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [5.2656e-05    6.272e-05   6.1616e-05]; % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600

