HW.Grad.ExtGradSN = 30;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.000342823,0.00331074,-0.0002937403,-0.002939028]; % Offset Current 14-May-2020 14:44:27

HW.Grad.PaUin2PaIout=([0.3317603,0.3349532,0.3327892,0.329336]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 14-May-2020 14:44:27

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6]; % Time delay of grad amp
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
