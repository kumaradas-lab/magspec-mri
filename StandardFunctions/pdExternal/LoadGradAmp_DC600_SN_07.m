HW.Grad.ExtGradSN = 7;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled, set 1; if voltage controlled, set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.001748187,0.0007014147,0.006626038,-0.00250315]; % Offset Current 03-Aug-2017 15:01:12 @ 42 deg C

HW.Grad.PaUin2PaIout=([0.3332023,0.3349715,0.3401493,0.3309973]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio 03-Aug-2017 15:03:37 @ 42 deg C

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Output impedance

HW.Grad.tRamp = 50e-6;                                  % minimum ramp time
HW.Grad.tEC = 50e-6;                                    % Eddy current Setting time;
HW.Grad.SystemTimeDelay(1:3) = [4.1456e-05, 5.1472e-05, 5.172e-05]; % Time delay of grad amp
HW.Grad.MaxAmpSlice = 0.1;                              % maximum gradient amplitude for slice selection

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600
