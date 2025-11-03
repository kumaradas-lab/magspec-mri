HW.Grad.ExtGradSN = 9;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[1.178049e-05,-0.00272666,-0.001894867,0.004445753]; % Offset Current

HW.Grad.PaUin2PaIout=([0.3356205,0.3322447,0.3330938,0.3380657]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6*5;                                  % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
% HW.Grad.SystemTimeDelay = [49.8e-6, 57e-6, 58.28e-6, 0e-6]; % Time delay of grad amp
HW.Grad.SystemTimeDelay(1:3) = [49.8e-6, 57e-6, 58.28e-6]; % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600
