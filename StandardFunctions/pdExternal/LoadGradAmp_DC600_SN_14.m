HW.Grad.ExtGradSN = 14;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.003034678,-0.0006314275,0.001174667,0.001662802]; % Offset Current 13-Oct-2017 15:05:51

HW.Grad.PaUin2PaIout=([0.3358092,0.332277,0.3327882,0.333347]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 13-Oct-2017 15:21:00

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
%HW.Grad.SystemTimeDelay(1:3) =[4.7776e-05, 5.6128e-05, 5.644e-05, 5.55e-05]; % Time delay of grad amp
HW.Grad.SystemTimeDelay(1:3) = [4.7776e-05, 5.6128e-05, 5.644e-05]; % Time delay of grad amp
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600
