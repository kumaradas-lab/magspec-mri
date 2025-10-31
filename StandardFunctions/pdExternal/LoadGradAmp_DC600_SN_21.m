HW.Grad.ExtGradSN = 21;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.002355968,0.001336827,0.0004297138,-0.001263458]; % Offset Current 13-May-2019 10:28:16
 
HW.Grad.PaUin2PaIout=([0.3304762,0.3344283,0.3327338,0.3316342]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 13-May-2019 10:28:16
 
HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [1.9584e-05   1.9648e-05   2.2192e-05]; % Time delay of grad amp mit Magnet #135
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600

