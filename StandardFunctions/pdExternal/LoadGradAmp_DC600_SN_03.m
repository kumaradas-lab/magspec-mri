HW.Grad.ExtGradSN = 3;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset Voltage
HW.Grad.PaOffsetI=[0.003021,-0.0040142,-0.0012837,0.00013454]; % Offset Current

HW.Grad.PaUin2PaUout=[0,0,0,0];                         % Voltage Gain
HW.Grad.PaUin2PaIout=([0.33677,0.33099,0.33327,0.33236]-HW.Grad.PaOffsetI)./1; % Input voltage to output current ratio

HW.Grad.PaPmaxInt=[117,117,117,117];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = 19.5e-6;                 % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600
