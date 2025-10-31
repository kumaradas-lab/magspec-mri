HW.Grad.ExtGradSN = 2;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            %Offset voltage
HW.Grad.PaOffsetI=[-0.0008233,-0.0014797,-0.00027688,-0.0006372]; % Offset Current

HW.Grad.PaUin2PaIout=([0.33243,0.33313,0.33282,0.33314]-HW.Grad.PaOffsetI)./1; % Input voltage to output current ratio

HW.Grad.PaPmaxInt=[90,100,120,120];                     % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = 70e-6;                 % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600
