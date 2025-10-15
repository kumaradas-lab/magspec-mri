HW.Grad.ExtGradSN=1001;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[0,0,0,0];                  % If current controlled set 1, if Voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                      % Offset voltage
HW.Grad.PaOffsetI=[0,0,0,0];    % Offset current 

HW.Grad.PaUin2PaUout=[1,1,1,1];                         % Voltage gain
HW.Grad.PaUin2PaIout=([0,0,0,0]-HW.Grad.PaOffsetI)./1; % Input voltage to output current ratio
HW.Grad.PaPmaxInt=[100,100,100,100];                      % Maximum internal power dissipation

HW.Grad.PaRout=[0,0,0,0];                 % Ouput impedance

HW.Grad.tRamp=18e-6*4;                                  % minimum ramp time
HW.Grad.tEC=36e-6*2;                                    % Setting time;
% HW.Grad.SystemTimeDelay = [4.8208e-05, 5.6208e-05, 5.452e-05, 1.95e-05]; % Time delay of grad amp
HW.Grad.SystemTimeDelay(1:3) = 00e-6;                   % Time delay of grad amp
% HW.Grad.SystemTimeDelay(1:3) =[4.1904e-05   6.1216e-05   6.1536e-05]; % Time delay of grad amp
% HW.Grad.SystemTimeDelay(1:3) =[9.6288e-05  0.000112736   0.00011904]; % Time delay of grad amp
HW.Grad.SystemTimeDelay(1:3) =[4.6296e-05   5.6568e-05   5.9056e-05]; % Time delay of grad amp

HW.Grad.Status1=0;                                      % Power supply ok of DC600
HW.Grad.Status2=0;                                      % Gradient and Temperatur ok of DC600

HW.Grad.CoilMaxDcCurrent=[0.75,0.75,0.75,0.75];         % A
HW.Grad.CoilCurrentSquareTime=[0.9,0.9,0.9,0.9];        % A^2*sec

