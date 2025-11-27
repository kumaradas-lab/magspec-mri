%% Settings for LNA #51 24.3 MHz

% Low power broadband LNA without switch built May 2023

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 51;

HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

% Gain settings:
% Use H1 gain for gain settings due to the fact that the H1 frequency is present
% in all double resonant coils. The signal gain of the specific nucleus needs to
% be calibrated as a ratio of signal strength in relation of H1. Gain
% values other than H1 are only for evaluation perposes.

HW.RX(iDevice).LNAGain = 10^(((-48.9241)-(-71.2525))/20);  % 22.3284 dB gain @ 24300449.2558 MHz F=0.74565 dB (-71.2525 dBm cal); H1 frequency
% HW.RX(iDevice).LNAGain = 10^(((-48.7684)-(-71.1181))/20);  % 22.3497 dB gain @ 22900423.4079 MHz F=0.97399 dB (-71.1181 dBm cal); F19 frequency
% HW.RX(iDevice).LNAGain = 10^(((-47.0248)-(-69.9535))/20);  % 22.9287 dB gain @ 9853182.6597 MHz F=1.9575 dB (-69.9535 dBm cal); P31 frequency
% HW.RX(iDevice).LNAGain = 10^(((-46.5966)-(-69.6535))/20);  % 23.0569 dB gain @ 6122113.3418 MHz F=2.1475 dB (-69.6535 dBm cal); C13 frequency
% HW.RX(iDevice).LNAGain = 10^(((-47.2423)-(-69.4435))/20);  % 22.2012 dB gain @ 2466045.5994 MHz F=2.5866 dB (-69.4435 dBm cal); N15 frequency

% The following gain values are also only for evaluation perposes.
% They show the noise level of the console which was used for the
% calibration.
% HW.RX(iDevice).LNAGain = 10^(((-71.2555)-(-71.2525))/20);  % -0.0029737 dB gain @ 24300449.0963 MHz F=3.6881 dB (-71.2525 dBm cal); H1 frequency
% HW.RX(iDevice).LNAGain = 10^(((-71.1199)-(-71.1181))/20);  % -0.0017713 dB gain @ 22900423.3288 MHz F=3.6224 dB (-71.1181 dBm cal); F19 frequency
% HW.RX(iDevice).LNAGain = 10^(((-69.956)-(-69.9535))/20);  % -0.002519 dB gain @ 9853182.6956 MHz F=4.16 dB (-69.9535 dBm cal); P31 frequency
% HW.RX(iDevice).LNAGain = 10^(((-69.6471)-(-69.6535))/20);  % 0.0064081 dB gain @ 6122113.4102 MHz F=4.4162 dB (-69.6535 dBm cal); C13 frequency
% HW.RX(iDevice).LNAGain = 10^(((-69.4443)-(-69.4435))/20);  % -0.00082715 dB gain @ 2466045.6111 MHz F=4.253 dB (-69.4435 dBm cal); N15 frequency
