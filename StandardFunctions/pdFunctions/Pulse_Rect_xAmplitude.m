function [pulseData] = Pulse_Rect_xAmplitude(HW, Center, BW, FlipAngle, maxpulseCount,  maxLength, Frequency, Phase)
%% create a rectangular RF pulse
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% x=2 * 6/5/0.86;
% x=2 * 12/9/0.86;
x=2 * 12/9/0.86;
% x=2*6/4/0.86;
x=1/10*HW.tFlip180Def/100e-6/1;
x=1/10*HW.tFlip180Def/100e-6/1*1.03;
x=1/10*HW.tFlip180Def/100e-6/1*1.03/(4100/3000);
x=1/10*HW.tFlip180Def/100e-6/1*1.03/(4100/2500);
x=1/10*HW.tFlip180Def/100e-6/1/(7/4);

x=1/5*HW.tFlip180Def/100e-6/1/2;

if nargin==2
    if ischar(Center)
        if strcmp(Center,'Amp')
            pulseData=1/x;
        elseif strcmp(Center,'Time')
            pulseData=x;
        else
            pulseData=nan;
    end
    else
        Pulse.MaxLength
        Pulse.Bandwidth
        Pulse.Center
        Pulse.FlipAngle
        Pulse.Frequency
        Pulse.Phase
        Pulse.MaxNumberOfSegments

    end

else

    tFlipPi=HW.TX.Amp2FlipPiIn1Sec/HW.TX.AmpDef;

    BlockLength=x/BW;

    gain=HW.TX.AmpDef*tFlipPi*(FlipAngle/pi)/BlockLength;

    if maxLength + 1/HW.TX(Pulse.iDevice).fSample < BlockLength
      error('PD:Pulse_Rect_xAmplitude:MaxLengthTooShort', ...
        'MaxLength of rf pulse is %.3f %cs too short.', ...
        (BlockLength - maxLength)*1e6, char(181));
    end

    if maxpulseCount<1
        error('maxpulseCount >= 1');
    end

    pulseData.Start=-BlockLength/2+Center;
    pulseData.Amplitude=gain;
    pulseData.Duration=BlockLength;
    pulseData.Frequency=Frequency;
    pulseData.Phase=Phase;

end

end
