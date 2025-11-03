function [pulseData] = Pulse_Rect_x2Amplitude(HW, Center, BW, FlipAngle, maxpulseCount,  maxLength, Frequency, Phase)
% create a rectangular RF pulse
% x=2 * 6/5/0.86;
% x=2 * 12/9/0.86;
x=2 * 12/9/0.86;
% x=2*6/4/0.86;
x=1/5*HW.tFlip180Def/100e-6/1;
x=1/5*HW.tFlip180Def/100e-6/1*1.16;
x=1/5*HW.tFlip180Def/100e-6/1;
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

    BlockLength=x/BW*0.999;

    gain=HW.TX.AmpDef*tFlipPi*(FlipAngle/pi)/(BlockLength/0.998);

    if maxLength<BlockLength;
        error('maxLength of HF Pulse to short')
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


%% ------------------------------------------------------------------------
% (C) Copyright 2019-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%------------------------------------------------------------------------
