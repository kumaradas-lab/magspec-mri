function [fCoil] = get_fCoil(HW)
Seq.plotSmith = 0;
Seq.plotNyquist = 0;
 Seq.plotReflection=1;
Seq.plotReflectionRaw=0;
Seq.Cal=0;

Seq.fCenter=HW.fLarmor;
Seq.fSpan=1e6;
Seq.fSteps=101;

[Network,~]=sequence_Network(HW, Seq);
if(20*log10(abs(Network.MinReflection)) < -6)
    HW.fLarmor = Network.fMinReflection;
    
    Seq.fCenter=HW.fLarmor;
    Seq.fSpan=0.05e6;
    Seq.fSteps=301;
    [Network,~]=sequence_Network(HW, Seq);
    if(20*log10(abs(Network.MinReflection)) < -6)
%         fCoil = Network.fMinReflection + Seq.fSpan/Seq.fSteps*rand(1);
        try        
            myP=polyfit((Network.Frequency-Network.fMinReflection)./(diff(Network.Frequency([1,end]))),abs(Network.Reflection)-abs(Network.MinReflection),12);
            fCoil = fminbnd(@(x) polyval(myP,x),-0.1,0.1)*diff(Network.Frequency([1,end]))+Network.fMinReflection;
        catch
            fCoil = Network.fMinReflection + Seq.fSpan/Seq.fSteps*(rand(1)-0.5);
        end
    else
        warning('Frequency not found.');
        fCoil = HW.fLarmor;
    end
else
    warning('Frequency not found.');
    fCoil = HW.fLarmor;
end
% fCoil = HW.fLarmor;
end

