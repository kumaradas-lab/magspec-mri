function [axisMinOut, axisMaxOut] = get_axisMinMax(axisMin, axisMax, axisMinOld, axisMaxOld)
%%
% axisMin=-0.001;
% axisMax=+0.011;
% axisMinOld=0.001;
% axisMaxOld=0.003;
axisSpan=axisMax-axisMin;
expo=floor(log10(axisSpan));
% axisSpanCeil=ceil(axisSpan*10^-expo)*10^expo;
axisMinOut=floor(axisMin*10^-expo)*10^expo;
axisMaxOut=ceil(axisMax*10^-expo)*10^expo;
if nargin == 4
  if axisMin*10^-expo < axisMinOld*10^-expo
    axisMinOut=axisMin;
  elseif axisMin*10^-expo <= axisMinOld*10^-expo +1
    axisMinOut=axisMinOld;
  else
    axisMinOut=axisMin;
  end

  if axisMax*10^-expo > axisMaxOld*10^-expo
    axisMaxOut=axisMax;
  elseif axisMax*10^-expo >= axisMaxOld*10^-expo -1
    axisMaxOut=axisMaxOld;
  else
    axisMaxOut=axisMax;
  end
end

end
