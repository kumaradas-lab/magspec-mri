function PIDStruct=PIDController(PIDStruct)
% PIDStruct.TargetValue=30;     % Target temperature
% PIDStruct.Actor=1;            % Actors number used [1:4]
% PIDStruct.ControlMax=1;       % Amperes of heating current at 48 V supply
% PIDStruct.Sensor=5;           % Sensor number used [1:8]
% PIDStruct.Kp=0.3;             % Proportional error ganin
% PIDStruct.Ki=0.01;            % Integral error ganin
% PIDStruct.Kd=2;               % Derivative error ganin
% PIDStruct.Plot=1;             % Plot the input and output values
% PIDStruct.nPlotMeasurments=500; % Plot the last n Measurments
% for n=1:1000
%   PIDStruct=PIDController(PIDStruct);
% end
% PIDStruct.exit=1;
% PIDStruct=PIDController(PIDStruct);

if ~isfield(PIDStruct,'USBTD'); PIDStruct.USBTD=[]; end; if isempty(PIDStruct.USBTD); PIDStruct.USBTD=get_USBTemperatureDevice(); end;
if ~isfield(PIDStruct,'exit'); PIDStruct.exit=[]; end; if isempty(PIDStruct.exit); PIDStruct.exit=0; end;

if ~isfield(PIDStruct,'TargetValue'); PIDStruct.TargetValue=[]; end; if isempty(PIDStruct.TargetValue); PIDStruct.TargetValue=0; end;
if ~isfield(PIDStruct,'nSmoothDerivativeError'); PIDStruct.nSmoothDerivativeError=[]; end; if isempty(PIDStruct.nSmoothDerivativeError); PIDStruct.nSmoothDerivativeError=5; end;
if ~isfield(PIDStruct,'nSmoothControlSignal'); PIDStruct.nSmoothControlSignal=[]; end; if isempty(PIDStruct.nSmoothControlSignal); PIDStruct.nSmoothControlSignal=1; end;
if ~isfield(PIDStruct,'Plot'); PIDStruct.Plot=[]; end; if isempty(PIDStruct.Plot); PIDStruct.Plot=1; end;
if ~isfield(PIDStruct,'nPlotMeasurments'); PIDStruct.nPlotMeasurments=[]; end; if isempty(PIDStruct.nPlotMeasurments); if PIDStruct.Plot; PIDStruct.nPlotMeasurments=60;else PIDStruct.nPlotMeasurments=2; end; end;
if ~isfield(PIDStruct,'Actor'); PIDStruct.Actor=[]; end; if isempty(PIDStruct.Actor); PIDStruct.Actor=1; end;
if ~isfield(PIDStruct,'ActorWeight'); PIDStruct.ActorWeight=[]; end; if isempty(PIDStruct.ActorWeight); PIDStruct.ActorWeight=ones(size(PIDStruct.Actor)); end;
if ~isfield(PIDStruct,'ActorMax'); PIDStruct.ActorMax=[]; end; if isempty(PIDStruct.ActorMax); PIDStruct.ActorMax=2.5*ones(size(PIDStruct.Actor)); end;
if ~isfield(PIDStruct,'ActorMin'); PIDStruct.ActorMin=[]; end; if isempty(PIDStruct.ActorMin); PIDStruct.ActorMin=zeros(size(PIDStruct.Actor)); end;
if ~isfield(PIDStruct,'ActorOffValue'); PIDStruct.ActorOffValue=[]; end; if isempty(PIDStruct.ActorOffValue); PIDStruct.ActorOffValue=zeros(size(PIDStruct.Actor)); end;
MyOnCleanup=onCleanup(@()AbortUSBTD(PIDStruct));
if ~isfield(PIDStruct,'Sensor'); PIDStruct.Sensor=[]; end; if isempty(PIDStruct.Sensor); PIDStruct.Sensor=1; end;
if ~isfield(PIDStruct,'SensorWeight'); PIDStruct.SensorWeight=[]; end; if isempty(PIDStruct.SensorWeight); PIDStruct.SensorWeight=ones(size(PIDStruct.Sensor)); end;
if ~isfield(PIDStruct,'Kp'); PIDStruct.Kp=[]; end; if isempty(PIDStruct.Kp); PIDStruct.Kp=1; end;
if ~isfield(PIDStruct,'Ki'); PIDStruct.Ki=[]; end; if isempty(PIDStruct.Ki); PIDStruct.Ki=0; end;
if ~isfield(PIDStruct,'Kd'); PIDStruct.Kd=[]; end; if isempty(PIDStruct.Kd); PIDStruct.Kd=0; end;
if ~isfield(PIDStruct,'Measurments'); PIDStruct.Measurments=[]; end; if isempty(PIDStruct.Measurments); PIDStruct.Measurments=max([PIDStruct.nPlotMeasurments(:);PIDStruct.nSmoothDerivativeError(:);PIDStruct.nSmoothControlSignal(:);2]); end;
if ~isfield(PIDStruct,'TargetValueMeasurments'); PIDStruct.TargetValueMeasurments=[]; end; if isempty(PIDStruct.TargetValueMeasurments); PIDStruct.TargetValueMeasurments=ones(PIDStruct.Measurments,1)*PIDStruct.TargetValue; end;
if ~isfield(PIDStruct,'StartTime'); PIDStruct.StartTime=[]; end; if isempty(PIDStruct.StartTime); PIDStruct.StartTime=now*24*3600; end;
if ~isfield(PIDStruct,'Time'); PIDStruct.Time=[]; end; if isempty(PIDStruct.Time); PIDStruct.Time=nan(PIDStruct.Measurments,1);end;
if ~isfield(PIDStruct,'InputValue'); PIDStruct.InputValue=[]; end; if isempty(PIDStruct.InputValue); PIDStruct.InputValue=nan(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'SensorValues'); PIDStruct.SensorValues=[]; end; if isempty(PIDStruct.SensorValues); PIDStruct.SensorValues=nan(PIDStruct.Measurments,numel(PIDStruct.Sensor)); end;
if ~isfield(PIDStruct,'ActorValue'); PIDStruct.ActorValue=[]; end; if isempty(PIDStruct.ActorValue); PIDStruct.ActorValue=nan(PIDStruct.Measurments,numel(PIDStruct.Actor)); end;
if ~isfield(PIDStruct,'ProportionalError'); PIDStruct.ProportionalError=[]; end; if isempty(PIDStruct.ProportionalError); PIDStruct.ProportionalError=zeros(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'IntegralError'); PIDStruct.IntegralError=[]; end; if isempty(PIDStruct.IntegralError); PIDStruct.IntegralError=zeros(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'DerivativeError'); PIDStruct.DerivativeError=[]; end; if isempty(PIDStruct.DerivativeError); PIDStruct.DerivativeError=zeros(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'SmoothDerivativeError'); PIDStruct.SmoothDerivativeError=[]; end; if isempty(PIDStruct.SmoothDerivativeError); PIDStruct.SmoothDerivativeError=nan(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'ControlSignal'); PIDStruct.ControlSignal=[]; end; if isempty(PIDStruct.ControlSignal); PIDStruct.ControlSignal=nan(PIDStruct.Measurments,1); end;
if ~isfield(PIDStruct,'SmoothControlSignal'); PIDStruct.SmoothControlSignal=[]; end; if isempty(PIDStruct.SmoothControlSignal); PIDStruct.SmoothControlSignal=nan(PIDStruct.Measurments,1); end;

if ~isfield(PIDStruct,'ControlMin'); PIDStruct.ControlMin=[]; end; if isempty(PIDStruct.ControlMin); PIDStruct.ControlMin=0; end;
if ~isfield(PIDStruct,'ControlMax'); PIDStruct.ControlMax=[]; end; if isempty(PIDStruct.ControlMax); PIDStruct.ControlMax=4; end;
if ~isfield(PIDStruct,'ControlSpan'); PIDStruct.ControlSpan=[]; end; if isempty(PIDStruct.ControlSpan); PIDStruct.ControlSpan=PIDStruct.ControlMax-PIDStruct.ControlMin; end;
if ~isfield(PIDStruct,'ControlMean'); PIDStruct.ControlMean=[]; end; if isempty(PIDStruct.ControlMean); PIDStruct.ControlMean=PIDStruct.ControlMin+(PIDStruct.ControlMax-PIDStruct.ControlMin)/2; end;
if ~isfield(PIDStruct,'ControlMinProportional'); PIDStruct.ControlMinProportional=[]; end; if isempty(PIDStruct.ControlMinProportional); PIDStruct.ControlMinProportional=PIDStruct.ControlMean-PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'ControlMaxProportional'); PIDStruct.ControlMaxProportional=[]; end; if isempty(PIDStruct.ControlMaxProportional); PIDStruct.ControlMaxProportional=PIDStruct.ControlMean+PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'ControlMinIntegral'); PIDStruct.ControlMinIntegral=[]; end; if isempty(PIDStruct.ControlMinIntegral); PIDStruct.ControlMinIntegral=PIDStruct.ControlMean-PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'ControlMaxIntegral'); PIDStruct.ControlMaxIntegral=[]; end; if isempty(PIDStruct.ControlMaxIntegral); PIDStruct.ControlMaxIntegral=PIDStruct.ControlMean+PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'ControlMinDerivative'); PIDStruct.ControlMinDerivative=[]; end; if isempty(PIDStruct.ControlMinDerivative); PIDStruct.ControlMinDerivative=PIDStruct.ControlMean-PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'ControlMaxDerivative'); PIDStruct.ControlMaxDerivative=[]; end; if isempty(PIDStruct.ControlMaxDerivative); PIDStruct.ControlMaxDerivative=PIDStruct.ControlMean+PIDStruct.ControlSpan*1.5; end;
if ~isfield(PIDStruct,'hAxesInput'); PIDStruct.hAxesInput=[]; end;
if ~isfield(PIDStruct,'hAxesOutput'); PIDStruct.hAxesOutput=[]; end;

if PIDStruct.exit==0;
    if PIDStruct.Plot
        if ~isfield(PIDStruct,'hFigure'); PIDStruct.hFigure=[]; end; 
        if isempty(PIDStruct.hFigure); 
            PIDStruct.hFigure=figure(23); 
            PIDStruct.hAxesInput=subplot(2,1,1);
            PIDStruct.hAxesOutput=subplot(2,1,2);
            linkaxes([PIDStruct.hAxesInput,PIDStruct.hAxesOutput],'x')
        end;
    end

    PIDStruct.TargetValueMeasurments=[PIDStruct.TargetValueMeasurments(2:end);PIDStruct.TargetValue];

    PIDStruct.SensorValues(1:end-1,:)=PIDStruct.SensorValues(2:end,:);
    for n=1:numel(PIDStruct.Sensor)
        PIDStruct.SensorValues(end,n)=PIDStruct.USBTD.ADT7320_GetTemperature(PIDStruct.Sensor(n));
    end
    PIDStruct.NormalizedSensorWeight=abs(PIDStruct.SensorWeight(:).')/sum(abs(PIDStruct.SensorWeight(:)));
    PIDStruct.InputValue=[PIDStruct.InputValue(2:end);sum(PIDStruct.SensorValues(end,:).*PIDStruct.NormalizedSensorWeight(:).',2)];
    PIDStruct.Time=[PIDStruct.Time(2:end);now*24*3600-PIDStruct.StartTime];

    PIDStruct.ProportionalError=[PIDStruct.ProportionalError(2:end);PIDStruct.TargetValueMeasurments(end)-PIDStruct.InputValue(end)];
    if and((PIDStruct.Time(end)-PIDStruct.Time(end-1))>0,sum(~isnan(PIDStruct.Time(:)))>1);
        PIDStruct.IntegralError=[PIDStruct.IntegralError(2:end);PIDStruct.IntegralError(end)+(PIDStruct.ProportionalError(end))*(PIDStruct.Time(end)-PIDStruct.Time(end-1))];
        PIDStruct.DerivativeError=[PIDStruct.DerivativeError(2:end);((PIDStruct.TargetValueMeasurments(end)-PIDStruct.InputValue(end))-(PIDStruct.TargetValueMeasurments(end-1)-PIDStruct.InputValue(end-1)))/(PIDStruct.Time(end)-PIDStruct.Time(end-1))];
    else
        PIDStruct.IntegralError=[PIDStruct.IntegralError(2:end);0];
        PIDStruct.DerivativeError=[PIDStruct.DerivativeError(2:end);0];
    end
    PIDStruct.SmoothDerivativeError=[PIDStruct.SmoothDerivativeError(2:end);mean(PIDStruct.DerivativeError(end-min(PIDStruct.nSmoothDerivativeError,sum(~isnan(PIDStruct.Time(:))))+1:end))];
    PIDStruct.ProportionalError(end)=max(min(PIDStruct.ProportionalError(end),PIDStruct.ControlMaxProportional/PIDStruct.Kp),PIDStruct.ControlMinProportional/PIDStruct.Kp);
    PIDStruct.IntegralError(end)=max(min(PIDStruct.IntegralError(end),PIDStruct.ControlMaxIntegral/PIDStruct.Ki),PIDStruct.ControlMinIntegral/PIDStruct.Ki);
    PIDStruct.SmoothDerivativeError(end)=max(min(PIDStruct.SmoothDerivativeError(end),PIDStruct.ControlMaxDerivative/PIDStruct.Kd),PIDStruct.ControlMinDerivative/PIDStruct.Kd);
    PIDStruct.ControlSignal=[PIDStruct.ControlSignal(2:end);PIDStruct.ProportionalError(end)*PIDStruct.Kp+PIDStruct.IntegralError(end)*PIDStruct.Ki+PIDStruct.SmoothDerivativeError(end)*PIDStruct.Kd];
    PIDStruct.ControlSignal(end)=max(min(PIDStruct.ControlSignal(end),PIDStruct.ControlMax),PIDStruct.ControlMin);
    PIDStruct.SmoothControlSignal=[PIDStruct.SmoothControlSignal(2:end);mean(PIDStruct.ControlSignal(end-min(PIDStruct.nSmoothControlSignal,sum(~isnan(PIDStruct.Time(:))))+1:end))];

    PIDStruct.NormalizedActorWeight=abs(PIDStruct.ActorWeight)/sum(abs(PIDStruct.ActorWeight(:)));
    PIDStruct.ActorValue(1:end-1,:)=PIDStruct.ActorValue(2:end,:);
    for n=1:numel(PIDStruct.Actor)
        PIDStruct.ActorValue(end,n)=max(min(PIDStruct.SmoothControlSignal(end)*PIDStruct.NormalizedActorWeight(n),PIDStruct.ActorMax(n)),PIDStruct.ActorMin(n));
        PIDStruct.USBTD.DAC8565_SetVoltage(PIDStruct.Actor(n),PIDStruct.ActorValue(end,n));
    end

    if ~isempty(PIDStruct.hAxesInput)
        plot(PIDStruct.hAxesInput,  PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.InputValue(end-PIDStruct.nPlotMeasurments+1:end),...
                                    PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.TargetValueMeasurments(end-PIDStruct.nPlotMeasurments+1:end));
        legend(PIDStruct.hAxesInput,{'InputValue','TargetValue'},'Location','West','Box', 'off','Color', 'none')
    end
    if ~isempty(PIDStruct.hAxesOutput)
        plot(PIDStruct.hAxesOutput,     PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.SmoothControlSignal(end-PIDStruct.nPlotMeasurments+1:end),...
                                        PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.ProportionalError(end-PIDStruct.nPlotMeasurments+1:end)*PIDStruct.Kp,...
                                        PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.IntegralError(end-PIDStruct.nPlotMeasurments+1:end)*PIDStruct.Ki,...
                                        PIDStruct.Time(end-PIDStruct.nPlotMeasurments+1:end),PIDStruct.SmoothDerivativeError(end-PIDStruct.nPlotMeasurments+1:end)*PIDStruct.Kd);
        legend(PIDStruct.hAxesOutput,{'ControlSignal','ProportionalError*Kp','IntegralError*Ki','DerivativeError*Kd'},'Location','West','Box', 'off','Color', 'none')
    end
    if or(~isempty(PIDStruct.hAxesInput),~isempty(PIDStruct.hAxesOutput))
        drawnow
    end
else
    for n=1:numel(PIDStruct.Actor)
        PIDStruct.USBTD.DAC8565_SetVoltage(PIDStruct.Actor(n),PIDStruct.ActorOffValue(n));
    end
end
end

function AbortUSBTD(PIDStruct)
    for nn=1:numel(PIDStruct.Actor)
        PIDStruct.USBTD.DAC8565_SetVoltage(PIDStruct.Actor(nn),PIDStruct.ActorOffValue(nn));
    end
end


