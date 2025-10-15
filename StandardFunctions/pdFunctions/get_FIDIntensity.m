function Intensity = get_FIDIntensity(Intensity)

if isemptyfield(Intensity,'FIDAbsCicTc');Intensity.FIDAbsCicTc=Intensity.FIDAbsCic;end
if isemptyfield(Intensity,'FIDTime'); warning('get_FIDIntensity:FIDTime','Intensity.FIDTime is empty'); Intensity.FIDTime=(1:size(Intensity.FIDAbsCicTc))';end
if isemptyfield(Intensity,'tFIDStart');Intensity.tFIDStart=0;end
if isemptyfield(Intensity,'tFIDEnd');Intensity.tFIDEnd=Intensity.FIDTime(end);end
if isemptyfield(Intensity,'nPolyfit');Intensity.nPolyfit=1;end
if isemptyfield(Intensity,'PolyfitPlot');Intensity.PolyfitPlot=30;end


iS=find(Intensity.tFIDStart<Intensity.FIDTime(:,1),1,'first');
iE=find(Intensity.FIDTime(:,1)<=Intensity.tFIDEnd,1,'last');
nIMs=Intensity.nInterleave;

if Intensity.PolyfitPlot
    figure(Intensity.PolyfitPlot);
    clf(Intensity.PolyfitPlot,'reset');
    set(Intensity.PolyfitPlot,'DefaultAxesColorOrder',hsv(size(Intensity.FIDAbsCicTc,2)))
    ax1=subplot(3,1,1,'Parent',Intensity.PolyfitPlot);
    ax2=subplot(3,1,2,'Parent',Intensity.PolyfitPlot);
    ax3=subplot(3,1,3,'Parent',Intensity.PolyfitPlot);
    plot(ax1,Intensity.FIDTime(iS:iE,1:nIMs:end),Intensity.FIDAbsCicTc(iS:iE,1:nIMs:end).*1e12)
end

if (Intensity.nPolyfit-1) || Intensity.PolyfitPlot
   Intensity.FIDIntensity=0;
   TempTime=Intensity.FIDTime(iS:iE,1);
   TempTimeOut=[linspace(0,TempTime(1),10)';TempTime];
   TempTime=[-TempTime(end:-1:1);TempTime];
   Intensity.polyfitFID=nan(numel(TempTimeOut),size(Intensity.FIDAbsCicTc,2));
   Intensity.polyfitFIDDiff=nan(size(Intensity.FIDAbsCicTc(iS:iE,:)));
   for tt=1:size(Intensity.FIDAbsCicTc,2)
       TempData1=Intensity.FIDAbsCicTc(iS:iE,tt);
       TempData=[TempData1(end:-1:1);TempData1];

       [p,S,mu] = polyfit(TempTime,TempData,Intensity.nPolyfit);
       Intensity.FIDIntensity(1,tt) = polyval(p,0,S,mu); 
       Intensity.polyfitFID(:,tt) = polyval(p,TempTimeOut,S,mu); 
       Intensity.polyfitFIDDiff(:,tt) = polyval(p,Intensity.FIDTime(iS:iE,1),S,mu)-Intensity.FIDAbsCicTc(iS:iE,tt); 
   end
   if Intensity.PolyfitPlot
       hold(ax1, 'all');
       plot(ax1,zeros(size(Intensity.polyfitFID(1,1:nIMs:end))),Intensity.FIDIntensity(:,1:nIMs:end).*1e12,'x')
       plot(ax1,repmat(TempTimeOut,1,size(Intensity.polyfitFID(:,1:nIMs:end),2)),Intensity.polyfitFID(:,1:nIMs:end).*1e12)
       hold(ax1, 'off');
       xlim(ax1,[0,Intensity.tFIDEnd]);
%        ylim(ax,[0,inf])
       xlabel(ax1,'time / s')
       ylabel(ax1,'Intensity / pT')
       grid(ax1,'on')
       plot(ax2,Intensity.FIDTime(iS:iE,1:nIMs:end),Intensity.polyfitFIDDiff(:,1:nIMs:end)./mean(Intensity.FIDIntensity(:,1:nIMs:end))*100)
       grid(ax2,'on')
       xlabel(ax2,'time / s')
       ylabel(ax2,'Intensity fit diff / %')
       linkaxes([ax1,ax2],'x');
       plot(ax3,(Intensity.FIDIntensity(:,1:nIMs:end)-mean(Intensity.FIDIntensity(:,1:nIMs:end)))./mean(Intensity.FIDIntensity(:,1:nIMs:end))*100,'-x')
       grid(ax3,'on')
       xlabel(ax3,'time / s')
       ylabel(ax3,'Intensity @ t=0 diff to mean / %')
   end

else
   Intensity.FIDIntensity=mean(Intensity.FIDAbsCic(iS:iE,1:1:end));
end
end
