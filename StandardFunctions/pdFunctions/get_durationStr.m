function durationStr = get_durationStr(durationInSec)
  if abs(durationInSec)>=60
    tt=durationInSec;
    dur=zeros(1,4);
    dur(1)=fix(tt/3600);
    tt=tt-dur(1)*3600;
    dur(2)=fix(tt/60);
    tt=tt-dur(2)*60;
    dur(3)=fix(tt/1);
    tt=tt-dur(3)*1;
    dur(4)=tt*1000;
    if dur(1)
      FormatStr='hh:mm:ss';
      dur(3)=dur(3)+round(dur(4)/1000);
    else  dur(2);FormatStr='mm:ss';
      if  dur(4);FormatStr=[FormatStr '.S'];end
      dur(3)=dur(3)+round(dur(4)/1000*10)/10;
    end
    dur(4)=[];
    durationStr=[char(duration(dur,'Format',FormatStr)) ' s'];
  elseif abs(durationInSec)>=1
    durationStr=[num2str(fix(durationInSec*1e0 ),'%d') '.' num2str(abs(round((durationInSec-fix(durationInSec*1e0)/1e0)*1e3   )),'%03.3d') ' s'];
  elseif abs(durationInSec)>=1e-3
    durationStr=[num2str(fix(durationInSec*1e3 ),'%d') '.' num2str(abs(round((durationInSec-fix(durationInSec*1e3)/1e3)*1e6   )),'%03.3d') ' ms'];
  elseif abs(durationInSec)>=1e-6
    durationStr=[num2str(fix(durationInSec*1e6 ),'%d') '.' num2str(abs(round((durationInSec-fix(durationInSec*1e6)/1e6)*1e9   )),'%03.3d') ' ' char(181) 's'];
  elseif abs(durationInSec)>=1e-9
    durationStr=[num2str(fix(durationInSec*1e9 ),'%d') '.' num2str(abs(round((durationInSec-fix(durationInSec*1e9)/1e9)*1e12  )),'%03.3d') ' ns'];
  else    
    durationStr=[num2str(fix(durationInSec*1e12),'%d') '.' num2str(abs(round((durationInSec-fix(durationInSec*1e12)/1e12)*1e15)),'%03.3d') ' ps'];
  end
end