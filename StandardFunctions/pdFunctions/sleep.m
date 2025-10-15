function  sleep(seconds)
% sleep(seconds)
if 1
    endTime=now*24*3600+seconds;
    while (endTime-now*24*3600)>0.01    
        java.lang.Thread.sleep(10);            % Note: sleep() accepts [mSecs] duration
    end 
    lastTime=(endTime-now*24*3600);
    if lastTime>0;java.lang.Thread.sleep(lastTime*1000);end
else
    pause(seconds);
end
end