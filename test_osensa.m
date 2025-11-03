clear all;
close all;

file_name='test'; % don't forget to change the name before next experiment

transmitter=enable_osensa("COM3"); % change COM port if needed
ts=1; %time step in seconds [s]
n=30; % total recording time ~ 30 s
data=zeros(n,2); % two columns: time, temperature
start=tic;
for i=1:n
    tic
    data(i,1)=toc(start);
    data(i,2)=transmitter.read_channel_temp();
    plot(data(1:n,1),data(1:n,2),'r*')
    ylim([20 37])
    pause(ts-toc)
end
save([file_name '.mat'], 'data')
csvwrite(file_name,data)
saveas(gcf,[file_name '.png'])
transmitter.close(); disp("Osensa Transmitter OFF")


