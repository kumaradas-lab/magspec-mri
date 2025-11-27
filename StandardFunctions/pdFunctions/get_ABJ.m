function [ABJ]=get_ABJ(ABJ)
if nargin == 0
  ABJ = 'COM10';
end
if nargin <= 1
  if ~isstruct(ABJ)
    temp=ABJ;
    SerObj = instrfind({'Port'},{temp});
    if length(SerObj)>=1
      fclose(SerObj(1));
      delete(SerObj(1));
      clear SerObj(1)
    end
    clear ABJ
    ABJ.COM=temp;
    ABJ.command='D09';
    ABJ.getWeight=0;
    ABJ.String='';
    ABJ.weight=zeros(0,1);
    ABJ.clearOldWeight=1;
    ABJ.weightTimeClear=now*24*3600;
    ABJ.weightTimeMin=[];
    ABJ.weightTimeMax=[];
    ABJ.weightTime=[];

  end
  if ~isfield(ABJ,'command')
    ABJ.command = '';
  end
  % if isempty(ABJ.command)
  %   ABJ.command = '';
  % end
  if ~isfield(ABJ,'getWeight')
    ABJ.getWeight = 0;
  end
  if ~isfield(ABJ, 'waitForWeight')
    ABJ.waitForWeight = 0;
  end
  if ~isfield(ABJ, 'clearOldWeight')
    ABJ.clearOldWeight = 0;
  end
  if ~isfield(ABJ,'Timeout')
    ABJ.Timeout = 1;
  end
end
if ~isfield(ABJ,'serial')
  ABJ.serial = [];
end
if isempty(ABJ.serial)
    ABJ.serial = serial(ABJ.COM,'BaudRate',1200,'DataBits',8,'StopBits',1,'Parity','none','Timeout',ABJ.Timeout,'Terminator','CR'); % ,'FlowControl','hardware'
    fopen(ABJ.serial);
end
if isemptyfield(ABJ, 'plot')
  ABJ.plot = 0;
end
if isemptyfield(ABJ, 'close')
  ABJ.close = 0;
end
% if isemptyfield(ABJ.status)
%   ABJ.status = 0;
% end

if ABJ.getWeight && ~ABJ.clearOldWeight
  if get(ABJ.serial, 'Timeout') ~= ABJ.Timeout
    set(ABJ.serial, 'Timeout', ABJ.Timeout);
  end
  t = numel(ABJ.weight);
  while get(ABJ.serial, 'BytesAvailable')
    t = t+1;
    try
      tempStr = fscanf(ABJ.serial);
    catch
      tempStr = '   NAN ';
      warning('ABJ: RS232 timout')
    end
    ABJ.String = strcat( ABJ.String,tempStr);
    tempStrVal = tempStr([1:8,10]); %DF.1
    % tempStrVal=tempStr(4:14); %DF.2
    % tempStrVal=tempStr(1:10); %DF.3
    % tempStrVal=tempStr(2:12); %DF.4
    ABJ.weight(t)=str2double(strrep(tempStrVal,' ',''))/1000;
    ABJ.weightTimeMax(t)=now*24*3600;
    if numel(ABJ.weightTimeMin)<numel(ABJ.weightTimeMax)
      ABJ.weightTimeMin(t)=ABJ.weightTimeMax(t);
    end
    ABJ.weightTime(t)=(ABJ.weightTimeMax(t)+ABJ.weightTimeMin(t))/2;
    sleep(0.2);
  end
elseif ABJ.clearOldWeight
  ABJ.String='';
  ABJ.weight=[];
  ABJ.weightTimeClear=now*24*3600;
  ABJ.weightTimeMin=[];
  ABJ.weightTimeMax=[];
  ABJ.weightTime=[];
  while get(ABJ.serial,'BytesAvailable')
    tt=get(ABJ.serial,'Timeout');
    set(ABJ.serial,'Timeout',0.1);
    fscanf(ABJ.serial);
    set(ABJ.serial,'Timeout',tt);
  end
end

switch ABJ.command
  case 'D01' % Fortlaufende Datenausgabe
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D02' % Fortlaufende Datenausgabe stabiler Wägewerte
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D03' % Status der Stabilitätsanzeige wird bei der fortlaufende Ausga-be den Daten angehängt. U: instabil S: stabil
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D05' % Einmalige Ausgabe
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D06' % Automatische Ausgabe
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D07' % Einmalige Ausgabe. Status der Stabilitätsanzeige wird bei der Ausgabe den Daten angehängt.
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D08' % Einmalige Ausgabe bei stabilem Wägewert
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  case 'D09' % Ausgabe abbrechen
    fprintf(ABJ.serial,ABJ.command);
    ABJ.waitForWeight=0;
  case 'BREAK' % BREAK Taste
    fprintf(ABJ.serial,ABJ.command);
  case 'CAL' % CAL Taste
    fprintf(ABJ.serial,ABJ.command);
  case 'TARE' % TARE Taste
    fprintf(ABJ.serial,ABJ.command);
  case 'PRINT' % PRINT Taste
    fprintf(ABJ.serial,ABJ.command);
    ABJ.weightTimeMin(numel(ABJ.weight)+1)=now*24*3600;
  otherwise

end

if ABJ.waitForWeight
  set(ABJ.serial,'Timeout',ABJ.Timeout);
  tempStr=fscanf(ABJ.serial);
  t=numel(ABJ.weight);
  while get(ABJ.serial,'BytesAvailable')||~isempty(tempStr)
    t=t+1;
    if isempty(tempStr)
      tempStr=fscanf(ABJ.serial);
      if isempty(tempStr);
        tempStr='      NAN          ';
      end
    end
    ABJ.String=strcat( ABJ.String,tempStr);
    tempStrVal=tempStr([1:8,10]); %DF.1
    % tempStrVal=tempStr(4:14); %DF.2
    % tempStrVal=tempStr(1:10); %DF.3
    % tempStrVal=tempStr(2:12); %DF.4
    ABJ.weight(t)=str2double(strrep(tempStrVal,' ',''))/1000;
    ABJ.weightTimeMax(t)=now*24*3600;
    ABJ.weightTime(t)=(ABJ.weightTimeMax(t)+ABJ.weightTimeMin(t))/2;
    tempStr='';
  end
end

if ABJ.plot
  disp([ 'The weight is ' num2str(ABJ.weight) ' g']);
end

if ABJ.close
  fclose(ABJ.serial);
  delete(ABJ.serial);
  clear ABJ.serial;
end

%% clear all serial
if 0
%%
SerObj=instrfind
while length(SerObj)>=1
  fclose(SerObj(1))
  delete(SerObj(1))
  clear SerObj(1)
  SerObj=instrfind;
end
clear ABJ
%%
end
