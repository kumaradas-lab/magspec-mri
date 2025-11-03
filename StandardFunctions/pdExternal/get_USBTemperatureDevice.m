function USBTemperatureDevice = get_USBTemperatureDevice()
% NET.addAssembly(strcat(mySave.reg.LibPath,'\PDTools.dll'));

USBTemperatureDevice = PDTools.USBTemperature();
if USBTemperatureDevice.isAvailable
  USBTemperatureDevice.Init();
  for t = 1:8
    USBTemperatureDevice.ADT7320_GetTemperature(t);
  end
  USBTemperatureDevice.ADS1220_GetTemperature();
  for t = 1:4
    USBTemperatureDevice.DAC8565_SetVoltage(t,0.0);
  end
end

end
