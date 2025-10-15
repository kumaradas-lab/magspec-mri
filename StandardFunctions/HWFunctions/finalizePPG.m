function sPPG = finalizePPG(HW, sPPG, config)
%% Add commands to PPG as returned by set_sequence to allow it to be executed
%
%   sPPG = finalizePPG(HW, sPPG, config)
%
% INPUT:
%
%   HW
%       HW object (or structure).
%
%   sPPG
%       .NET pulse program as returned in one of the fields of Seq.PPG by
%       set_sequence.
%
%   config
%       Structure with the following optional fields:
%
%     newCLTimeFirst
%       Update the first CLTime in the pulse program to this value. (Default: 0)
%
%     GetDataAtEnd
%       Boolean value to indicate whether package end should be appended to the
%       pulse sequence. (Default: false)
%
% OUTPUT:
%
%   sPPG
%       The input .NET pulse program with additional commands that allow it to
%       be executed by the console.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.PPGWrapper
%#function hex2uint64
%#function isemptyfield

