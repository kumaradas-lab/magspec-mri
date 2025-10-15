function structure = set_EmptyField(structure, fieldName, defaultValue)
%% Set field to default value if it is empty
%
%   structure = set_EmptyField(structure, fieldName, defaultValue)
%
% Potentially adds the field "fieldName" with the value "defaultValue" to the
% structure and returns it. If the field "fieldName already exists, the function
% returns the structure as is.
%
% See also: isemptyfield
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2017 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------    

if ~isfield(structure, fieldName) || isempty(structure.(fieldName))
  structure.(fieldName) = defaultValue;
end

end

