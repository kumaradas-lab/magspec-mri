function HW_obj = load_GradAmp(HW_obj, gradAmpPropFile, channels)
%% Load only specified channels from gradient amplifier properties file
%
%     HW = load_GradAmp(HW, gradAmpPropFile, channels)
%
% Only properties in HW.Grad will be loaded.
%
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   gradAmpPropFile
%       Name of the script that contains the properties of the gradient
%       amplifier.
%
%   channels
%       Indices of channels to be loaded. If elements that are defined in the
%       properties file contain four (4) elements, only the elements specified
%       in this vector will be assigned in the output HW. (Default: 1:4)
%
%
% OUTPUT:
%
%   HW
%       HW object or structure where elements in HW.Grad are updated with the
%       values defined in the gradient amplifier properties script.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% FIXME: Fields might need to be indexed differently if they are in xyzB
%        ordering.


%% Check input
if nargin < 2
  channels = 1:4;
end

if ~exist(gradAmpPropFile, 'file')
  error('PD:load_GradAmp:NoFile', ...
    'The file "%s" could not be found.', gradAmpPropFile);
end


%% Execute script to get properties in HW structure
eval(gradAmpPropFile);


%% recursively assign fields from HW structure to output object
for iDevice = 1:numel(HW.Grad)
  HW_obj.Grad(iDevice) = ...
    assign_fields(HW_obj.Grad(iDevice), HW.Grad(iDevice), channels);
end

end


function obj = assign_fields(obj, stru, channels)
%% Assign fields from struct to properties in object

fields = fieldnames(stru);

for iField = 1:numel(fields)
  if isstruct(stru.(fields{iField}))
    % recursively call function for sub-structures
    obj.(fields{iField}) = ...
      assign_fields(obj.(fields{iField}), stru.(fields{iField}), channels);
  elseif numel(stru.(fields{iField})) == 4
    if (isstruct(obj) && isemptyfield(obj, fields{iField})) || ...
        (isobject(obj) && isempty(obj.(fields{iField})))
      % Not yet set. So, we aren't overriding anything.
      obj.(fields{iField}) = stru.(fields{iField});
    else
      % override only selected channels
      obj.(fields{iField})(channels) = stru.(fields{iField})(channels);
    end
  else
    % override elements in field
    num_elem = numel(stru.(fields{iField}));
    obj.(fields{iField})(1:num_elem) = stru.(fields{iField});
  end
end

end
