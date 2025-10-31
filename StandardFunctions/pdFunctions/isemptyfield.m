function is_empty = isemptyfield(structure, field_name)
% ISEMPTYFIELD Check if field in structure is empty
%
%   is_empty = isemptyfield(structure, field_name)
%
% Returns true if structure has no field "field_name" or it is empty.
% Returns false if the field is not empty.
% If "field_name" is a cell, it is checked if the corresponding nested field
% exists and is not empty.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if iscell(field_name)
  for i_cell = 1:numel(field_name)
    is_empty = ~isfield(structure, field_name{i_cell}) || isempty(structure.(field_name{i_cell}));
    if is_empty
      return;
    end
    structure = structure.(field_name{i_cell});
  end
else
  is_empty = ~isfield(structure, field_name) || isempty(structure.(field_name));
end

end

