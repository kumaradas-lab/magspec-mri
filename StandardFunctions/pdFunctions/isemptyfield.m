function is_empty = isemptyfield(structure, field_name)
% Check if field in structure is empty
%
%   is_empty = isemptyfield(structure, field_name)
%
% Returns true if structure has no field "field_name" or it is empty.
% Returns false if the field is not empty.
% If "field_name" is a cell, it is checked if the corresponding nested field
% exists and is not empty. In case any of the nested fields are an array of
% structures, it is assumed that all of those structures are equivalent to
% the first element up to the highest checked nesting level.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if iscell(field_name)
  for i_cell = 1:numel(field_name)
    if numel(structure) > 1
      % We need to do that conditionally because some classes don't allow
      % indexing with ().
      % FIXME: This assumes that all elements in the array of struct are
      % equivalent.
      structure = structure(1);
    end
    is_empty =  isempty(structure) ...
      || ~isfield(structure, field_name{i_cell}) ...
      || isempty(structure.(field_name{i_cell}));
    if is_empty
      return;
    end
    structure = structure.(field_name{i_cell});
  end
else
  if numel(structure) > 1
    % We need to do that conditionally because some classes don't allow
    % indexing with ().
    structure = structure(1);
  end
  is_empty = isempty(structure) ...
    || ~isfield(structure, field_name) ...
    || isempty(structure.(field_name));
end

end

