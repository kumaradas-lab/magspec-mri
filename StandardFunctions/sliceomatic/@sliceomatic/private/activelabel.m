function hText = activelabel(hax, label, string)
% ACTIVELABEL(HAX, LABEL, STRING) - Create a label on HAX which is
%     active.  LABEL is the property of GCA whose label you are
%     setting.  STRING is the initial text string for the label.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

hText = get(hax, label);

set(hText, 'String', string);
set(hText, 'ButtonDownFcn', @ActiveLabelButtonDown_Callback);

end

function ActiveLabelButtonDown_Callback(hText, eventData)
% Callback when one of the active labels is clicked on.

set(hText, 'Editing', 'on');

end
