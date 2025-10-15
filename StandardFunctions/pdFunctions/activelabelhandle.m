function activelabelhandle(ah, label, string)
% activelabelhandle(ah, LABEL, STRING) - Create a label on ah which is
%     active.  LABEL is the property of ah whose label you are
%     setting.  STRING is the initial text string for the label.

  l = get(ah,label);
  
  set(l,'string',string);
  set(l,'buttondownfcn',@activelabelbuttondown);
  
function activelabelbuttondown(obj, action)
% Callback when one of our active labels is clicked on.
  
  set(obj,'edit','on');
  
   