function textureizeslice(slice, onoff)
% Convert a regular slice into a texture map slice, or a texture
% slice into a regular slice.

for k = 1:numel(slice)

  % d = getappdata(slice(k), 'textureoptimizeations');

  switch onoff
    case 'on'
      % d.xdata = get(slice(k), 'XData');
      % d.ydata = get(slice(k), 'YData');
      % d.zdata = get(slice(k), 'ZData');
      % setappdata(slice(k), 'textureoptimizeations', d);
      % if max(size(d.xdata)==1)
      %   nx = [d.xdata(1) d.xdata(end)];
      % else
      %   nx = [d.xdata(1,1)   d.xdata(1,end);
      %         d.xdata(end,1) d.xdata(end,end)];
      % end
      % if max(size(d.ydata)==1)
      %   ny = [d.ydata(1) d.ydata(end)];
      % else
      %   ny = [d.ydata(1,1)   d.ydata(1,end);
      %         d.ydata(end,1) d.ydata(end,end)];
      % end
      % if max(size(d.zdata)==1)
      %   nz = [d.zdata(1) d.zdata(end)];
      % else
      %   nz = [d.zdata(1,1)   d.zdata(1,end);
      %         d.zdata(end,1) d.zdata(end,end)];
      % end
      % set(slice(k), 'XData', nx, 'YData', ny, 'ZData', nz, ...
      %   'FaceColor', 'texturemap');
      if ischar(get(slice(k), 'FaceAlpha'))
        set(slice(k), 'FaceAlpha', 'texturemap');
      end
      if ischar(get(slice(k), 'FaceColor'))
        set(slice(k), 'FaceColor', 'texturemap');
      end
    case 'off'
      % if ~isempty(d)
      %   set(slice(k), 'XData', d.xdata, 'YData', d.ydata, 'ZData', d.zdata);
      %   setappdata(slice(k), 'textureoptimizeations', []);
      % end
      if ischar(get(slice(k), 'FaceAlpha')) && strcmp(get(slice(k), 'FaceAlpha'), 'texturemap')
        set(slice(k), 'FaceAlpha', 'flat');
      end
      if ischar(get(slice(k), 'FaceColor')) && strcmp(get(slice(k), 'FaceColor'), 'texturemap')
        set(slice(k), 'FaceColor', 'flat');
      end
  end
end

end
