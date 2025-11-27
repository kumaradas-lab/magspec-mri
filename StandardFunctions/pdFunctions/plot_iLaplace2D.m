function Seq = plot_iLaplace2D(data, Seq)
%% Plot results of inverse Laplace 2D transform
%
%       Seq = plot_iLaplace2D(data, Seq)
%
% This function can be used in combination with get_iLaplace2D to open figures
% showing the results of the inverse Laplace 2D transform.
%
%
% INPUT:
%
%   data
%       Structure with data with the inverse Laplace 2D transform as returned by
%       get_iLaplace2D. This structure must contain at least the field iLaplace
%       which in turn is a structure that must contain at least the following
%       fields:
%
%     SpectrumAmplitude
%         2D matrix with the amplitudes of the iLaplace2D result.
%
%     T1
%         Vector with the T1 times corresponding to the second dimension of
%         SpectrumAmplitude.
%
%     T2
%         Vector with the T2 times corresponding to the first dimension of
%         SpectrumAmplitude.
%
%     DataAmplitude
%         2D matrix with the amplitudes of the input to the iLaplace2D
%         transform.
%
%     tau1
%         Vector with the preparation times corresponding to T1 corresponding to
%         the second dimension of DataAmplitude.
%
%     tau2
%         Vector with the echo times corresponding to T2 corresponding to
%         the first dimension of DataAmplitude.
%
%     FitAmplitude
%         2D matrix with the amplitudes for the fit that corresponds to the
%         iLaplace2D result.
%
%   Seq
%       Structure with settings for the plot. That structure must contain at
%       least the field iLaplace2D which in turn is a structure with the
%       following (optional) fields. If some of those fields are omitted or
%       empty, default values are used:
%
%     plotT1T2Map
%         Scalar number for the figure with the T1-T2 map or a graphics handle
%         that can be a valid parent for the T1-T2 map.
%         (Default: 84)
%
%     plotAmp
%         Scalar number for the figure showing the input amplitudes for the
%         iLaplace2D fit or a graphics handle that can be a valid parent for
%         that plot.
%         (Default: 83)
%
%     plotFitAmp
%         Scalar number for the figure showing the amplitudes for the fit that
%         corresponds to the iLaplace2D result.
%         (Default: 85)
%
%     plotResAmp
%         Scalar number for the figure showing the residuals of the input
%         amplitudes with respect to the fit amplitudes.
%         (Default: 86)
%
%     FullScaleAmplitude
%         The input to the iLaplace2D transform can be normalized by
%         get_iLaplace2D. This value corresponds to the full scale amplitude
%         before the normalization. The value is informational only and is shown
%         in the title of T1-T2 map if it is non-empty.
%         (Default: [])
%
%
% OUTPUT:
%
%   Seq
%       Same as input argument Seq but potentially with updated fields (figure
%       handles).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2024-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input parameters

% plot parents
if isemptyfield(Seq, {'iLaplace2D', 'plotT1T2Map'})
  % figure for T1-T2 map
  Seq.iLaplace2D.plotT1T2Map = 84;
end
if isemptyfield(Seq.iLaplace2D, 'plotAmp')
  % figure for measured amplitude
  Seq.iLaplace2D.plotAmp = 83;
end
if isemptyfield(Seq.iLaplace2D, 'plotFitAmp')
  % figure for fitted amplitude
  Seq.iLaplace2D.plotFitAmp = 85;
end
if isemptyfield(Seq.iLaplace2D, 'plotResAmp')
  % figure for residual amplitude
  Seq.iLaplace2D.plotResAmp = 86;
end

if ~isfield(Seq.iLaplace2D, 'FullScaleAmplitude')
  Seq.iLaplace2D.FullScaleAmplitude = []; %max(abs(data.iLaplace2D.DataAmplitude(:)));
end

tau1 = data.iLaplace2D.tau1;
tau2 = data.iLaplace2D.tau2;


%% plot "measured" signals
if (isnumeric(Seq.iLaplace2D.plotAmp) && Seq.iLaplace2D.plotAmp > 0) ...
    || ishghandle(Seq.iLaplace2D.plotAmp, 'figure') ...
    || ishghandle(Seq.iLaplace2D.plotAmp, 'uipanel')
  if isnumeric(Seq.iLaplace2D.plotAmp)
    hf = figure(Seq.iLaplace2D.plotAmp);
    clf(hf);
  else
    hf = Seq.iLaplace2D.plotAmp;
    delete(get(hf, 'Children'));
  end
  hax = axes(hf);
  surf(hax, tau2, tau1, data.iLaplace2D.DataAmplitude.', 'LineStyle', 'none');
  colormap(hax, 'jet');
  if isempty(Seq.iLaplace2D.FullScaleAmplitude)
    title(hax, 'Measured Signal Amplitude');
  else
    title(hax, ['Measured Signal Amplitude [% FS] (FS = ' num2str(Seq.iLaplace2D.FullScaleAmplitude*1e9) ' nT)']);
  end
  ylabel(hax, 'Preparation Time [s]');
  xlabel(hax, 'Echo Time [s]');
end


%% plot final result
if Seq.iLaplace2D.plotFitAmp > 0
  hf = figure(Seq.iLaplace2D.plotFitAmp);
  clf(hf);
  hax = axes(hf);
  surf(hax, tau2, tau1, data.iLaplace2D.FitAmplitude.', 'LineStyle', 'none');
  colormap(hf, 'jet');
  title(hax, 'Fitted Signal Amplitude');
  ylabel(hax, 'Preparation Time [s]');
  xlabel(hax, 'Echo Time [s]');
end

if Seq.iLaplace2D.plotResAmp > 0
  hf = figure(Seq.iLaplace2D.plotResAmp);
  clf(hf);
  hax = axes(hf);
  surf(hax, tau2, tau1, data.iLaplace2D.DataAmplitude.' - data.iLaplace2D.FitAmplitude.', ...
    'LineStyle', 'none');
  colormap(hf, 'jet');
  title(hax, {'Residuals', 'Measured Signal Amplitude - Fitted Signal Amplitude'});
  ylabel(hax, 'Preparation Time [s]');
  xlabel(hax, 'Echo Time [s]');
end


%% T1-T2-map
if (isnumeric(Seq.iLaplace2D.plotT1T2Map) && Seq.iLaplace2D.plotT1T2Map > 0) ...
    || ishghandle(Seq.iLaplace2D.plotT1T2Map, 'figure') ...
    || ishghandle(Seq.iLaplace2D.plotT1T2Map, 'uipanel')
  % full-scale amplitude of the T1-T2 spectrum (in an integral sense)
  FS = sum(data.iLaplace2D.SpectrumAmplitude(:));

  T1 = data.iLaplace2D.T1;
  T2 = data.iLaplace2D.T2;

  if isnumeric(Seq.iLaplace2D.plotT1T2Map)
    hf = figure(Seq.iLaplace2D.plotT1T2Map);
    set(hf, 'Name', 'T1-T2-map');
  else
    hf = Seq.iLaplace2D.plotT1T2Map;
  end
  ax(1) = subplot(10,10,[3:10,13:20,23:30,33:40,43:50,53:60,63:70,73:80], 'Parent', hf);
  s = data.iLaplace2D.SpectrumAmplitude.';
  s(~conv2(double(s~=0), ones(3,3), 'same')) = NaN; % nans are transparent
  pcolor(ax(1), T2, T1, s * 100/FS);
  set(ax(1), 'XScale', 'log', 'YScale', 'log', 'Tag', 'main');

  shading(ax(1), 'interp');
  cmap = colormap(ax(1), 'jet');
  cmap = [1,1,1; cmap];
  colormap(ax(1), cmap);
  axis(ax(1), [min(T2) max(T2) min(T1) max(T1)]);

  if isempty(Seq.iLaplace2D.FullScaleAmplitude)
    title(ax(1), 'Fitted Spectrum');
  else
    title(ax(1), ['Fitted Spectrum [% FS] (FS = ' num2str(Seq.iLaplace2D.FullScaleAmplitude*1e9) ' nT)']);
  end
  set(ax(1), 'XTickLabel', []);
  set(ax(1), 'YTickLabel', []);
  grid(ax(1), 'on');

  ax(2) = subplot(10,10,2:10:72, 'Parent', hf);
  semilogy(ax(2), sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,2), T1, 'LineWidth', 2);
  grid(ax(2), 'on');
  set(ax(2), 'YTickLabel', [], 'XDir', 'reverse', 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'Tag', 'tau1_spec');

  ax(3) = subplot(10,10,83:90, 'Parent', hf, 'Tag', 'tau2_spec');
  semilogx(ax(3), T2, sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,1), 'LineWidth', 2);
  grid(ax(3), 'on');
  set(ax(3), 'XTickLabel', [], 'YAxisLocation', 'right');

  ax(4) = subplot(10,10,1:10:71, 'Parent', hf);
  semilogy(ax(4), cumsum(sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,2),1), T1, 'LineWidth', 2);
  ylabel(ax(4), 'T1 [s]');
  xlim(ax(4), [0, 100]);
  grid(ax(4), 'on');
  set(ax(4), 'XDir', 'reverse', 'XAxisLocation', 'top', 'XTickLabelRotation', 90, 'Tag', 'tau1_accum');

  ax(5) = subplot(10,10,93:100, 'Parent', hf, 'YAxisLocation', 'right');
  semilogx(ax(5), T2, cumsum(sum(data.iLaplace2D.SpectrumAmplitude.'*100/FS,1),2), 'LineWidth', 2);
  xlabel(ax(5), 'T2 [s]');
  ylim(ax(5), [0, 100]);
  grid(ax(5), 'on');
  set(ax(5), 'YAxisLocation', 'right', 'Tag', 'tau2_accum');

  % Another callback function might have been queued that deletes these axes.
  % Call "drawnow" here (instead of implicitly in "linkaxes" to avoid using
  % references to deleted graphics object.
  drawnow();
  if ~all(ishghandle(ax, 'axes'))
    return;
  end

  linkaxes(ax([1,2,4]), 'y');
  linkaxeskeep(ax([1,3,5]), 'x');
  set(ax(1), 'XLim', [min(T2(:)), max(T2(:))], 'YLim', [min(T1(:)), max(T1(:))]);
end


end
