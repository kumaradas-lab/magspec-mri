classdef IgnoreWarning < handle
  %% Disable warnings and reset on deletion
  %
  %   protectWarning = PD.IgnoreWarning(ignoredWarnIds)
  %
  % This class can be used to reset the warning state and the values that are
  % returned by the "lastwarn" function after a section of code has run. During
  % that section of code, the given warnings are ignored.
  % For this, create an object of this class at the start of the section of code
  % and delete it after the end of that section.
  %
  % INPUT:
  %
  %   ignoredWarnIds
  %     String or cell array of strings with warning ids which are ignored while
  %     this object is valid. If the warning id returned by "lastwarn()" when
  %     this object is deleted matches any of these warning ids, the last
  %     warning message and the last warning id are reset to the values that
  %     they had when this object was created. Note that intermediate return
  %     values for lastwarn from other (not ignored) warnings might have been
  %     overridden by ignored warnings at that point. So the information about
  %     not ignored warnings might be unrecoverable when this object is deleted.
  %
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------

  properties (GetAccess = private, SetAccess = private)
    ignoredWarningIds

    lastWarnMsg
    lastWarnId
    warnState
  end


  methods

    function this = IgnoreWarning(ignoredWarnIds)
      %% save current warning state and disable warnings

      if ~iscell(ignoredWarnIds)
        ignoredWarnIds = {ignoredWarnIds};
      end

      this.ignoredWarningIds = ignoredWarnIds;

      [this.lastWarnMsg, this.lastWarnId] = lastwarn('', '');

      this.warnState = warning();
      for iWarn = 1:numel(ignoredWarnIds)
        warning('off', ignoredWarnIds{iWarn});
      end
    end


    function delete(this)
      %% reset warning state and reset lastwarn if appropriate

      warning(this.warnState);

      % If the new last warning id matches any of the disabled ones, reset to
      % the old last warning.
      [~, newLastWarnId] = lastwarn();
      if any(strcmp(this.ignoredWarningIds, newLastWarnId))
        lastwarn(this.lastWarnMsg, this.lastWarnId);
      end
    end

  end

end
