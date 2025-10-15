classdef UnwindProtectGuard < handle
  % Execute function handle passed to constructor in destructor if
  % property "isUnwindProtect" is true.
  % Very similar to Matlab's onCleanup class

  properties
    isUnwindProtect = true;
    fcnHandle = @defaultFcn;
  end

  methods
    % Contructor
    function this = UnwindProtectGuard(fHandle)
      if ~isa(fHandle, 'function_handle')
        error('Only argument must be a function handle to a function with no arguments')
      end

      this.fcnHandle = fHandle;
    end

    % Destructor
    function delete(this)
      if this.isUnwindProtect
        this.fcnHandle();
      end
    end
  end
end

function defaultFcn()
% do nothing
end
