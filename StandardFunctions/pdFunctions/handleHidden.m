classdef handleHidden < handle
    %% HANDLEHIDDEN Derived from class handle with hidden default methods.
    % The hidden methods are still accessible and can be overloaded.
    % However, they do not show or auto-complete.
    %
    % handleHidden methods:
    %   addlistener  - Add listener for event.
    %   delete       - Delete a handle object.
    %   eq           - Test handle equality.
    %   findobj      - Find objects with specified property values.
    %   findprop     - Find property of MATLAB handle object.
    %   ge           - Greater than or equal relation.
    %   gt           - Greater than relation.
    %   isvalid      - Test handle validity. (Is sealed and cannot be
    %                  hidden)
    %   le           - Less than or equal relation for handles.
    %   lt           - Less than relation for handles.
    %   ne           - Not equal relation for handles.
    %   notify       - Notify listeners of event.
    %
    % See also:
    % handle
    %
    
    methods (Hidden)
        %% These hidden functions are to hide those inherited from handle
        % handle methods:
        %   addlistener  - add listener for event
        %   delete       - delete a handle object.
        %   eq           - test handle equality.
        %   findobj      - find objects with specified property values.
        %   findprop     - Find property of MATLAB handle object.
        %   ge           - Greater than or equal relation.
        %   gt           - Greater than relation.
        %   isvalid      - Test handle validity. (Is sealed and cannot be
        %                  hidden)
        %   le           - Less than or equal relation for handles.
        %   lt           - Less than relation for handles.
        %   ne           - Not equal relation for handles.
        %   notify       - Notify listeners of event.
        %
        
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function delete(varargin)
          % destructor of parent class is called automatically
        end
        function result = eq(varargin)
            result = eq@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function result = ge(varargin)
            result = ge@handle(varargin{:});
        end
        function result = gt(varargin)
            result = gt@handle(varargin{:});
        end
        function result = le(varargin)
            result = le@handle(varargin{:});
        end
        function result = lt(varargin)
            result = lt@handle(varargin{:});
        end
        function result = ne(varargin)
            result = ne@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
    end
    
end
