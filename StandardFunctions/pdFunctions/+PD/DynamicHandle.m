classdef DynamicHandle < dynamicprops
  %% Handle class that behaves like a structure
  %
  %   PD.DynamicHandle(inStruct)
  %
  % The class itself has no properties. All properties are dynamically added to
  % the object on the fly.
  % Nested struct references in strings are automatically expanded, i.e. it is
  % possible to index an object of this class with "a.('b.c')".
  %
  % INPUT:
  %   inStruct    Optional input structure. All fields with their respective
  %               values are added to the returned object. If inStruct is an
  %               object of class "PD.DynamicHandle", that exact reference is
  %               returned (no copy constructor).
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------


  % Programming note:
  % Using the nested field syntax doesn't work for functions inside this class.
  % Call subsref or subsasgn directly instead.


  properties
    % no "static" properties
  end


  methods

    function this = DynamicHandle(inStruct)
      %% constructor
      if nargin < 1, return; end
      if isa(inStruct, 'PD.DynamicHandle')
        % Do not copy
        this = inStruct;
      elseif isa(inStruct, 'struct')
        fnames = fieldnames(inStruct);
        for iField = 1:numel(fnames)
          this.SetProp(fnames{iField}, inStruct.(fnames{iField}));
        end
      else
        error('PD:DynamicHandle:Convert', ...
          'Cannot convert %s to PD.DynamicHandle', class(inStruct));
      end
    end


    function SetProp(this, prop, value)
      %% dynamically add properties if not yet existent
      if ~isprop(this, prop)
        this.addprop(prop);
      end
      this.(prop) = value;
    end


    function this = subsasgn(this, s, varargin)
      %% overload subsasgn to get "structure-like" syntax
      if ~strcmp(s(1).type, '.')
        error('Class %s must be indexed with ".".', class(this));
      end
      % replace subs containing '.' by nested subs
      s = this.NestifySubs(s);
      % add property if necessary and assign values
      if isscalar(s)
        this.SetProp(s.subs, varargin{:});
      else
        if isfield(this, s(1).subs)
          temp = this.(s(1).subs);
        else
          temp = [];
        end
        % create (nested) structure by calling built-in subsasgn
        temp = subsasgn(temp, s(2:end), varargin{:});
        this.SetProp(s(1).subs, temp);
      end
    end


    function varargout = subsref(this, s)
      %% overload subsref to expand nested structs in string indexing
      if ~strcmp(s(1).type, '.')
        error('Class %s must be indexed with ".".', class(this));
      end
      % replace subs containing '.' by nested subs
      s = this.NestifySubs(s);
      % get value of (nested struct in) property by calling built-in subsref
      [varargout{1:nargout}] = builtin('subsref', this, s);
    end


    function dp = addprop(this, propName)
      %% Add a dynamic property to the class
      dp = addprop@dynamicprops(this, propName);
    end


    function val = isprop(this, propName)
      %% Check whether nested property exists
      subs = regexp(propName, '([^\.]*)', 'match');
      if isempty(subs)
        val = ~isempty(findprop(this, propName));
        return
      end
      if isempty(findprop(this, subs{1}))
        val = false;
      else
        if numel(subs) < 2
          val = true;
        else
          val = isfield_recursive(this.(subs{1}), subs(2:end));
        end
      end
    end


    function delprop(this, propName)
      %% Delete the dynamic property or nested field "propName"
      subs = regexp(propName, '([^\.]*)', 'match');
      dp = findprop(this, subs{1});
      if ~isempty(dp)
        nSubs = numel(subs);
        if nSubs < 2
          delete(dp);
        else
          [s(1:nSubs-1).subs] = deal(subs{1:end-1});
          [s(:).type] = deal('.');
          this.subsasgn(s, rmfield(this.subsref(s), subs{end}));
        end
      else
        warning('PD:DynamicHandle:NoProp', 'DynamicHandle: No such property "%s" (for "%s").', subs{1}, propName);
      end
    end


    function Absorb(this, dynHdl)
      %% Absorb the (nested) fields of second DynamicHandle dynHdl
      nestedFields = dynHdl.GetNestedFields();
      for iNF = 1:numel(nestedFields)
        s.type = '.';
        s.subs = nestedFields{iNF};
        this.subsasgn(s, dynHdl.subsref(s));
      end
    end


    function allProps = GetNestedFields(this)
      %% Get a cell with all nested fields in object
      allProps = properties(this);
      iProp = 1;
      while iProp <= numel(allProps)
        s.type = '.';
        s.subs = allProps{iProp};
        if isstruct(this.subsref(s))
          % add fields of nested struct to end of list
          allProps = [allProps; strcat([allProps{iProp}, '.'], fieldnames(this.subsref(s)))];
          % and remove nested struct from list
          allProps(iProp) = [];
          iProp = iProp - 1;
        end
        iProp = iProp + 1;
      end
    end


    function strct = ToStruct(this)
      %% return the properties of this object in a struct
      protectWarning = PD.IgnoreWarning('MATLAB:structOnObject');
      strct = struct(this);
      delete(protectWarning);
    end


  end


  methods (Hidden)
    function val = isfield(this, fieldName)
      %% Wrapper of "isprop"

      % warning('PD:DynamicHandle:isfield', ...
      %   '"isfield" was used where you should use "isprop". Consider replacing it.');
      val = isprop(this, fieldName);
    end


    function this = rmfield(this, fieldName)
      %% Wrapper of "delprop"
      % ATTENTION: Unlike "rmfield" on a struct, this method does NOT create a
      % copy of the object.

      % Issue a warning to track down where "rmfield" is inappropriately used in
      % the code.
      % warning('PD:DynamicHandle:rmfield', ...
      %   '"rmfield" was used where you should use "delprop". Consider replacing it.')
      this.delprop(fieldName);
    end
  end


  methods (Static)

    function s = NestifySubs(sIn)
      %% replace subs containing '.' by nested subs
      subsTypeDot = strcmp({sIn(:).type}, '.');
      subsContainingDot = cellfun(@any, strfind({sIn(subsTypeDot).subs}, '.'));
      if any(subsContainingDot)
        s = struct('subs', {}, 'type', {});
        for iSubs = 1:numel(sIn)
          if subsTypeDot(iSubs)
            subs = regexp(sIn(iSubs).subs, '[^\.]*', 'match');
            s(end+(1:numel(subs))) = struct('subs', subs(:), 'type', '.');
          else
            s(end+1) = sIn(iSubs);
          end
        end
      else
        s = sIn;
      end
    end

  end

end


function val = isfield_recursive(strct, fieldCell)
%% recursively check fields in struct

val = false;
if isfield(strct, fieldCell{1})
  if numel(fieldCell) > 1
    val = isfield_recursive(strct.(fieldCell{1}), fieldCell(2:end));
  else
    val = true;
  end
end

end
