function res = pdNamedMMF()  %#ok<STOUT>
%% Interface to named memory mapped files
%
%   res = pdNamedMMF(file_name, id, max_size)
%   res = pdNamedMMF(field_name, data)
%   res = pdNamedMMF(field_name)
%   res = pdNamedMMF('release')
%   res = pdNamedMMF('unlock')
%
%
% This function is used to instantiate and interface to a named memory mapped
% file.
% Currently, only one single memory mapped file can be used at the same time.
% Currently, this function is only implemented for Windows.
% A memory map is created which can be accessed using a C struct with the
% following fields:
%     int tag;
%     int SR;
%     int dl;
%     int info;
%     double data[max_size];
%
%
% INPUT:
%
%   The behaviour of the function dependends on the number of input arguments.
%   With three input arguments a named memory mapped file can be instantiated.
%   In this case, the arguments are as follows:
%
%     file_name
%       The name of the existing(!) file on disc for which a memory map will be
%       instantiated.
%
%     id
%       A (unique) identifier that can be used by other programs to access the
%       memory mapped file.
%
%     max_size
%       The (maximum) size of the memory map in bytes. The existing file must be
%       large enough to hold a structure with that size. I.e.:
%       3*4 (for the 3 int fields) + 4 (for alignment of the subsequent double
%       array) + max_size*8 (for the double array)
%
%   With two input arguments, the value of a field in the map can be set.
%   In this case, the arguments are as follows:
%
%     field_name
%       Name of the field to be set in the structure. This must be 'tag', 'SR',
%       'data', or 'info. (The value of 'dl' is set automatically matching the
%       number of elements in data.)
%
%     data
%       New value for selected field. This must be a numeric scalar of type
%       `int132` or `double` for the fields 'tag' or 'SR'; or it must be a
%       double vector not larger than `max_size` for 'data'.
%
%   With one input argument, the currently set value of a field in the map can
%   be querried or "special" tasks can be triggered. In this case, the argument
%   is as follows:
%
%     field_name
%       Name of the field in the structure to be querried. This must be 'tag',
%       'SR', 'data', or 'info'.
%
%     'release'
%       Release the memory associated to the memory mapped file. After this, the
%       behavior of other programs accessing the memory mapped file will be
%       undefined.
%
%     'unlock'
%       Unlock the memory that hold the Matlab function containing handles
%       associated to the memory mapped file. Without this, the corresponding
%       .mex* file cannot be modified. After this, accessing the memory mapped
%       file will result in undefined behavior.
%
%
% OUTPUT:
%
%   res
%     For the two input argument syntax, this returns the value currently set in
%     the structure corresponding to the memory mapped file. For the other
%     syntaxes, the function returns 1 on success.
%
%
% For more details on memory mapped files on Windows, see:
% https://learn.microsoft.com/en-us/windows/win32/memory/using-the-memory-management-functions
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com

end
