function matlabStruct = py2mat(pythonObj)
% matlabStruct = py2mat(pythonStruct)
%
% Converts a Python data object (a list, a dictionary, a numpy array, etc.)
% into a matlab structure.
%
% Args:
%   pythonObj (obj, required, positional): a shape-(1, 1) Python object.
%     Supported Python objects:
%     py.list
%     py.dict
%     py.numpy.ndarray
%     py.int
%     py.str
%     py.pandas.core.frame.DataFrame
%     py.pandas.core.series.Series
%
% Returns:
%   matlabStruct (struct): a Matlab structure containing the converted data
%     from the Python object.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  pythonObj
end

% Convert the data container into a Matlab structure
if isa(pythonObj, 'py.list') % list
  matlabStruct = cell(pythonObj);
elseif isa(pythonObj, 'py.dict') % dictionary
  matlabStruct = struct(pythonObj);
elseif isa(pythonObj, 'py.numpy.ndarray') % numpy array
  try
    matlabStruct = double(pythonObj);
  catch
    matlabStruct = string(cell(pythonObj.tolist));
  end
elseif isa(pythonObj, 'py.int') % integer
  matlabStruct = double(pythonObj);
elseif isa(pythonObj, 'py.str') % string
  matlabStruct = string(pythonObj);
elseif isa(pythonObj, 'py.pandas.core.frame.DataFrame') % data frame
  matlabStruct = table(pythonObj);
elseif isa(pythonObj, 'py.pandas.core.series.Series') % data series
  matlabStruct = table2array(table(pythonObj));
else
  matlabStruct = pythonObj;
  warning('Not an acceptable Python object type. The supplied data type has not been converted.');
end