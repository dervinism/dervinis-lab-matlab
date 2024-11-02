function dataStruct = pkl2mat(pklFilename, options)
% [dataStruct, matFilename] = pkl2mat(pklFilename, <filename2save>)
%
% Converts data in python PKL format to a Matlab structure.
%
% Args:
%   pklFilename (char, required, positional): a shape-(1, N) character
%     array with full data file path.
%   filename2save (char, optional, keyword): a shape-(1, L) character
%     array with the MAT file name to save the converted data. If left
%     empty, converted data will not be saved.
%
% Returns:
%   dataStruct (struct): a Matlab structure containing the converted data
%     loaded from the PKL file.
%
% Comments:
%   Conversion may not be full. Therefore, examine the function output
%   structure and perform  any further conversion manually.
%
% Dependencies:
%   Python
%   (https://uk.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
%   numpy
%   pandas
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  pklFilename (1,:) {mustBeVector,mustBeText,endsWith(pklFilename,'.pkl')}
  options.filename2save (1,:) {mustBeVector,mustBeText,endsWith(options.filename2save,'.mat')}
end

% Import the pickle module
pickle = py.importlib.import_module('pickle');

% Open the pickle file
fid = py.open(pklFilename, 'rb');

% Load the data from the pickle file
dataStruct = pickle.load(fid);

% Close the file
fid.close();

% Convert the Python data container into a Matlab structure
dataStruct = py2mat(dataStruct);

% Convert fields of the data structure
if isstruct(dataStruct) % First data layer
  fieldsL1 = fieldnames(dataStruct);
  if ~isempty(fieldsL1)
    for iFieldL1 = 1:numel(fieldsL1)
      dataStruct.(fieldsL1{iFieldL1}) = py2mat(dataStruct.(fieldsL1{iFieldL1}));

      if iscell(dataStruct.(fieldsL1{iFieldL1})) && ~isempty(dataStruct.(fieldsL1{iFieldL1})) % Second data layer
        for iCell = 1:numel(dataStruct.(fieldsL1{iFieldL1}))
          dataStruct.(fieldsL1{iFieldL1}){iCell} = py2mat(dataStruct.(fieldsL1{iFieldL1}){iCell});

          if isstruct(dataStruct.(fieldsL1{iFieldL1}){iCell}) % Third data layer
            fieldsL2 = fieldnames(dataStruct.(fieldsL1{iFieldL1}){iCell});
            if ~isempty(fieldsL2)
              for iFieldL2 = 1:numel(fieldsL2)
                dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2}) = ...
                  py2mat(dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2}));

                if isstruct(dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2})) % Fourth data layer
                  fieldsL3 = fieldnames(dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2}));
                  if ~isempty(fieldsL3)
                    for iFieldL3 = 1:numel(fieldsL3)
                      dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2}).(fieldsL3{iFieldL3}) = ...
                        py2mat(dataStruct.(fieldsL1{iFieldL1}){iCell}.(fieldsL2{iFieldL2}).(fieldsL3{iFieldL3}));
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end

% Save the converted data
if ~isempty(options.filename2save)
  save(options.filename2save, 'dataStruct', '-v7.3');
end