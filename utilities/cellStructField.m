function fieldValues = cellStructField(cellStruct, field)
% ieldValues = cellStructField(cellStruct, field)
%
% cellStructField concatenates field values of individual structures within
% a cell array of such structures.
%
% Args:
%   cellStruct (cell, required, positional): a shape-(M, N) cell array of
%     matlab structures.
%   field (char, required, positional): a shape-(1, L) character array
%     specifying the structure field of interest.
%
% Returns:
%   fieldValues (cell): a shape-(M*N,1) cell array of individual field
%     values.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  cellStruct (:,:) {mustBeA(cellStruct,'cell')}
  field (1,:) {mustBeA(field,'char')}
end

nCells = numel(cellStruct);
fieldValues = cell(nCells, 1);
for iCell = 1:nCells
  fieldValues{iCell} = cellStruct{iCell}.(field);
end