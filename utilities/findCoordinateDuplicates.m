function dupGroups = findCoordinateDuplicates(names, xyz, options)
% dupGroups = findCoordinateDuplicates(names, xyz, <options>)
%
% Function finds items that share the same 3D coordinate. It is a
% lightweight quality check for electrode-localization files and electrode
% tables: two distinct contacts at an identical coordinate almost always
% means a localization mistake (the same point clicked twice, or an
% unlocalised contact that defaulted to a neighbour).
%
% Args:
%   names (string | cellstr | char, required, positional): a shape-(n, 1)
%     list of item names, one per coordinate row.
%   xyz (numeric, required, positional): a shape-(n, 3) numeric array of
%     [x, y, z] coordinates. Rows with any non-finite value (NaN/Inf) are
%     ignored.
%   Tolerance (numeric, optional, keyword): a shape-(1, 1) non-negative
%     scalar (default 0). With 0, coordinates must match bit-exactly. With
%     a positive value, coordinates are snapped to a grid of that spacing
%     before comparison, so points closer than roughly the tolerance on
%     every axis are grouped (a cheap approximation of clustering, adequate
%     for a QC warning).
%
% Returns:
%   dupGroups (cell): a shape-(m, 1) cell array; each element is a struct
%     with fields:
%       names - a string array of the >= 2 names sharing the coordinate
%       xyz   - the shared [x, y, z] (the first member's coordinate)
%     Returned as a 0x1 cell when there are no duplicates.
%
% Dependencies:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  names {mustBeVector}
  xyz (:,3) {mustBeNumeric}
  options.Tolerance (1,1) {mustBeNumeric, mustBeNonnegative} = 0
end

names = string(names);
names = names(:);
if numel(names) ~= size(xyz, 1)
  error('findCoordinateDuplicates:sizeMismatch', ...
    'names (%d) and xyz (%d rows) must have the same length.', ...
    numel(names), size(xyz, 1));
end

dupGroups = cell(0, 1);

finiteRows = find(all(isfinite(xyz), 2));
if numel(finiteRows) < 2
  return
end

P = xyz(finiteRows, :);
if options.Tolerance > 0
  P = round(P ./ options.Tolerance) .* options.Tolerance;
end

[~, ~, groupIndex] = unique(P, 'rows', 'stable');
for iGroup = 1:max(groupIndex)
  members = finiteRows(groupIndex == iGroup);
  if numel(members) > 1
    s = struct();
    s.names = names(members);
    s.xyz = xyz(members(1), :);
    dupGroups{end+1, 1} = s; %#ok<AGROW>
  end
end
