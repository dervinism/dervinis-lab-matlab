function [reducedNumericArray, reducedNumericArrayIdx] = selectArrayValues(numericArray, cutoffs, options)
% reducedNumericArray = selectArrayValues(numericArray, cutoffs, <cutoffType>)
%
% function selects numeric array values within given cutoffs. It also
% provides indices of these values.
%
% Args:
%   numericArray (numeric, required, positional): a shape-(K, L) numeric
%     array of values. Typically, the input would be a vector of spike
%     times.
%   cutoffs (numeric, required, positional): a shape-(L, 2) numeric array
%     of cutoffs for selecting array values. L represents the number of
%     cutoff intervals, while the second dimension corresponds to the lower
%     and the upper cutoff values. Only elements falling within these
%     cutoff values will remain. The rest of the values will be removed.
%   cutoffType (char, optional, keyword): a shape-(1, N) character array
%     describing the type of cutoffs to be used:
%     'inclusive' - includes the cutoff end values (default);
%     'exclusive' - excludes the cutoff end values.
%
% Returns:
%   reducedNumericArray (numeric): a shape-(1, M) numeric array of values
%     (e.g., spike times). This vector contains only selected values from
%     numericArray after applying the cutoffs.
%   reducedNumericArrayIdx (numeric): a shape-(1, M) numeric array of
%     linear indices to numericArray producing reducedNumericArray
%     (reducedNumericArray = numericArray(reducedNumericArrayIdx)).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  numericArray (:,:) {mustBeNumeric,mustBeNonempty}
  cutoffs (:,2) {mustBeNumeric}
  options.cutoffType (1,:) {mustBeMember(options.cutoffType,{'inclusive','exclusive'})} = 'inclusive'
end

% Find indices within given cutoffs
if ~isempty(cutoffs)
  nCutoffs = size(cutoffs,1);
  reducedNumericArrayIdx = [];
  for cutoff = 1:nCutoffs
    if strcmpi(options.cutoffType, 'inclusive')
      cutoffIdx = find(numericArray >= cutoffs(cutoff,1) ...
        & numericArray <= cutoffs(cutoff,2));
    elseif strcmpi(options.cutoffType, 'exclusive')
      cutoffIdx = find(numericArray > cutoffs(cutoff,1) ...
        & numericArray < cutoffs(cutoff,2));
    end
    cutoffIdx = cutoffIdx(:)';
    reducedNumericArrayIdx = [reducedNumericArrayIdx cutoffIdx]; %#ok<*AGROW>
  end
  reducedNumericArrayIdx = sort(unique(reducedNumericArrayIdx));
else
  reducedNumericArrayIdx = 1:numel(numericArray);
end

% Select values within given cutoffs
reducedNumericArray = numericArray(reducedNumericArrayIdx);