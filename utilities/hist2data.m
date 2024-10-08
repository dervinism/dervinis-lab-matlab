function data = hist2data(counts, bins)
% Converts a histogram into data points.
%
% Args:
%   counts (numeric, required, positional): a shape-(1, N) numeric array of
%     histogram counts.
%   bins (numeric, required, positional): a shape-(1, N) numeric array of
%     histogram bin values (centres, not edges).
%
% Returns:
%   data (numeric): a shape-(1, sum(counts)) numeric array of data points
%     generated from the supplied histogram.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

data = [];
for i = 1:length(counts)
  binData = repmat(bins(i), 1, counts(i));
  data = [data, binData]; %#ok<AGROW>
end