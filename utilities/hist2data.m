function data = hist2data(counts, bins)
% Converts a histogram into data points.
%
% Args:
%   counts (numeric, required, positional): a shape-(M, N) numeric array of
%     histogram counts. Matrix rows correspond to individual histograms.
%   bins (numeric, required, positional): a shape-(1, N) numeric array of
%     histogram bin values (centres, not edges).
%
% Returns:
%   data (cell): a shape-(M, 1) cell array of data points generated from
%     the supplied histograms.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  counts (:,:) {mustBeNumeric,mustBeNonnegative}
  bins (1,:) {mustBeNumeric}
end

data = {};
for iHisto = 1:size(counts,1)
  histoData = [];
  for iBin = 1:length(counts)
    binData = repmat(bins(iBin), 1, counts(iHisto,iBin));
    histoData = [histoData, binData]; %#ok<*AGROW>
  end
  data{iHisto} = histoData;
end