function [binCounts, binLocs, totalCounts, significantFractions, ...
  phaseMeans, phaseSDs] = phaseHistrogram(phase, options)
% [binCounts, binLocs, totalCounts, significantFractions, phaseMeans, ...
%   phaseSDs] = phaseHistrogram(phase, <options>)
%
% Function produces a phase histogram (bin counts of phase values, not an
% actual graph).
%
% Args:
%   phase (numeric, required, positional): a shape-(M, N) numeric array of
%     phase values in radians. The histogram bins values column-wise.
%   centre (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing the centre of the new phase range for the phase
%     histogram (default = 0).
%   nBins (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     indicating the numer of histogram bins (default = 10).
%
% Returns:
%   binCounts (numeric): a shape-(nBins, N) numeric array with phase value
%     counts in each bin.
%   binLocs (numeric): a shape-(1, nBins) numeric array with phase
%     histogram bin locations in radians.
%   totalCounts (numeric): a shape-(1, N) numeric array of non-NaN counts
%     (over phase matrix columns).
%   significantFractions (numeric): a shape-(1, N) numeric array of non-NaN
%     count fractions (over phase matrix columns): A fraction of
%     non-NaN (significant) phase valuesout of all values supplied.
%   phaseMeans (numeric): a shape-(1, N) numeric array of phase means (over
%     phase matrix columns).
%   phaseSDs (numeric): a shape-(1, N) numeric array of phase standard
%     deviations (over phase matrix columns).
%
% Dependencies
%   
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%   dervinism/circStatNP (https://github.com/dervinism/circStatNP).
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  phase (:,:) {mustBeNumeric,mustBeNonempty}
  options.centre (1,1) {mustBeNumeric,mustBeNonNan,mustBeReal} = 0
  options.nBins (1,1) {mustBeNumeric,mustBePositive} = 10
end

% Recentre phase values
[phase, phaseRange] = recentrePhase(phase, options.centre);

% Work out bin locations
binSize = (phaseRange(2) - phaseRange(1))/options.nBins;
binLocs = phaseRange(1)+binSize/2:binSize:phaseRange(2);

% Bin phase values
if size(phase,1) == 1
  binCounts = hist([phase; phase], binLocs)./2; % deal with the singular case
else
  binCounts = hist(phase, binLocs); %#ok<*HIST>
end

% Calculate phase means and all values per column (excludes NaNs)
[phaseMeans, ~, totalCounts, phaseSDs] = datamean(phase, 'circularNP');
significantFractions = totalCounts./size(phase,1);