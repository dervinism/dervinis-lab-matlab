function [resampledSpikeCounts, timeBins] = resampleSpikeCounts(spikeCounts, options)
% [resampledSpikeCounts, timeBins] = resampleSpikeCounts(spikeCounts, <options>)
%
% Function resamples spike count vectors in a row matrix. You can only
% downsample, not upsample.
%
% Args:
%   spikeCounts (numeric, required, positional): a shape-(M, N) numeric
%     array containing spike counts per sample point. Individual spike
%     count vectors correspond to individual rows in a matrix.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the sampling interval used in the spikeCounts vector
%     (default = 0.002).
%   newStepsize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the new sampling interval (default = 0.05).
%
% Returns:
%   resampledSpikeCounts (numeric) a shape-(M, L) numeric array containing
%     resampled spike count vectors (rows in a matrix).
%   timeBins (numeric) a shape-(1, L) numeric array time bins corresponding
%     to the columns of resampledSpikeCounts matrix. The time bins are
%     centred rather than being positioned at the bin edges.
%
% Dependencies:
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeCounts (:,:) {mustBeNumeric,mustBeNonempty,mustBeNonnegative}
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.newStepsize (1,1) {mustBeNumeric,mustBePositive} = 0.05
end

% parse input
if issparse(spikeCounts)
  spikeCounts = full(spikeCounts);
end
sr = 1/options.stepsize;
newsr = 1/options.newStepsize;
r = rem(sr, newsr);
if r
  warning(['Old sampling rate is not a multiplicative of the new sampling rate. ' ...
    'For loop will be used to downsample instead of vectorisation considerably slowing down calculations.']);
end

% Re-bin the spikes
if r
  oldTimeBins = (1:size(spikeCounts,2)).*options.stepsize;
  duration = oldTimeBins(end);
  newTimeBins = 0:options.newStepsize:duration;
  if newTimeBins(end) < duration
    newTimeBins = [newTimeBins newTimeBins(end)+options.newStepsize];
  end

  resampledSpikeCounts = zeros(size(spikeCounts,1),numel(newTimeBins)-1);
  for iBin = 1:numel(newTimeBins)-1
    %disp(['Progress: ' num2str((100*iBin)/numel(newTimeBins)) '%'])
    if mod(iBin, 2) == 0
      [~, binInds] = selectArrayValues(oldTimeBins, ...
        newTimeBins(iBin:iBin+1), cutoffType='exclusive');
    else
      [~, binInds] = selectArrayValues(oldTimeBins, ...
        newTimeBins(iBin:iBin+1), cutoffType='inclusive');
    end
    resampledSpikeCounts(:,iBin) = sum(spikeCounts(:,binInds),2);
  end

  timeBins = newTimeBins(1:end-1) + 0.5/newsr;

else
  durationNew = ceil(size(spikeCounts,2)/(sr/newsr));
  if size(spikeCounts,2) < durationNew*(sr/newsr)
    spikeCounts = [spikeCounts zeros(size(spikeCounts,1), durationNew*(sr/newsr) - size(spikeCounts,2))];
  end
  rasterHist = reshape(spikeCounts', [(sr/newsr) durationNew size(spikeCounts,1)]);

  resampledSpikeCounts = squeeze(sum(rasterHist,1))';
  timeBins = (1:durationNew).*(1/newsr) + 0.5/newsr;
end