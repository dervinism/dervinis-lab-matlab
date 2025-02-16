function [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikes(spikeTimes, options)
% [resampledSpikes, spikeTimeBins, resampledSpikeTimes] = resampleSpikes(spikeTimes, <options>)
%
% Function resamples spike times series and converts it into a spike count
% series.
%
% Args:
%   spikeTimes (numeric, required, positional): a shape-(1, N) numeric
%     array of spike times where N corresponds to spike times.
%   stepsize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the new sampling interval (default = 0.002).
%   startTime (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     corresponding to the start time bin (default = stepsize).
%
% Returns:
%   resampledSpikes (numeric): a shape-(1, M) numeric array of resampled
%     spike counts with the time bin of the same size as the input stepsize
%     variable.
%   spikeTimeBins (numeric): a shape-(1, M) numeric array with time bins
%     corresponding to resampledSpikes output variable.
%   resampledSpikeTimes (numeric): a shape-(1, L) numeric array of
%     resampled spike times.
%
% Dependencies:
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes (1,:) {mustBeNumeric,mustBeNonempty,mustBeNonnegative}
  options.stepsize (1,1) {mustBeNumeric,mustBePositive} = 0.002
  options.startTime (1,1) {mustBeNumeric} = 0
end


%% Generating continuous resampled spike container
spikeTimeBins = max([1 floor(options.startTime/options.stepsize)])*options.stepsize:options.stepsize:(getMaxSpikeTime(spikeTimes)+options.stepsize);
resampledSpikes = zeros(1,numel(spikeTimeBins));


%% Filling in the spike container with actual spikes
idx = round(spikeTimes/options.stepsize) - round(spikeTimeBins(1)/options.stepsize) + 1; % Assign bin indices to spikes
idx(idx == 0) = 1;
[spkCounts,idx2] = groupcounts(idx'); % Count spikes within bin indices
resampledSpikes(idx2') = spkCounts; % Mark spike presentations


%% Convert the continuous spike vector into spike times vector
resampledSpikeTimes = spikeTimeBins(idx);