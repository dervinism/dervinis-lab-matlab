function [waveformMean, waveformSd, nSpikesUsed] = computeUnitWaveformsH5(nwbPath, acqKey, spikeIndices, acqColumns, options)
% [waveformMean, waveformSd, nSpikesUsed] = computeUnitWaveformsH5(nwbPath, acqKey, spikeIndices, acqColumns, <options>)
%
% Computes per-unit mean and standard deviation spike waveforms directly
% from an NWB file's raw acquisition data, via low-level HDF5 (no MatNWB).
% Used by cellexplorer2nwb.m to fill /units waveform_mean and waveform_sd
% in physical units (volts). These replace the CellExplorer waveforms,
% which were extracted from Kilosort's whitened binary (temp_wh.dat) and
% are therefore not in physical units, and whose SD window (52 samples)
% did not match the mean's (78 samples).
%
% For each unit, up to nSpikes spikes are drawn at random (fixed seed, so
% the result is reproducible) from its full spike train. A padded stretch
% of the unit's peak channel is read around each spike (spikes close
% together share one read), band-pass filtered (3rd-order Butterworth,
% zero-phase - the same filter nwb2binary.m applies before sorting), and a
% window of [samplesBefore samplesAfter] samples around the spike sample
% is kept. The mean and SD are taken across spikes, sample by sample.
% Spikes are not re-aligned: the stored waveform is what the stored spike
% times point to.
%
% Args:
%   nwbPath (char, required, positional): path to the NWB file.
%   acqKey (char, required, positional): acquisition group name, e.g.
%     'TimeSeries_32000_Hz'.
%   spikeIndices (cell, required, positional): a shape-(1, N) cell array
%     of 1-based sample indices into the acquisition data, one per unit.
%   acqColumns (numeric, required, positional): a shape-(1, N) vector of
%     1-based channel (column) positions within the acquisition data, one
%     per unit (the unit's peak channel).
%   nSpikes (numeric, optional, keyword): maximum number of spikes per unit
%     (default 2000).
%   samplesBefore (numeric, optional, keyword): samples before the spike
%     sample (default 40, so the spike sits at 0-based index 40).
%   samplesAfter (numeric, optional, keyword): samples after the spike
%     sample (default 37; 40 + 1 + 37 = 78 samples in total).
%   padSamples (numeric, optional, keyword): extra samples read on either
%     side of each stretch so filter edge effects never reach a waveform
%     (default 1600, i.e. 50 ms at 32 kHz).
%   frequencyRange (numeric, optional, keyword): band-pass range in Hz
%     (default [300 6000]).
%   sampleRate (numeric, optional, keyword): sampling rate in Hz used for
%     filtering (default 32000).
%   seed (numeric, optional, keyword): random seed (default 0).
%
% Returns:
%   waveformMean (numeric): a shape-(N, nWindow) array of mean waveforms
%     in volts.
%   waveformSd (numeric): a shape-(N, nWindow) array of waveform standard
%     deviations in volts.
%   nSpikesUsed (numeric): a shape-(1, N) vector of spikes averaged per unit.
%
% Dependencies:
%   bandpassFilterTimeSeries.m (this repo).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  nwbPath (1,:) {mustBeA(nwbPath,'char'), mustBeVector, mustBeFile}
  acqKey (1,:) {mustBeA(acqKey,'char'), mustBeVector}
  spikeIndices (1,:) cell
  acqColumns (1,:) {mustBeNumeric, mustBeInteger, mustBePositive}
  options.nSpikes (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 2000
  options.samplesBefore (1,1) {mustBeNumeric, mustBeNonnegative, mustBeInteger} = 40
  options.samplesAfter (1,1) {mustBeNumeric, mustBeNonnegative, mustBeInteger} = 37
  options.padSamples (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 1600
  options.frequencyRange (1,2) {mustBeNumeric, mustBePositive} = [300 6000]
  options.sampleRate (1,1) {mustBeNumeric, mustBePositive} = 32000
  options.seed (1,1) {mustBeNumeric} = 0
end

assert(numel(spikeIndices) == numel(acqColumns), ...
  'computeUnitWaveformsH5:sizeMismatch', 'spikeIndices and acqColumns must have one entry per unit.');

dataPath = ['/acquisition/' acqKey '/data'];
dataSize = local_h5_dataset_size(nwbPath, dataPath);
[timeAxis, nSamplesTotal] = local_time_axis(dataSize, max(acqColumns));
try
  conversion = double(h5readatt(nwbPath, dataPath, 'conversion'));
catch
  conversion = 1;
end
try
  offsetValue = double(h5readatt(nwbPath, dataPath, 'offset'));
catch
  offsetValue = 0;
end

nUnits = numel(spikeIndices);
nWindow = options.samplesBefore + 1 + options.samplesAfter;
waveformMean = nan(nUnits, nWindow);
waveformSd = nan(nUnits, nWindow);
nSpikesUsed = zeros(1, nUnits);
maxGap = 2*options.padSamples; % spikes closer than this share one read
lowest = options.samplesBefore + options.padSamples + 1;
highest = nSamplesTotal - options.samplesAfter - options.padSamples;

for iUnit = 1:nUnits
  inds = sort(double(spikeIndices{iUnit}(:))); % reads below assume ascending order
  inds = inds(inds >= lowest & inds <= highest);
  if isempty(inds)
    warning('computeUnitWaveformsH5:noSpikes', 'Unit %d has no usable spikes; waveform left as NaN.', iUnit);
    continue
  end
  rng(options.seed + iUnit);
  if numel(inds) > options.nSpikes
    inds = inds(sort(randperm(numel(inds), options.nSpikes)));
  else
    inds = sort(inds);
  end

  snippets = zeros(numel(inds), nWindow);
  groupStarts = [1; find(diff(inds) > maxGap) + 1];
  groupEnds = [groupStarts(2:end) - 1; numel(inds)];
  for iGroup = 1:numel(groupStarts)
    groupInds = inds(groupStarts(iGroup):groupEnds(iGroup));
    readStart = groupInds(1) - options.samplesBefore - options.padSamples;
    readCount = groupInds(end) + options.samplesAfter + options.padSamples - readStart + 1;
    if timeAxis == 2
      raw = h5read(nwbPath, dataPath, [acqColumns(iUnit), readStart], [1, readCount]);
    else
      raw = h5read(nwbPath, dataPath, [readStart, acqColumns(iUnit)], [readCount, 1]);
    end
    stretch = double(raw(:)')*conversion + offsetValue;
    stretch = bandpassFilterTimeSeries(stretch, sampleRate=options.sampleRate, ...
      frequencyRange=options.frequencyRange);
    for iSpike = 1:numel(groupInds)
      centre = groupInds(iSpike) - readStart + 1;
      snippets(groupStarts(iGroup)+iSpike-1, :) = ...
        stretch(centre-options.samplesBefore:centre+options.samplesAfter);
    end
  end

  waveformMean(iUnit,:) = mean(snippets, 1);
  waveformSd(iUnit,:) = std(snippets, 0, 1);
  nSpikesUsed(iUnit) = numel(inds);
  fprintf('computeUnitWaveformsH5: unit %d/%d - %d spikes, trough %.1f uV\n', ...
    iUnit, nUnits, nSpikesUsed(iUnit), 1e6*min(waveformMean(iUnit,:)));
end
end


%% Helper functions
function [timeAxis, nSamplesTotal] = local_time_axis(dataSize, maxColumn)
% LOCAL_TIME_AXIS  Which MATLAB-view dimension of the acquisition data is
%   time. The channel dimension is the smaller one and must hold every
%   requested column.
if dataSize(1) < dataSize(2)
  timeAxis = 2;
else
  timeAxis = 1;
end
nSamplesTotal = dataSize(timeAxis);
assert(maxColumn <= dataSize(3-timeAxis), 'computeUnitWaveformsH5:badColumn', ...
  'Requested channel column %d exceeds the %d channels in the data.', maxColumn, dataSize(3-timeAxis));
end


function sz = local_h5_dataset_size(h5Path, datasetPath)
% LOCAL_H5_DATASET_SIZE  Dataset dimensions via low-level HDF5, in h5read's
%   (MATLAB) order. Avoids h5info, which inspects the full chunk layout and
%   can take many minutes on these finely chunked acquisition datasets.
fileId = H5F.open(h5Path, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fileCleanup = onCleanup(@() H5F.close(fileId)); %#ok<NASGU>
dsetId = H5D.open(fileId, datasetPath);
dsetCleanup = onCleanup(@() H5D.close(dsetId)); %#ok<NASGU>
spaceId = H5D.get_space(dsetId);
spaceCleanup = onCleanup(@() H5S.close(spaceId)); %#ok<NASGU>
[~, dims, ~] = H5S.get_simple_extent_dims(spaceId);
sz = fliplr(dims);
end
