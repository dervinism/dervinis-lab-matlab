function timestampsAtIndex = resolveTimestampsAtIndices(timeseriesData, indices)
% resolveTimestampsAtIndices(timeseriesData, indices)
%
% Function resolves the true NWB timestamp for each of INDICES (1-based
% sample indices into a TimeSeries' full data array). Kilosort spike
% sorting assumes a fixed sampling rate that might not always agree with
% the timestamps actually recorded in the NWB file, so this looks up the
% real per-sample timestamp when the TimeSeries has an explicit
% 'timestamps' array, or falls back to starting_time + rate arithmetic
% when it doesn't.
%
% This is used by both cellexplorer2nwb.m (to place spike times on the
% true NWB timeline) and annotatedSz2nwb's load_zeroed_periods.m (to
% apply the exact same correction to zeroed/artifact period boundaries).
% Kept as one shared implementation so the two pipelines cannot
% independently drift out of sync with each other.
%
% Args:
%   timeseriesData (required, positional): a matnwb TimeSeries (or
%     subclass, e.g. ElectricalSeries) object, as returned by
%     nwb.acquisition.get(...).
%   indices (numeric, required, positional): a shape-(n, 1) numeric
%     column vector of 1-based sample indices into timeseriesData's full
%     data array.
%
% Returns:
%   timestampsAtIndex (numeric): a shape-(n, 1) numeric column vector of
%     the corrected NWB timestamps (seconds) for each index.
%
% Dependencies:
%   matnwb (https://neurodatawithoutborders.github.io/matnwb/)
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  timeseriesData
  indices (:,1) {mustBeNumeric}
end

timestampsAtIndex = nan(size(indices));
if isempty(indices)
  return
end

if ~isempty(timeseriesData.timestamps)
  timestampsAtIndex = local_load_timestamps_at_indices(timeseriesData.timestamps, indices);
else
  t0 = timeseriesData.starting_time;        % seconds, relative to 12AM of the session_start_time day
  rate = timeseriesData.starting_time_rate; % Hz (samples per second)
  timestampsAtIndex = t0 + (double(indices) - 1)/rate;
end


%% Helper functions
function timestampsAtIndex = local_load_timestamps_at_indices(timestampsDataStub, indices, options)
% LOCAL_LOAD_TIMESTAMPS_AT_INDICES  Gather only the per-sample timestamps at
%   INDICES (1-based, into the full session-length dataset) instead of
%   materializing the whole array. A session-length, per-sample timestamps
%   array (tens of kHz over hours) can be tens of GB - too large to even
%   preallocate in memory on long recordings, let alone read in one
%   h5read call. Spike sample indices are a tiny fraction of the full
%   extent, so this sweeps the dataset once in blocks and keeps only the
%   requested samples from each block.
%
%   Uses low-level HDF5 calls with a SINGLE persistent file/dataset handle
%   for the whole sweep, rather than DataStub/DataPipe's .load() (which
%   goes through h5read and reopens the file on every call). On a huge
%   file living on a slow externally-attached disk, hundreds of repeated
%   file-opens - not the actual data volume - is what made this
%   prohibitively slow; a persistent handle removes that overhead.

arguments
  timestampsDataStub
  indices (:,1) {mustBeNumeric}
  options.chunkSize (1,1) {isnumeric} = 2e8 % ~1.6GB per chunk of doubles
end

timestampsAtIndex = nan(size(indices));
if isempty(indices)
  return
end

[filename, datasetPath, n] = local_resolve_h5_dataset(timestampsDataStub);

fileId = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fileCleanup = onCleanup(@() H5F.close(fileId)); %#ok<NASGU>
datasetId = H5D.open(fileId, datasetPath);
datasetCleanup = onCleanup(@() H5D.close(datasetId)); %#ok<NASGU>
fileSpaceId = H5D.get_space(datasetId);
spaceCleanup = onCleanup(@() H5S.close(fileSpaceId)); %#ok<NASGU>

nChunks = ceil(n / options.chunkSize);
ticStart = tic;
for iChunk = 1:nChunks
  s = (iChunk-1)*options.chunkSize + 1;
  c = min(options.chunkSize, n - s + 1);
  inChunk = indices >= s & indices < s + c;
  if any(inChunk)
    H5S.select_hyperslab(fileSpaceId, 'H5S_SELECT_SET', s-1, [], c, []);
    memSpaceId = H5S.create_simple(1, c, c);
    chunkData = H5D.read(datasetId, 'H5ML_DEFAULT', memSpaceId, fileSpaceId, 'H5P_DEFAULT');
    H5S.close(memSpaceId);
    timestampsAtIndex(inChunk) = chunkData(indices(inChunk) - s + 1);
  end
  fprintf('local_load_timestamps_at_indices: chunk %d/%d (%.1f%%), elapsed %.0f s\n', ...
    iChunk, nChunks, 100*iChunk/nChunks, toc(ticStart));
end


function [filename, datasetPath, n] = local_resolve_h5_dataset(dataStubOrPipe)
% LOCAL_RESOLVE_H5_DATASET  Extract the on-disk filename, absolute in-file
%   HDF5 dataset path, and element count (along dimension 1) from either a
%   types.untyped.DataStub or a file-bound types.untyped.DataPipe, so a
%   low-level H5F/H5D handle can be opened directly against the dataset.
if isa(dataStubOrPipe, 'types.untyped.DataStub')
  filename = dataStubOrPipe.filename;
  datasetPath = dataStubOrPipe.path;
  n = dataStubOrPipe.dims(1);
elseif isa(dataStubOrPipe, 'types.untyped.DataPipe')
  internal = dataStubOrPipe.internal;
  assert(isa(internal, 'types.untyped.datapipe.BoundPipe'), ...
    'NWB:ResolveTimestampsAtIndices:UnboundDataPipe', ...
    ['Expected a file-bound DataPipe when reading timestamps from an ' ...
    'existing NWB file, but found an unbound (in-memory only) pipe.']);
  filename = internal.filename;
  datasetPath = internal.path;
  n = internal.dims(1);
else
  error('NWB:ResolveTimestampsAtIndices:UnsupportedTimestampsType', ...
    'Unsupported timestamps object type: %s', class(dataStubOrPipe));
end
