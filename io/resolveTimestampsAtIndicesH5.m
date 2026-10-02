function timestampsAtIndex = resolveTimestampsAtIndicesH5(nwbPath, acqKey, indices)
% resolveTimestampsAtIndicesH5(nwbPath, acqKey, indices)
%
% Low-level-HDF5 equivalent of resolveTimestampsAtIndices.m: resolves the
% true NWB timestamp for each of INDICES (1-based sample indices into an
% acquisition TimeSeries' full data array), without requiring a matnwb
% TimeSeries object (or MatNWB on the path at all). Reads directly from
% nwbPath via h5info/h5readatt/H5D, matching the exact algorithm of
% resolveTimestampsAtIndices.m's local_load_timestamps_at_indices (a single
% persistent file/dataset handle, swept once in chunks, so a session-length
% per-sample timestamps array is never materialized in full) - the only
% difference is how the timestamps dataset (or starting_time/rate fallback)
% is located.
%
% Args:
%   nwbPath (char, required, positional): path to the NWB file.
%   acqKey (char, required, positional): name of the group under
%     /acquisition, e.g. 'TimeSeries_32000_Hz'.
%   indices (numeric, required, positional): a shape-(n, 1) numeric column
%     vector of 1-based sample indices into the TimeSeries' full data array.
%
% Returns:
%   timestampsAtIndex (numeric): a shape-(n, 1) numeric column vector of
%     the corrected NWB timestamps (seconds) for each index.
%
% See also: resolveTimestampsAtIndices (the matnwb-object-based original;
%   this file mirrors its algorithm and both must be kept in sync)
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  nwbPath (1,:) {mustBeA(nwbPath,'char'), mustBeVector, mustBeFile}
  acqKey (1,:) {mustBeA(acqKey,'char'), mustBeVector}
  indices (:,1) {mustBeNumeric}
end

timestampsAtIndex = nan(size(indices));
if isempty(indices)
  return
end

base = ['/acquisition/' acqKey];
timestampsPath = [base '/timestamps'];

if local_h5_exists(nwbPath, timestampsPath)
  n = local_dataset_size1(nwbPath, timestampsPath);
  timestampsAtIndex = local_load_timestamps_at_indices(nwbPath, timestampsPath, n, indices);
else
  t0   = double(h5read(nwbPath, [base '/starting_time']));
  rate = double(h5readatt(nwbPath, [base '/starting_time'], 'rate'));
  timestampsAtIndex = t0 + (double(indices) - 1)/rate;
end


%% Helper functions
function timestampsAtIndex = local_load_timestamps_at_indices(nwbPath, datasetPath, n, indices, options)
% Sweeps the dataset once in blocks and keeps only the requested samples
% from each block, via a single persistent H5F/H5D handle for the whole
% sweep - same rationale as resolveTimestampsAtIndices.m's helper of the
% same purpose (repeated file-opens, not data volume, dominate on a slow
% network share).

arguments
  nwbPath (1,:) char
  datasetPath (1,:) char
  n (1,1) {isnumeric}
  indices (:,1) {mustBeNumeric}
  options.chunkSize (1,1) {isnumeric} = 2e8 % ~1.6GB per chunk of doubles
end

timestampsAtIndex = nan(size(indices));

fileId = H5F.open(nwbPath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
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
  fprintf('resolveTimestampsAtIndicesH5: chunk %d/%d (%.1f%%), elapsed %.0f s\n', ...
    iChunk, nChunks, 100*iChunk/nChunks, toc(ticStart));
end


function n = local_dataset_size1(nwbPath, datasetPath)
% Dimension-1 size via H5D.get_space + H5S.get_simple_extent_dims, NOT
% h5info - h5info eagerly inspects a dataset's full storage/chunk layout
% and can take minutes on a huge, finely-chunked dataset like this one
% (see annotatedSz2nwb/h5_dataset_size.m, which documents the same issue;
% duplicated here rather than called, to keep this library free of a
% dependency on that downstream project).
fileId = H5F.open(nwbPath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fileCleanup = onCleanup(@() H5F.close(fileId)); %#ok<NASGU>
datasetId = H5D.open(fileId, datasetPath);
datasetCleanup = onCleanup(@() H5D.close(datasetId)); %#ok<NASGU>
spaceId = H5D.get_space(datasetId);
spaceCleanup = onCleanup(@() H5S.close(spaceId)); %#ok<NASGU>
[~, dims, ~] = H5S.get_simple_extent_dims(spaceId);
dims = fliplr(dims);
n = dims(1);


function tf = local_h5_exists(nwbPath, linkPath)
% Existence check that walks each path component (H5L.exists can error on a
% multi-level path whose parent is missing) - mirrors annotatedSz2nwb/
% check_nwb.m's local_h5_exists.
fid = H5F.open(nwbPath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
c = onCleanup(@() H5F.close(fid)); %#ok<NASGU>
parts = strsplit(linkPath, '/');
cur = ''; tf = true;
for i = 2:numel(parts)
  cur = [cur '/' parts{i}]; %#ok<AGROW>
  if ~H5L.exists(fid, cur, 'H5P_DEFAULT'), tf = false; return; end
end
