function outputPath = writeUnitsTableH5(nwbPath, unitsData)
% writeUnitsTableH5(nwbPath, unitsData)
%
% Creates /units inside an existing NWB file, in place, via low-level HDF5
% - no file copy, no MatNWB, no nwbRead/nwbExport. This is cellexplorer2nwb.m's
% write step, factored out and ported off MatNWB's general-purpose export
% path (see cellexplorer2nwb.m for the data-preparation half that builds
% unitsData).
%
% The on-disk layout - the typed group plus its colnames/description/
% namespace/neurodata_type/object_id attributes, an id ElementIdentifiers
% dataset, one VectorData dataset per scalar/string column, the
% spike_times/spike_times_index ragged-array pair (a VectorIndex whose
% 'target' attribute is an HDF5 object reference to the spike_times
% dataset), the electrode_group column (one HDF5 object reference per unit
% into /general/extracellular_ephys/<group_name>), and the 2-D
% waveform_mean/waveform_sd matrices - was verified byte-for-byte against a
% real MatNWB-written /units group (sub-MSEL02793 ..._ieeg.nwb, 17 units)
% before this was written, including the object-reference machinery (empty
% test file: two dummy groups referenced from a 3-element array dataset and
% from a scalar attribute, round-tripped through H5R.create/H5R.get_name)
% and the 2-D dataspace dimension-order convention for H5D.write (a MATLAB
% [num_samples, num_units] array requires the low-level H5S dims to be
% given as [num_units, num_samples] - i.e. fliplr(size(array)) - with the
% array passed to H5D.write UNTRANSPOSED; confirmed empirically, not by
% inference from documentation).
%
% Writing the waveform matrices via a real 2-D dataspace (rather than
% relying on MatNWB's generic VectorData export path, which flattens a
% singleton dimension to 1-D) means the single-unit case never needs the
% post-hoc patch cellexplorer2nwb.m's local_force_units_dim used to apply.
%
% If /units already exists it is deleted and rewritten (idempotent re-run,
% matching cellexplorer2nwb.m's existing behaviour).
%
% Args:
%   nwbPath (char, required, positional): path to the existing NWB file
%     (edited in place).
%   unitsData (struct, required, positional): a shape-(1,1) scalar struct,
%     one row per unit throughout (N = number of units):
%     clusterId          - Nx1 double, unique cluster id
%     type               - Nx1 cellstr, 'unit' or 'mua'
%     peakChannelIndex   - Nx1 double, 0-based row index into the
%                          electrodes table (NWB convention). Files injected
%                          before 2026-09-29 stored it 1-based.
%     peakChannelId      - Nx1 cellstr, e.g. 'LTAMicro_0003'
%     x, y, z            - Nx1 double, unit coordinates
%     area               - Nx1 cellstr, e.g. 'Left_Hippocampus'
%     leadId             - Nx1 cellstr, e.g. 'LTA'
%     electrodeGroupName - Nx1 cellstr, the electrode group each unit
%                          belongs to (must match a
%                          /general/extracellular_ephys/<name> group
%                          already present in the file)
%     spikeTimes         - 1xN cell array; spikeTimes{i} is a column vector
%                          of unit i's spike times (seconds, on the NWB
%                          timeline)
%     waveformMean       - [nSamplesMean x N] double
%     waveformSd         - [nSamplesSd x N] double
%
% Returns:
%   outputPath - nwbPath (char), for call-site convenience.
%
% See also: cellexplorer2nwb, resolveTimestampsAtIndicesH5,
%   annotatedSz2nwb/write_time_intervals_h5.m (the same low-level strategy,
%   applied to /intervals tables)
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  nwbPath (1,:) {mustBeA(nwbPath,'char'), mustBeVector, mustBeFile}
  unitsData (1,1) struct
end

outputPath = nwbPath;

req = {'clusterId','type','peakChannelIndex','peakChannelId','x','y','z', ...
    'area','leadId','electrodeGroupName','spikeTimes','waveformMean','waveformSd'};
for i = 1:numel(req)
  if ~isfield(unitsData, req{i})
    error('writeUnitsTableH5:missingField', 'unitsData is missing field ''%s''.', req{i});
  end
end

n = numel(unitsData.clusterId);
colnames = {'cluster_id','type','peak_channel_index','peak_channel_id','x','y','z', ...
    'area','lead_id','electrode_group','spike_times','waveform_mean','waveform_sd'};

% ── Open file read-write ────────────────────────────────────────────────
fileId = H5F.open(nwbPath, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
cleanupFile = onCleanup(@() H5F.close(fileId));

% ── Replace any existing /units (idempotent re-run) ────────────────────
if H5L.exists(fileId, '/units', 'H5P_DEFAULT')
  H5L.delete(fileId, '/units', 'H5P_DEFAULT');
end

grpId = H5G.create(fileId, '/units', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
cleanupGrp = onCleanup(@() H5G.close(grpId));

local_write_str_attr(grpId, 'namespace', 'core');
local_write_str_attr(grpId, 'neurodata_type', 'Units');
local_write_str_attr(grpId, 'description', 'Units table');
local_write_str_attr(grpId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_array_attr(grpId, 'colnames', colnames);

% ── id (ElementIdentifiers): 0..N-1, no description attribute ─────────
idId = local_create_int64_dataset(grpId, 'id', int64(0:n-1)');
local_write_str_attr(idId, 'namespace', 'hdmf-common');
local_write_str_attr(idId, 'neurodata_type', 'ElementIdentifiers');
local_write_str_attr(idId, 'object_id', char(matlab.lang.internal.uuid()));
H5D.close(idId);

% ── Scalar/string VectorData columns ────────────────────────────────────
local_write_vector_data(grpId, 'cluster_id', double(unitsData.clusterId(:)), ...
    'Unique cluster id');
local_write_vector_data(grpId, 'type', unitsData.type(:), ...
    'Cluster type: unit vs mua');
local_write_vector_data(grpId, 'peak_channel_index', double(unitsData.peakChannelIndex(:)), ...
    'Peak channel row index in the electrode table');
local_write_vector_data(grpId, 'peak_channel_id', unitsData.peakChannelId(:), ...
    'Unique ID of the channel with the largest cluster waveform amplitude');
local_write_vector_data(grpId, 'x', double(unitsData.x(:)), 'x coordinate');
local_write_vector_data(grpId, 'y', double(unitsData.y(:)), 'y coordinate');
local_write_vector_data(grpId, 'z', double(unitsData.z(:)), 'z coordinate');
local_write_vector_data(grpId, 'area', unitsData.area(:), ...
    'Brain area where the unit is located.');
local_write_vector_data(grpId, 'lead_id', unitsData.leadId(:), ...
    'Lead id where the unit is located');

% ── spike_times (VectorData) + spike_times_index (VectorIndex) ─────────
spikeTimesFlat = cat(1, unitsData.spikeTimes{:});
if isempty(spikeTimesFlat), spikeTimesFlat = zeros(0,1); end
stId = local_create_double_dataset(grpId, 'spike_times', spikeTimesFlat);
local_write_str_attr(stId, 'namespace', 'hdmf-common');
local_write_str_attr(stId, 'neurodata_type', 'VectorData');
local_write_str_attr(stId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(stId, 'description', 'Single unit activity');
H5D.close(stId);

counts = cellfun(@numel, unitsData.spikeTimes(:));
cumulativeIndex = uint64(cumsum(counts));
idxId = local_create_uint64_dataset(grpId, 'spike_times_index', cumulativeIndex);
local_write_str_attr(idxId, 'namespace', 'hdmf-common');
local_write_str_attr(idxId, 'neurodata_type', 'VectorIndex');
local_write_str_attr(idxId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(idxId, 'description', 'indexes data');
local_write_ref_attr(idxId, 'target', fileId, '/units/spike_times');
H5D.close(idxId);

% ── electrode_group: one HDF5 object reference per unit ─────────────────
egId = local_create_ref_array_dataset(grpId, fileId, 'electrode_group', ...
    unitsData.electrodeGroupName(:));
local_write_str_attr(egId, 'namespace', 'hdmf-common');
local_write_str_attr(egId, 'neurodata_type', 'VectorData');
local_write_str_attr(egId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(egId, 'description', 'Recording channel groups');
H5D.close(egId);

% ── waveform_mean / waveform_sd: real 2-D datasets from the start ──────
wmId = local_create_matrix_dataset(grpId, 'waveform_mean', double(unitsData.waveformMean));
local_write_str_attr(wmId, 'namespace', 'hdmf-common');
local_write_str_attr(wmId, 'neurodata_type', 'VectorData');
local_write_str_attr(wmId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(wmId, 'description', ['Mean waveforms on the probe channel with the ' ...
    'largest waveform amplitude. The order that waveforms are stored match the order of ' ...
    'units in the unit table.']);
local_write_str_attr(wmId, 'unit', 'volts');
H5D.close(wmId);

wsId = local_create_matrix_dataset(grpId, 'waveform_sd', double(unitsData.waveformSd));
local_write_str_attr(wsId, 'namespace', 'hdmf-common');
local_write_str_attr(wsId, 'neurodata_type', 'VectorData');
local_write_str_attr(wsId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(wsId, 'description', 'Standard deviation of waveforms.');
local_write_str_attr(wsId, 'unit', 'volts');
H5D.close(wsId);

clear cleanupGrp cleanupFile;  % H5G.close then H5F.close

fprintf('  Wrote /units (%d units)\n', n);


%% ── Local helpers ───────────────────────────────────────────────────────

function local_write_vector_data(grpId, name, values, description)
if iscell(values)
  values = cellfun(@(x) char(string(x)), values, 'UniformOutput', false);
  dsId = local_create_vlen_string_dataset(grpId, name, values);
else
  dsId = local_create_double_dataset(grpId, name, values);
end
local_write_str_attr(dsId, 'namespace', 'hdmf-common');
local_write_str_attr(dsId, 'neurodata_type', 'VectorData');
local_write_str_attr(dsId, 'object_id', char(matlab.lang.internal.uuid()));
local_write_str_attr(dsId, 'description', description);
H5D.close(dsId);


function dsetId = local_create_double_dataset(locId, name, values)
n = numel(values);
typeId = H5T.copy('H5T_NATIVE_DOUBLE');
spaceId = H5S.create_simple(1, n, n);
dsetId = H5D.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, typeId, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', double(values(:)));
H5S.close(spaceId); H5T.close(typeId);


function dsetId = local_create_int64_dataset(locId, name, values)
n = numel(values);
typeId = H5T.copy('H5T_NATIVE_INT64');
spaceId = H5S.create_simple(1, n, n);
dsetId = H5D.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, typeId, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int64(values(:)));
H5S.close(spaceId); H5T.close(typeId);


function dsetId = local_create_uint64_dataset(locId, name, values)
% MatNWB writes VectorIndex data (e.g. spike_times_index) as UNSIGNED
% 64-bit - confirmed against a real MatNWB-written file (sub-MSEL02799
% ..._ieeg.nwb, /units/spike_times_index: H5T_INTEGER, size 8, sign
% unsigned). A signed int64 of the same width fails pynwb-validate
% ("incorrect type - expected 'uint8' [i.e. unsigned], got 'int64'").
n = numel(values);
typeId = H5T.copy('H5T_NATIVE_UINT64');
spaceId = H5S.create_simple(1, n, n);
dsetId = H5D.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, typeId, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', uint64(values(:)));
H5S.close(spaceId); H5T.close(typeId);


function dsetId = local_create_vlen_string_dataset(locId, name, values)
n = numel(values);
typeId = H5T.copy('H5T_C_S1');
H5T.set_size(typeId, 'H5T_VARIABLE');
H5T.set_cset(typeId, H5ML.get_constant_value('H5T_CSET_UTF8'));
spaceId = H5S.create_simple(1, n, n);
dsetId = H5D.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, typeId, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', values(:));
H5S.close(spaceId); H5T.close(typeId);


function dsetId = local_create_matrix_dataset(locId, name, data)
% data is a MATLAB [nRows x nCols] array (e.g. [num_samples, num_units]).
% Empirically-verified convention: declare the low-level H5S dataspace
% dims as fliplr(size(data)) and pass data to H5D.write UNTRANSPOSED - see
% this file's top-level comment for how that was confirmed.
sz = size(data);
if numel(sz) < 2, sz = [sz 1]; end
dims = fliplr(sz);
typeId = H5T.copy('H5T_NATIVE_DOUBLE');
spaceId = H5S.create_simple(2, dims, dims);
dsetId = H5D.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, typeId, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
H5S.close(spaceId); H5T.close(typeId);


function dsetId = local_create_ref_array_dataset(grpId, fileId, name, groupNames)
% One HDF5 object reference per row, into /general/extracellular_ephys/<name>.
n = numel(groupNames);
refType = H5T.copy('H5T_STD_REF_OBJ');
refSize = H5T.get_size(refType);
refs = zeros(refSize, n, 'uint8');
for i = 1:n
  targetPath = ['/general/extracellular_ephys/' groupNames{i}];
  if ~H5L.exists(fileId, targetPath, 'H5P_DEFAULT')
    error('writeUnitsTableH5:missingElectrodeGroup', ...
        'Electrode group ''%s'' referenced by unit %d does not exist in the file (%s).', ...
        groupNames{i}, i, targetPath);
  end
  refs(:,i) = H5R.create(fileId, targetPath, 'H5R_OBJECT', -1);
end
spaceId = H5S.create_simple(1, n, n);
dsetId = H5D.create(grpId, name, refType, spaceId, 'H5P_DEFAULT');
H5D.write(dsetId, refType, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', refs);
H5S.close(spaceId); H5T.close(refType);


function local_write_ref_attr(objId, name, fileId, targetPath)
% Scalar HDF5 object-reference attribute.
if ~H5L.exists(fileId, targetPath, 'H5P_DEFAULT')
  error('writeUnitsTableH5:missingRefTarget', ...
      'Reference target ''%s'' does not exist in the file.', targetPath);
end
refType = H5T.copy('H5T_STD_REF_OBJ');
spaceId = H5S.create('H5S_SCALAR');
ref = H5R.create(fileId, targetPath, 'H5R_OBJECT', -1);
attrId = H5A.create(objId, name, refType, spaceId, 'H5P_DEFAULT');
H5A.write(attrId, refType, ref);
H5A.close(attrId); H5S.close(spaceId); H5T.close(refType);


function local_write_str_attr(locId, name, value)
% Scalar variable-length UTF-8 string attribute (matches MatNWB's own
% datatype/cset convention - see annotatedSz2nwb/write_time_intervals_h5.m).
typeId = H5T.copy('H5T_C_S1');
H5T.set_size(typeId, 'H5T_VARIABLE');
H5T.set_cset(typeId, H5ML.get_constant_value('H5T_CSET_UTF8'));
spaceId = H5S.create('H5S_SCALAR');
attrId = H5A.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5A.write(attrId, typeId, {char(value)});
H5A.close(attrId); H5S.close(spaceId); H5T.close(typeId);


function local_write_str_array_attr(locId, name, values)
% 1-D array of variable-length UTF-8 strings (used for 'colnames').
n = numel(values);
typeId = H5T.copy('H5T_C_S1');
H5T.set_size(typeId, 'H5T_VARIABLE');
H5T.set_cset(typeId, H5ML.get_constant_value('H5T_CSET_UTF8'));
spaceId = H5S.create_simple(1, n, n);
attrId = H5A.create(locId, name, typeId, spaceId, 'H5P_DEFAULT');
H5A.write(attrId, typeId, cellstr(values(:)));
H5A.close(attrId); H5S.close(spaceId); H5T.close(typeId);
