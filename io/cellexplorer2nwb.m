function cellexplorer2nwb(spikeData, nwbFile, options)
% cellexplorer2nwb(spikeData, nwbFile, <options>)
%
% Function converts spike sorting output saved in the CellExplorer
% (https://cellexplorer.org/) format (spikes.cellinfo.mat) to NWB format by
% attaching spike times timeseries data to an exisiting raw NWB file.
%
% Args:
%   spikeData (struct, required, positional): a shape-(1,1) scalar
%     structure with the following fields:
%     files - a cell array of CellExplorer files (spikes.cellinfo.mat)
%       containing spike times of the indicated units in the field
%       existingUnitIDs.
%     startTimes - spike time displacements relative to the start of the
%       full recording file.
%     existingUnitIDs - unit IDs stored in CellExplorer files.
%     newGlobalUnitID - new unit IDs.
%     newGlobalUnitCh - new unit channels (channels can slightly vary
%       between chunk files).
%     chLabels - channel labels (e.g., LCmicro0001).
%     leadLabels - lead labels (e.g., LC).
%     areaLabels - area labels (e.g., Left_Hippocampus).
%   nwbFile (char, required, positional): a shape-(1, n) character array
%     containing the input NWB file name with a full path that needs
%     updating. It must end with '.nwb'.
%   outputFile (char, optional, keyword): a shape-(1, m) character array
%     containing the name of the updated NWB file in case it is different
%     to the original input file. It must end with '.nwb'.
%   loadCoordsFromUnitTable (logical, optional, keyword): a shape-(1, 1)
%     logical scalar indicating whether to load channel coordinates from a
%     unit table file instead of the NWB raw data file (default is false).
%
% Returns:
%   None.
%
% Dependencies:
%   matnwb (https://neurodatawithoutborders.github.io/matnwb/)
%   resolveTimestampsAtIndices.m (this repo's io folder)
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeData (1,1) {mustBeA(spikeData,'struct')}
  nwbFile (1,:) {mustBeA(nwbFile,'char'),mustBeVector,endsWith(nwbFile,'.nwb')}
  options.outputFile (1,:) {mustBeA(options.outputFile,'char'),mustBeVector,endsWith(options.outputFile,'.nwb')} = '';
  options.loadCoordsFromUnitTable (1,1) {islogical} = false
end

% Parse input
if isempty(options.outputFile)
  options.outputFile = nwbFile;
end

% Load spiking data and convert it
nUnits = numel(spikeData.existingUnitIDs);
if ~nUnits
  warning('Empty spikeData supplied. Exiting cellexplorer2nwb function.');
  return
end
loadedFile = '';
spikesConv = struct();
for iUnit = 1:nUnits
  fileToLoad = spikeData.files{iUnit};
  if ~strcmpi(loadedFile,fileToLoad)
    if exist('trackCheckedUnits','var')
      if numel(trackCheckedUnits) ~= ...
          numel(spikes.cluID(cellfun(@(x) strcmpi(x, 'unit'), spikes.labels))) && ...
          numel(trackCheckedUnits) ~= ...
          numel(spikes.cluID(cellfun(@(x) strcmpi(x, 'mua'), spikes.labels)))
        error(['Not all units from the previous CellExplorer file were checked: ' ...
          loadedFile]);
      end
    end
    load(fileToLoad); %#ok<*LOAD>
    loadedFile = fileToLoad;
    trackCheckedUnits = [];
  end
  unitIndExisting = find(spikeData.existingUnitIDs(iUnit) == spikes.cluID);
  trackCheckedUnits = [trackCheckedUnits spikes.cluID(unitIndExisting)]; %#ok<AGROW>
  if isempty(unitIndExisting)
    error(['Unit ' num2str(spikeData.existingUnitIDs(iUnit)) ...
      ' does not exist in the file ' fileToLoad]);
  end
  unitIndNew = spikeData.newGlobalUnitIDs(iUnit);
  if iUnit == 1 || unitIndNew > numel(spikesConv.cluID)
    spikesConv.times{unitIndNew} = spikes.times{unitIndExisting} + spikeData.startTimes(iUnit);
    spikesConv.cluID{unitIndNew} = unitIndNew;
    spikesConv.maxWaveformCh1{unitIndNew} = spikeData.newGlobalUnitCh(iUnit)+1;
    spikesConv.filtWaveform{unitIndNew} = ...
      spikes.rawWaveform_all{unitIndExisting}(spikesConv.maxWaveformCh1{unitIndNew},:);
    spikesConv.filtWaveformSD{unitIndNew} = spikes.rawWaveform_std{unitIndExisting};
    spikesConv.labels{unitIndNew} = spikes.labels{unitIndExisting};
    spikesConv.chLabels{unitIndNew} = spikeData.chLabels{iUnit};
    spikesConv.leadLabels{unitIndNew} = spikeData.leadLabels{iUnit};
    spikesConv.areaLabels{unitIndNew} = spikeData.areaLabels{iUnit};
    if options.loadCoordsFromUnitTable
      spikesConv.x{unitIndNew} = spikeData.x(iUnit);
      spikesConv.y{unitIndNew} = spikeData.y(iUnit);
      spikesConv.z{unitIndNew} = spikeData.z(iUnit);
    end
  elseif unitIndNew <= numel(spikesConv.cluID)
    spikesConv.times{unitIndNew} = [spikesConv.times{unitIndNew}; ...
      spikes.times{unitIndExisting} + spikeData.startTimes(iUnit)];
    spikesConv.cluID{unitIndNew} = unitIndNew;
    spikesConv.maxWaveformCh1{unitIndNew} = spikeData.newGlobalUnitCh(iUnit)+1;
    spikesConv.filtWaveform{unitIndNew} = ...
      spikes.rawWaveform_all{unitIndExisting}(spikesConv.maxWaveformCh1{unitIndNew},:);
    spikesConv.filtWaveformSD{unitIndNew} = spikes.rawWaveform_std{unitIndExisting};
    spikesConv.labels{unitIndNew} = spikes.labels{unitIndExisting};
    spikesConv.chLabels{unitIndNew} = spikeData.chLabels{iUnit};
    spikesConv.leadLabels{unitIndNew} = spikeData.leadLabels{iUnit};
    spikesConv.areaLabels{unitIndNew} = spikeData.areaLabels{iUnit};
    if options.loadCoordsFromUnitTable
      spikesConv.x{unitIndNew} = spikeData.x(iUnit);
      spikesConv.y{unitIndNew} = spikeData.y(iUnit);
      spikesConv.z{unitIndNew} = spikeData.z(iUnit);
    end
  else
    error('Unit IDs are not expected to decrease.');
  end
end

% Delete the existing units table
file_id = H5F.open(options.outputFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
if H5L.exists(file_id, '/units', 'H5P_DEFAULT')
    H5L.delete(file_id, '/units', 'H5P_DEFAULT');
end
H5F.close(file_id);

% Load NWB data
nwb = nwbRead(nwbFile);
electrodesTable = nwb.general_extracellular_ephys_electrodes.toTable();
timeseriesData = nwb.acquisition.get(['TimeSeries_' num2str(spikes.sr) '_Hz']);

% Remap spike times based on NWB timestamps because kilosort spike sorter
% uses fixed sampling rate which might not always agree with NWB timestamps.
% Only the specific per-sample timestamps a unit's spikes actually land on
% are fetched/computed - materializing the full per-sample timestamps
% vector for the whole session (previously done unconditionally here) ran
% out of memory on long, high sampling-rate recordings.
nUnits = numel(spikesConv.cluID);
spikesConv.timesNWB = cell(1,nUnits);
spikeIndsPerUnit = cell(1,nUnits);
for iUnit = 1:nUnits
  spikeIndsPerUnit{iUnit} = round(spikesConv.times{iUnit}./(1/spikes.sr));
end
allSpikeInds = cat(1, spikeIndsPerUnit{:});
if ~isempty(allSpikeInds)
  allTimestamps = resolveTimestampsAtIndices(timeseriesData, allSpikeInds);
  splitPoint = 0;
  for iUnit = 1:nUnits
    nSpikes = numel(spikeIndsPerUnit{iUnit});
    spikesConv.timesNWB{iUnit} = allTimestamps(splitPoint+1:splitPoint+nSpikes);
    splitPoint = splitPoint + nSpikes;
  end
end

% Obtain additional channel info
for iUnit = 1:nUnits
  channelInd = find(ismember(electrodesTable.ChName, spikesConv.chLabels{iUnit}));
  if isempty(channelInd)
    electrodesTable.ChName = cellfun(@(x) strrep(x, '_000', ''), electrodesTable.ChName, 'UniformOutput', false);
    electrodesTable.ChName = cellfun(@(x) strrep(x, '_00', ''), electrodesTable.ChName, 'UniformOutput', false);
    electrodesTable.ChName = cellfun(@(x) strrep(x, '_0', ''), electrodesTable.ChName, 'UniformOutput', false);
    electrodesTable.ChName = cellfun(@(x) strrep(x, '_', ''), electrodesTable.ChName, 'UniformOutput', false);
    electrodesTable.ChName = cellfun(@(x) strrep(x, 'Ch', ''), electrodesTable.ChName, 'UniformOutput', false);
    channelInd = find(ismember(electrodesTable.ChName, spikesConv.chLabels{iUnit}));
    assert(~isempty(channelInd));
  end
  spikesConv.channelInds{iUnit} = channelInd;
  if ~options.loadCoordsFromUnitTable
    spikesConv.x{iUnit} = electrodesTable.x(spikesConv.channelInds{iUnit});
    spikesConv.y{iUnit} = electrodesTable.y(spikesConv.channelInds{iUnit});
    spikesConv.z{iUnit} = electrodesTable.z(spikesConv.channelInds{iUnit});
  end
  spikesConv.group(iUnit) = electrodesTable.group(spikesConv.channelInds{iUnit});
end

% Create units table
% waveform_mean/waveform_sd below are transposed to [num_samples, num_units]:
% matnwb's DynamicTable row-height check for matrix-valued columns
% (types.util.dynamictable.internal.getColumnHeight) takes the LAST MATLAB
% dimension as the row count, and its HDF5 writer reverses dimension order
% on export - so a [num_samples, num_units] in-memory array is what ends up
% on disk as the schema-required [num_units, num_samples]. Building these
% as [num_units, num_samples] instead (the untransposed, more "natural"
% orientation) silently passes the height check only when num_samples
% happens to equal num_units, and errors otherwise.
dataDescription = 'Single unit activity';
[spike_times_vector, spike_times_index] = util.create_indexed_column( ...
  spikesConv.timesNWB', dataDescription);

nwb.units = types.core.Units( ...
  'colnames', { ... % Provide the column order. All column names have to be defined below
    'cluster_id','type', 'peak_channel_index','peak_channel_id','x','y','z', ...
    'area','lead_id','electrode_group','spike_times','spike_times_index'}, ...
  'description', 'Units table', ...
  'id', types.hdmf_common.ElementIdentifiers( ...
    'data', int64(0:numel(spikesConv.cluID) - 1)), ...
  'cluster_id', types.hdmf_common.VectorData( ...
    'data', cell2mat(spikesConv.cluID), ...
    'description', 'Unique cluster id'), ...
  'type', types.hdmf_common.VectorData( ...
    'data', spikesConv.labels, ...
    'description', 'Cluster type: unit vs mua'), ...
  'peak_channel_index', types.hdmf_common.VectorData( ...
    'data', cell2mat(spikesConv.channelInds), ...
    'description', 'Peak channel row index in the electrode table'), ...
  'peak_channel_id', types.hdmf_common.VectorData( ...
    'data', spikesConv.chLabels, ...
    'description', 'Unique ID of the channel with the largest cluster waveform amplitude'), ...
  'x', types.hdmf_common.VectorData( ...
    'data', cell2mat(spikesConv.x), ...
    'description', 'x coordinate'), ...
  'y', types.hdmf_common.VectorData( ...
    'data', cell2mat(spikesConv.y), ...
    'description', 'y coordinate'), ...
  'z', types.hdmf_common.VectorData( ...
    'data', cell2mat(spikesConv.z), ...
    'description', 'z coordinate'), ...
  'area', types.hdmf_common.VectorData( ...
    'data', spikesConv.areaLabels, ...
    'description', 'Brain area where the unit is located.'), ...
  'lead_id', types.hdmf_common.VectorData( ...
    'data', spikesConv.leadLabels, ...
    'description', 'Lead id where the unit is located'), ...
  'spike_times', spike_times_vector, ...
  'spike_times_index', spike_times_index, ...
  'electrode_group', types.hdmf_common.VectorData( ...
    'data', spikesConv.group, ...
    'description', 'Recording channel groups'), ...
  'waveform_mean', types.hdmf_common.VectorData( ...
    'data', cat(1, spikesConv.filtWaveform{:}).', ...
    'description', ['Mean waveforms on the probe channel with the largest waveform amplitude. ' ...
    'The order that waveforms are stored match the order of units in the unit table.']), ...
  'waveform_sd', types.hdmf_common.VectorData( ...
    'data', cat(1, spikesConv.filtWaveformSD{:}).', ...
    'description', 'Standard deviation of waveforms.'));

% Save the updated file
nwbExport(nwb, options.outputFile);

% matnwb's generic VectorData export path writes any array with a
% singleton dimension (e.g. a lone unit's [1, n_samples] waveform) as flat
% 1-D, violating the NWB schema's required [num_units, num_samples] shape
% for waveform_mean/waveform_sd whenever there is exactly one unit. Patch
% those two (tiny) datasets back to 2-D in place — see
% local_force_units_dim for the full explanation.
local_force_units_dim(options.outputFile, '/units/waveform_mean');
local_force_units_dim(options.outputFile, '/units/waveform_sd');



%% Helper functions
function n = getTimeDim(data, options)
% A helper function for extracting the total number of timeseries samples

arguments
  data
  options.whichDim (1,1) {isnumeric} = 2 % Columns are sample points (2)
end

if isa(data, 'types.untyped.DataStub')
  n = data.dims(options.whichDim);             % lazily-read dataset
elseif isa(data, 'types.untyped.DataPipe')
  internal = data.internal;                    % BlueprintPipe or BoundPipe
  if isprop(internal, 'data') && ~isempty(internal.data)
    n = size(internal.data, options.whichDim); % in-memory data not yet written
  else
    n = internal.dims(options.whichDim);    % already bound to file
  end
else
  n = size(data, options.whichDim);            % plain MATLAB array
end


function local_force_units_dim(nwbPath, dsPath)
% LOCAL_FORCE_UNITS_DIM  Ensure a Units waveform_mean/waveform_sd dataset
%   keeps its [num_units, num_samples] 2-D shape when there is only one
%   unit.
%
%   matnwb's generic VectorData export path (+io/mapData2H5.m) writes any
%   MATLAB array with a singleton dimension as a flat 1-D HDF5 dataset
%   unless its internal 'forceMatrix' flag is set — a flag that
%   types.core.Units' export() never passes for waveform_mean/waveform_sd
%   (unlike the dedicated ClusterWaveforms neurodata type, which does pass
%   it for its own, structurally identical, fields). A single-unit waveform
%   array is exactly [1, num_samples] — one dimension trivially 1 — so this
%   silently collapses to flat 1-D on export, violating the NWB schema
%   (core/nwb.misc.yaml: dims [num_units, num_samples]). This only ever
%   happens when num_units == 1: with more than one unit neither MATLAB
%   dimension is 1, so matnwb already writes proper 2-D data and this
%   function is a no-op.
%
%   Deleting and recreating the dataset to change its shape also drops its
%   existing HDF5 attributes (namespace/neurodata_type/object_id/unit/
%   description — all plain strings for a VectorData object), so those are
%   captured before the delete and rewritten after; nothing but the
%   dataspace shape actually changes.
if ~local_h5_link_exists(nwbPath, dsPath)
  return
end

info = h5info(nwbPath, dsPath);
if numel(info.Dataspace.Size) >= 2
  return % already >=2-D; nothing to do
end

data = reshape(h5read(nwbPath, dsPath), 1, []); % -> [num_units=1, num_samples]

attrNames = {info.Attributes.Name};
attrValues = cell(size(attrNames));
for iAttr = 1:numel(attrNames)
  attrValues{iAttr} = char(h5readatt(nwbPath, dsPath, attrNames{iAttr}));
end

fileID = H5F.open(nwbPath, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
H5L.delete(fileID, dsPath, 'H5P_DEFAULT');
H5F.close(fileID);

h5create(nwbPath, dsPath, size(data), 'Datatype', class(data));
h5write(nwbPath, dsPath, data);

fileID = H5F.open(nwbPath, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
datasetID = H5D.open(fileID, dsPath);
for iAttr = 1:numel(attrNames)
  local_write_vlen_string_attr(datasetID, attrNames{iAttr}, attrValues{iAttr});
end
H5D.close(datasetID);
H5F.close(fileID);


function tf = local_h5_link_exists(nwbPath, dsPath)
% LOCAL_H5_LINK_EXISTS  Check whether an HDF5 link exists without erroring.
fileID = H5F.open(nwbPath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
tf = H5L.exists(fileID, dsPath, 'H5P_DEFAULT');
H5F.close(fileID);


function local_write_vlen_string_attr(objID, name, value)
% LOCAL_WRITE_VLEN_STRING_ATTR  Write a scalar variable-length UTF-8 string
%   attribute (mirrors annotatedSz2nwb/inject_implantation_date_into_nwb.m's
%   helper of the same name, which confirmed this datatype/cset combination
%   matches matnwb's own attribute convention against a real file).
typeID = H5T.copy('H5T_C_S1');
H5T.set_size(typeID, 'H5T_VARIABLE');
H5T.set_cset(typeID, H5ML.get_constant_value('H5T_CSET_UTF8'));
spaceID = H5S.create('H5S_SCALAR');

attrID = H5A.create(objID, name, typeID, spaceID, 'H5P_DEFAULT');
H5A.write(attrID, typeID, {char(value)});

H5A.close(attrID);
H5S.close(spaceID);
H5T.close(typeID);