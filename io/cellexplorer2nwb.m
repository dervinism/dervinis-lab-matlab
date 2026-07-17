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
if ~isempty(timeseriesData.timestamps)
  timestamps = timeseriesData.timestamps.load();
else
  t0 = timeseriesData.starting_time;        % seconds, relative to 12AM of the session_start_time day
  rate = timeseriesData.starting_time_rate; % Hz (samples per second)
  n = getTimeDim(timeseriesData.data);      % number of samples along the TIME dimension
  timestamps = t0 + (0:double(n)-1)'/rate;  % column vector of timestamps in seconds
end

% Remap spike times based on NWB timestamps because kilosort spike sorter
% uses fixed sampling rate which might not always agree with NWB timestamps
nUnits = numel(spikesConv.cluID);
spikesConv.timesNWB = cell(1,nUnits);
if ~isempty(timestamps)
  for iUnit = 1:nUnits
    spikeInds = round(spikesConv.times{iUnit}./(1/spikes.sr));
    spikesConv.timesNWB{iUnit} = timestamps(spikeInds);
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
    'data', cell2mat(cellfun(@(w) w(:), spikesConv.filtWaveform, 'UniformOutput', false)), ...
    'description', ['Mean waveforms on the probe channel with the largest waveform amplitude. ' ...
    'The order that waveforms are stored match the order of units in the unit table.']), ...
  'waveform_sd', types.hdmf_common.VectorData( ...
    'data', cell2mat(cellfun(@(w) w(:), spikesConv.filtWaveformSD, 'UniformOutput', false)), ...
    'description', 'Standard deviation of waveforms.'));

% Save the updated file
nwbExport(nwb, options.outputFile);



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