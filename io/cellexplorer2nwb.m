function unitsData = cellexplorer2nwb(spikeData, nwbFile, options)
% unitsData = cellexplorer2nwb(spikeData, nwbFile, <options>)
%
% Function converts spike sorting output saved in the CellExplorer
% (https://cellexplorer.org/) format (spikes.cellinfo.mat) to NWB format by
% attaching spike times timeseries data to an exisiting raw NWB file.
%
% Writes /units directly via low-level HDF5 (H5F/H5G/H5T/H5S/H5D/H5A/H5R -
% see writeUnitsTableH5.m), in place, with no nwbRead/nwbExport and no
% MatNWB dependency at write time. The previous implementation went through
% nwbRead()+nwbExport(), MatNWB's general-purpose (and, for this table,
% comparatively slow) export path - the one that produced two bugs this
% repo has already had to patch after the fact: 'spike_times_index' leaking
% into /units/colnames (breaks pynwb reads entirely - see
% annotatedSz2nwb/fix_units_colnames.m) and a single-unit waveform_mean/
% waveform_sd collapsing to a 1-D dataset instead of the schema-required
% 2-D shape (see local_force_units_dim, formerly in this file). The
% low-level writer gets both right from the start instead of patching
% around them, and was verified byte-for-byte (including the electrode_group
% and spike_times_index HDF5 object references, and the waveform matrices'
% dataspace orientation) against a real MatNWB-written /units group before
% this file was changed - see io/test_writeUnitsTableH5.m.
% resolveTimestampsAtIndicesH5.m replaces resolveTimestampsAtIndices.m for
% the same reason (no matnwb TimeSeries object required); its own chunked-
% read algorithm is unchanged, just re-pointed at nwbPath+acqKey instead of
% a matnwb DataStub.
%
% If options.outputFile differs from nwbFile, nwbFile is copied there first
% (once) and the copy is edited in place - matching every other injector in
% annotatedSz2nwb (e.g. inject_seizures_epic_into_nwb.m). The default
% (outputFile unset) edits nwbFile in place, no copy.
%
% Args:
%   spikeData (struct, required, positional): a shape-(1,1) scalar
%     structure with the following fields:
%     files - a cell array of CellExplorer files (spikes.cellinfo.mat)
%       containing spike times of the indicated units in the field
%       existingUnitIDs.
%     startTimes - spike time displacements relative to the start of the
%       full recording file (informational only; not used to place spikes).
%     startIndices - 1-based sample indices, into the full NWB acquisition
%       TimeSeries, of the first sample of each unit's binary chunk (see
%       extractUnitTableInfo.m). Each spike's NWB sample index is its
%       0-based Kilosort sample within the chunk plus this start index; its
%       time is then the NWB timestamp at that index. Spikes must NOT be
%       placed via round((localTime + startTime)*sr): startTime is on the
%       NWB timestamp clock, which drifts relative to the nominal sampling
%       rate, so that conversion lands later chunks' spikes the wrong number
%       of samples away (e.g. 2661 samples / 83 ms late in chunk 2 of
%       sub-MSEL02750 ses-nonstim01).
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
%   dryRun (logical, optional, keyword): a shape-(1, 1) logical scalar. If
%     true, everything is computed from nwbFile but nothing is copied or
%     written - use the returned unitsData to verify spike times against
%     the raw data before injecting (default is false).
%   acqGroup (char, optional, keyword): name of the NWB acquisition
%     TimeSeries the spikes were sorted from. Default '' derives it from
%     the CellExplorer sampling rate ('TimeSeries_<sr>_Hz'), which is wrong
%     when the group was renamed after sorting - e.g. recordings whose MED
%     header rate was corrected (MSEL02685 ses-nonstim01: sorted as 14223 Hz,
%     now TimeSeries_32000_Hz). Sample indices are unaffected.
%
% Returns:
%   unitsData (struct): the input struct passed to writeUnitsTableH5.m,
%     including spikeTimes (NWB timestamps) and spikeIndices (1-based NWB
%     sample indices) per unit.
%
% Dependencies:
%   writeUnitsTableH5.m, resolveTimestampsAtIndicesH5.m (this repo's io
%     folder). MatNWB is NOT required.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeData (1,1) {mustBeA(spikeData,'struct')}
  nwbFile (1,:) {mustBeA(nwbFile,'char'),mustBeVector,endsWith(nwbFile,'.nwb')}
  options.outputFile (1,:) {mustBeA(options.outputFile,'char'),mustBeVector,endsWith(options.outputFile,'.nwb')} = '';
  options.loadCoordsFromUnitTable (1,1) {islogical} = false
  options.dryRun (1,1) {islogical} = false
  options.acqGroup (1,:) char = ''
end

% Parse input
if ~isfield(spikeData, 'startIndices') || numel(spikeData.startIndices) ~= numel(spikeData.existingUnitIDs)
  error('cellexplorer2nwb:missingStartIndices', ...
    ['spikeData.startIndices (one 1-based chunk start sample index per unit) ' ...
    'is required - regenerate spikeData with the current extractUnitTableInfo.m.']);
end
if options.dryRun
  options.outputFile = nwbFile; % read-only use below; nothing is written
elseif isempty(options.outputFile)
  options.outputFile = nwbFile;
elseif ~strcmp(options.outputFile, nwbFile)
  copyfile(nwbFile, options.outputFile);
end

% Load spiking data and convert it
nUnits = numel(spikeData.existingUnitIDs);
if ~nUnits
  warning('Empty spikeData supplied. Exiting cellexplorer2nwb function.');
  return
end
loadedFile = '';
loadedType = '';
spikesConv = struct();
for iUnit = 1:nUnits
  fileToLoad = spikeData.files{iUnit};
  typeToLoad = spikeData.types{iUnit};
  % A "streak" is (file, type): extractUnitTableInfo.m lays spikeData out
  % as [every chunk's manually-curated 'unit' rows, in chunk order] then
  % [every chunk's auto-swept 'mua' rows, in chunk order] - so a given
  % chunk's own 'unit' rows and its own 'mua' rows can land back-to-back in
  % this array (same fileToLoad, no other chunk's rows in between)
  % whenever every later chunk in the 'unit' pass contributed zero rows
  % (e.g. a chunk with only MUA, no accepted single units - a legitimate,
  % not-uncommon case). Keying the streak boundary on file alone conflated
  % that chunk's 'unit' count and 'mua' count into one combined total that
  % matches neither individually, raising a false "not all units checked"
  % error even though every row was in fact accounted for. Keying on
  % (file, type) keeps each streak homogeneous so the check below compares
  % against the single correct expected count.
  if ~strcmpi(loadedFile,fileToLoad) || ~strcmpi(loadedType,typeToLoad)
    if exist('trackCheckedUnits','var')
      expectedCount = numel(spikes.cluID(cellfun(@(x) strcmpi(x, loadedType), spikes.labels)));
      if numel(trackCheckedUnits) ~= expectedCount
        error(['Not all ' loadedType ' units from the previous CellExplorer file were checked: ' ...
          loadedFile]);
      end
    end
    if ~strcmpi(loadedFile,fileToLoad)
      load(fileToLoad); %#ok<*LOAD>
      loadedFile = fileToLoad;
    end
    loadedType = typeToLoad;
    trackCheckedUnits = [];
  end
  unitIndExisting = find(spikeData.existingUnitIDs(iUnit) == spikes.cluID);
  trackCheckedUnits = [trackCheckedUnits spikes.cluID(unitIndExisting)]; %#ok<AGROW>
  if isempty(unitIndExisting)
    error(['Unit ' num2str(spikeData.existingUnitIDs(iUnit)) ...
      ' does not exist in the file ' fileToLoad]);
  end
  unitIndNew = spikeData.newGlobalUnitIDs(iUnit);
  % Kilosort spike times are 0-based chunk sample numbers divided by the
  % nominal rate, so round(times*sr) recovers the exact sample number.
  chunkSpikeInds = round(spikes.times{unitIndExisting}.*spikes.sr) + spikeData.startIndices(iUnit);
  if iUnit == 1 || unitIndNew > numel(spikesConv.cluID)
    spikesConv.inds{unitIndNew} = chunkSpikeInds;
    spikesConv.cluID{unitIndNew} = unitIndNew;
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
    spikesConv.inds{unitIndNew} = [spikesConv.inds{unitIndNew}; chunkSpikeInds];
    spikesConv.cluID{unitIndNew} = unitIndNew;
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

% Load electrode table + acquisition sampling rate directly via low-level
% HDF5 (no matnwb) - only the columns actually needed.
electrodesTable = local_read_electrodes_table(options.outputFile, ...
  {'ChName','label','x','y','z','group_name'});
if isempty(options.acqGroup)
  acqKey = ['TimeSeries_' num2str(spikes.sr) '_Hz'];
else
  acqKey = options.acqGroup;
end

% Remap spike sample indices to NWB timestamps because kilosort spike sorter
% uses fixed sampling rate which might not always agree with NWB timestamps.
% Only the specific per-sample timestamps a unit's spikes actually land on
% are fetched/computed - materializing the full per-sample timestamps
% vector for the whole session (previously done unconditionally here) ran
% out of memory on long, high sampling-rate recordings.
nUnits = numel(spikesConv.cluID);
spikesConv.timesNWB = cell(1,nUnits);
% Sort each unit's spikes: they are concatenated chunk by chunk in unit-table
% column order, which is not guaranteed to be chronological, and NWB spike
% times must be ascending.
spikeIndsPerUnit = cellfun(@sort, spikesConv.inds, 'UniformOutput', false);
allSpikeInds = cat(1, spikeIndsPerUnit{:});
if ~isempty(allSpikeInds)
  allTimestamps = resolveTimestampsAtIndicesH5(options.outputFile, acqKey, allSpikeInds);
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
  % Store the file's own electrode label, so peak_channel_id always matches
  % the electrodes table (unit tables may spell it differently, e.g.
  % 'LBmicro7' for an electrode labelled 'LB_micro7').
  spikesConv.chLabels{iUnit} = electrodesTable.label{channelInd};
  if ~options.loadCoordsFromUnitTable
    spikesConv.x{iUnit} = electrodesTable.x(spikesConv.channelInds{iUnit});
    spikesConv.y{iUnit} = electrodesTable.y(spikesConv.channelInds{iUnit});
    spikesConv.z{iUnit} = electrodesTable.z(spikesConv.channelInds{iUnit});
  end
  spikesConv.groupName{iUnit} = electrodesTable.group_name{spikesConv.channelInds{iUnit}};
end

% Mean and SD waveforms on each unit's peak channel, from the NWB raw data
% in volts (see computeUnitWaveformsH5.m). The acquisition group's own
% electrodes dataset maps its data columns to (0-based) electrode rows.
acqElectrodeRows = double(h5read(options.outputFile, ['/acquisition/' acqKey '/electrodes'])) + 1;
acqColumns = zeros(1, nUnits);
for iUnit = 1:nUnits
  column = find(acqElectrodeRows == spikesConv.channelInds{iUnit});
  assert(isscalar(column), 'cellexplorer2nwb:peakChannelNotInAcquisition', ...
    'Peak channel %s of unit %d is not recorded in %s.', spikesConv.chLabels{iUnit}, iUnit, acqKey);
  acqColumns(iUnit) = column;
end
[waveformMean, waveformSd] = computeUnitWaveformsH5(options.outputFile, acqKey, ...
  spikeIndsPerUnit, acqColumns, sampleRate=spikes.sr);

% Assemble writeUnitsTableH5's input struct and write /units in place.
unitsData = struct();
unitsData.clusterId          = cell2mat(spikesConv.cluID)';
unitsData.type                = spikesConv.labels(:);
unitsData.peakChannelIndex    = cell2mat(spikesConv.channelInds)' - 1; % 0-based electrodes table row (NWB convention)
unitsData.peakChannelId       = spikesConv.chLabels(:);
unitsData.x                   = cell2mat(spikesConv.x)';
unitsData.y                   = cell2mat(spikesConv.y)';
unitsData.z                   = cell2mat(spikesConv.z)';
unitsData.area                = spikesConv.areaLabels(:);
unitsData.leadId               = spikesConv.leadLabels(:);
unitsData.electrodeGroupName  = spikesConv.groupName(:);
unitsData.spikeTimes          = spikesConv.timesNWB;
unitsData.waveformMean        = waveformMean.';
unitsData.waveformSd          = waveformSd.';

if ~options.dryRun
  writeUnitsTableH5(options.outputFile, unitsData);
end
unitsData.spikeIndices = spikeIndsPerUnit;


%% Helper functions
function electrodesTable = local_read_electrodes_table(nwbPath, columns)
% LOCAL_READ_ELECTRODES_TABLE  Read select /general/extracellular_ephys/
%   electrodes columns via h5read, without matnwb. Duplicated from (rather
%   than calling) annotatedSz2nwb/read_nwb_electrodes_table.m to keep this
%   library free of a dependency on that downstream project - dervinis-lab-
%   matlab is the lower-level shared library annotatedSz2nwb depends on,
%   never the reverse.
groupPath = '/general/extracellular_ephys/electrodes';
data = struct();
for i = 1:numel(columns)
  raw = h5read(nwbPath, [groupPath '/' columns{i}]);
  if iscell(raw)
    raw = cellfun(@(x) strtrim(char(x)), raw, 'UniformOutput', false);
  elseif ischar(raw) || isstring(raw)
    raw = cellstr(raw);
  end
  data.(columns{i}) = raw(:);
end
electrodesTable = data;
