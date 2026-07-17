function spikeData = extractUnitTableInfo(binaryFileBasename, unitTableFile, options)
% extractUnitTableInfo(binaryFileBasename, unitTableFile, <chunkDuration>)
%
% Function extracts some useful info about recorded units from a custom
% Excel file.
%
% Args:
%   binaryFileBasename (char, required, positional): a shape-(1, n)
%     character array containing basename of CellExplorer spikes.cellinfo
%     files.
%   unitTableFile (char, required, positional): a shape-(1, m) character
%     array containing the path to the unit table Excel file.
%   chunkDuration (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar containing the data chunk duration in seconds corresponding to
%     a single spike sorted binary file. Default is 3125 seconds.
%
% Returns:
%   spikeData (struct): a shape-(1, 1) Matlab scalar structure with the
%     following fields:
%       existingUnitIDs - old unit IDs from Phy;
%       newGlobalUnitIDs - newly-assigned global unit IDs;
%       newGlobalUnitCh - unit channel number;
%       chLabels - unit channel labels;
%       leadLabels - lead labels;
%       areaLabels - area where units were recorded;
%       files - names of binary files with recorded units;
%       startTimes - start times of those binary files in seconds.
%       x, y, z - 3D atlas coordinates
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  binaryFileBasename (1,:) {mustBeA(binaryFileBasename,'char'),mustBeVector}
  unitTableFile (1,:) {mustBeA(unitTableFile,'char'),mustBeVector}
  options.chunkDuration (1,1) {mustBePositive} = 3125
end

% Load the unit table
if isempty(unitTableFile)
  unitTable = [];
else
  unitTable = readtable(unitTableFile);
end

% Initialise storage variables
spikeData.existingUnitIDs = [];
spikeData.newGlobalUnitIDs = [];
spikeData.newGlobalUnitCh = [];
spikeData.chLabels = {};
spikeData.leadLabels = {};
spikeData.areaLabels = {};
spikeData.files = {};
spikeData.startTimes = [];
spikeData.x = [];
spikeData.y = [];
spikeData.z = [];
spikeData.types = {};

spikeData_muas.existingUnitIDs = [];
spikeData_muas.newGlobalUnitCh = [];
spikeData_muas.newGlobalUnitIDs = [];
spikeData_muas.files = {};
spikeData_muas.startTimes = [];
spikeData_muas.types = {};

% Extract data
if isempty(unitTable)
  matFilename = [binaryFileBasename filesep 'temp_wh.spikes.cellinfo.mat'];
  load(matFilename) %#ok<*LOAD>
  spikeData.existingUnitIDs = spikes.cluID;
  spikeData.newGlobalUnitCh = spikes.maxWaveformCh;
  for iUnit = 1:spikes.numcells
    spikeData.newGlobalUnitIDs = [spikeData.newGlobalUnitIDs iUnit];
    spikeData.files = [spikeData.files matFilename];
    spikeData.startTimes = [spikeData.startTimes 0];
  end
  spikeData.types = spikes.labels;
else
  nBinFiles = 0;
  nColumns = size(unitTable,2);
  labels = unitTable.Properties.VariableNames;
  for iColumn = 1:nColumns
    label = labels{iColumn};
    columnData = unitTable.(label)';
    columnData = columnData(1:end-1);
    if isnumeric(columnData)
      valueMask = ~isnan(columnData);
    else
      valueMask = true(1,numel(columnData));
    end
    nUnits = sum(valueMask);
    columnData = columnData(valueMask);
    valueMask = [valueMask false]; %#ok<AGROW>
    if (startsWith(label, 'x') || startsWith(label, '_')) && (endsWith(label, 'Id') || endsWith(label, 'id'))
      nBinFiles = nBinFiles + 1;
      spikeData.existingUnitIDs = [spikeData.existingUnitIDs columnData];
      spikeData.newGlobalUnitIDs = [spikeData.newGlobalUnitIDs unitTable.id(valueMask)'];
      spikeData.newGlobalUnitCh = [spikeData.newGlobalUnitCh unitTable.ch(valueMask)'];
      spikeData.chLabels = [spikeData.chLabels unitTable.channelLabel(valueMask)'];
      spikeData.leadLabels = [spikeData.leadLabels unitTable.leadLabel(valueMask)'];
      spikeData.areaLabels = [spikeData.areaLabels unitTable.areaLabel(valueMask)'];
      try
        matFilename = [binaryFileBasename label(2:8) filesep 'temp_wh.spikes.cellinfo.mat'];
      catch
        matFilename = [binaryFileBasename filesep 'temp_wh.spikes.cellinfo.mat'];
      end
      for iUnit = 1:nUnits
        spikeData.files = [spikeData.files matFilename];
        spikeData.startTimes = [spikeData.startTimes (nBinFiles-1)*options.chunkDuration];
        spikeData.types = [spikeData.types 'unit'];
      end
      try
        spikeData.x = [spikeData.x unitTable.x(valueMask)'];
        spikeData.y = [spikeData.y unitTable.y(valueMask)'];
        spikeData.z = [spikeData.z unitTable.z(valueMask)'];
      catch
        % do nothing
      end

      load(matFilename) %#ok<*LOAD>
      for iUnit = 1:spikes.numcells
        if strcmpi(spikes.labels{iUnit}, 'mua')
          spikeData_muas.existingUnitIDs = [spikeData_muas.existingUnitIDs spikes.cluID(iUnit)];
          spikeData_muas.newGlobalUnitIDs = [spikeData_muas.newGlobalUnitIDs numel(spikeData_muas.newGlobalUnitIDs)+1];
          spikeData_muas.newGlobalUnitCh = [spikeData_muas.newGlobalUnitCh spikes.maxWaveformCh(iUnit)];
          spikeData_muas.files = [spikeData_muas.files matFilename];
          spikeData_muas.startTimes = [spikeData_muas.startTimes (nBinFiles-1)*options.chunkDuration];
          spikeData_muas.types = [spikeData_muas.types 'mua'];
        end
      end
    end
  end

  if ~isempty(spikeData_muas.existingUnitIDs)
    nUnits = max(spikeData.newGlobalUnitIDs);
    nEntries = numel(spikeData.newGlobalUnitIDs);
    for iUnit = 1:numel(spikeData_muas.existingUnitIDs)
      spikeData.existingUnitIDs = [spikeData.existingUnitIDs spikeData_muas.existingUnitIDs(iUnit)];
      spikeData.newGlobalUnitIDs = [spikeData.newGlobalUnitIDs spikeData_muas.newGlobalUnitIDs(iUnit)+nUnits];
      spikeData.newGlobalUnitCh = [spikeData.newGlobalUnitCh spikeData_muas.newGlobalUnitCh(iUnit)];
      spikeData.files{nEntries+iUnit} = spikeData_muas.files{iUnit};
      spikeData.startTimes = [spikeData.startTimes spikeData_muas.startTimes(iUnit)];
      spikeData.types{nEntries+iUnit} = spikeData_muas.types{iUnit};
      chInds = find(ismember(spikeData.newGlobalUnitCh(1:end-1), spikeData.newGlobalUnitCh(end)));
      if ~isempty(chInds)
        chInd = chInds(1);
        spikeData.chLabels{nEntries+iUnit} = spikeData.chLabels{chInd};
        spikeData.leadLabels{nEntries+iUnit} = spikeData.leadLabels{chInd};
        spikeData.areaLabels{nEntries+iUnit} = spikeData.areaLabels{chInd};
        if ~isempty(spikeData.x)
          spikeData.x = [spikeData.x spikeData.x(chInd)];
          spikeData.y = [spikeData.y spikeData.y(chInd)];
          spikeData.z = [spikeData.z spikeData.z(chInd)];
        end
      else
        spikeData.chLabels{nEntries+iUnit} = 'unknown';
        spikeData.leadLabels{nEntries+iUnit} = 'unknown';
        spikeData.areaLabels{nEntries+iUnit} = 'unknown';
        if ~isempty(spikeData.x)
          spikeData.x = [spikeData.x 0];
          spikeData.y = [spikeData.y 0];
          spikeData.z = [spikeData.z 0];
        end
      end
    end
  end
end