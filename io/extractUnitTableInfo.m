function spikeData = extractUnitTableInfo(binaryFileBasename, unitTableFile)
% extractUnitTableInfo(binaryFileBasename, unitTableFile)
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
%
% Each chunk's start time is read from its own companion conversion .mat
% file (saved by nwb2binary.m alongside the chunk's binary data), which
% is required to exist - this replaces the previous chunkDuration-based
% assumption that every chunk spans the same duration, which silently
% produced wrong start times whenever chunk folders span different
% numbers of raw segments.
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

  % Locate the trailing 'Total:' summary row (first column reads 'Total:')
  % and treat every row from it downwards as non-data. Previously the code
  % simply dropped the final table row, which only worked while 'Total:'
  % happened to be the last row of Excel's used range - any blank row or
  % stray cell below it left the Total row in place, and its per-chunk unit
  % COUNTS were then misread as Phy cluster IDs.
  tableHeight = height(unitTable);
  firstCol = unitTable.(labels{1});
  if ~iscell(firstCol)
    firstCol = cellstr(string(firstCol));
  end
  totalRowIdx = find(~cellfun('isempty', regexpi(strtrim(firstCol), '^total', 'once')), 1);
  if isempty(totalRowIdx)
    nDataRows = tableHeight - 1; % preserve previous behaviour: drop the last row
  else
    nDataRows = totalRowIdx - 1;
  end

  for iColumn = 1:nColumns
    label = labels{iColumn};
    columnData = unitTable.(label)';
    columnData = columnData(1:nDataRows);
    if iscell(columnData)
      % A column reads back as a cell array of char instead of numeric
      % double whenever Excel/readtable sees mixed cell formatting within
      % that column (e.g. some numeric-looking values entered/pasted as
      % text) - even a single such cell makes readtable treat the whole
      % column this way. Blank cells then come back as '' rather than
      % NaN, so they must be masked out explicitly here instead of being
      % (incorrectly) treated as valid entries.
      isEntryMissing = cellfun(@(x) isempty(x) || (ischar(x) && isempty(strtrim(x))), columnData);
      numericColumnData = nan(1, numel(columnData));
      numericColumnData(~isEntryMissing) = str2double(columnData(~isEntryMissing));
      columnData = numericColumnData;
      valueMask = ~isnan(columnData);
    elseif isnumeric(columnData)
      valueMask = ~isnan(columnData);
    else
      valueMask = true(1,numel(columnData));
    end
    nUnits = sum(valueMask);
    columnData = columnData(valueMask);
    valueMask = [valueMask false(1, tableHeight - nDataRows)]; %#ok<AGROW>
    if (startsWith(label, 'x') || startsWith(label, '_')) && (endsWith(label, 'Id') || endsWith(label, 'id'))
      nBinFiles = nBinFiles + 1;
      spikeData.existingUnitIDs = [spikeData.existingUnitIDs columnData];
      spikeData.newGlobalUnitIDs = [spikeData.newGlobalUnitIDs unitTable.id(valueMask)'];
      spikeData.newGlobalUnitCh = [spikeData.newGlobalUnitCh unitTable.ch(valueMask)'];
      spikeData.chLabels = [spikeData.chLabels unitTable.channelLabel(valueMask)'];
      spikeData.leadLabels = [spikeData.leadLabels unitTable.leadLabel(valueMask)'];
      spikeData.areaLabels = [spikeData.areaLabels unitTable.areaLabel(valueMask)'];
      try
        chunkFolder = [binaryFileBasename label(2:8)];
        matFilename = [chunkFolder filesep 'temp_wh.spikes.cellinfo.mat'];
      catch
        chunkFolder = binaryFileBasename;
        matFilename = [binaryFileBasename filesep 'temp_wh.spikes.cellinfo.mat'];
      end

      % Determine this chunk's start time from nwb2binary.m's own
      % conversion record (the chunk folder's companion .mat file) rather
      % than assuming every chunk spans the same chunkDuration. Segment-
      % range folders are not guaranteed to be equal width (e.g. a
      % recording's chunks may span 40/20/20/20/7 raw segments
      % respectively), in which case a single chunkDuration cannot give
      % correct start times for every chunk. There is no safe fallback
      % for a missing conversion file - guessing via chunkDuration is
      % exactly the assumption that silently breaks on non-uniform chunks
      % - so this errors instead.
      [~, chunkFolderName] = fileparts(chunkFolder);
      conversionMatFile = fullfile(chunkFolder, [chunkFolderName '.mat']);
      if ~isfile(conversionMatFile)
        error('extractUnitTableInfo:missingConversionFile', ...
          ['No companion conversion .mat file found for chunk folder %s ' ...
          '(expected %s). This file is produced by nwb2binary.m and is ' ...
          'required to determine the chunk''s true start time.'], ...
          chunkFolder, conversionMatFile);
      end
      conversionInfo = load(conversionMatFile, 'sessionStartTime'); %#ok<*LOAD>
      chunkStartTime = conversionInfo.sessionStartTime;

      for iUnit = 1:nUnits
        spikeData.files = [spikeData.files matFilename];
        spikeData.startTimes = [spikeData.startTimes chunkStartTime];
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
          spikeData_muas.startTimes = [spikeData_muas.startTimes chunkStartTime];
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