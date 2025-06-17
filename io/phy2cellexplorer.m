function phy2cellexplorer(phyFolder, binaryFile, outputFolder, options)
% phy2cellexplorer(kilosortFolder, binaryFile, outputFolder, options)
%
% Function converts Phy (https://github.com/cortex-lab/phy) manual curation
% output to a format that is compatible with and can be read by
% CellExplorer (https://cellexplorer.org/). As a result the following files
% are saved in the same folder as the Phy output:
% noiseLevel.channelInfo.mat, session.mat, spikes.cellinfo.mat.
%
% Args:
%   phyFolder (char, required, positional): a shape-(1, n) character array
%     containing the full path to the Phy output folder.
%   binaryFile (char, required, positional): a shape-(1, m) character array
%     containing the full path to the binary file that was used for spike
%     sorting.
%   outputFolder (char, required, positional): a shape-(1, l) character
%     array containing the full path to the folder where files for
%     CellExplorer are saved. Typically, it should be the same folder as
%     phyFolder.
%   nChannels (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing the number of data recording channels (default=16).
%   samplingRate (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing the data sampling rate (default=32000).
%   precision (char, optional, keyword): a shape-(1, k) character array
%     describing the format that the binary file is saved in
%     (default='int16'). CellExplorer only reads 'int16' format binary
%     files, so make sure that the spikesorting prior to Phy examination is
%     performed on 'int16' binary files. Available formats include 'int8',
%     'uint8', 'int16', 'uint16', 'int32', 'uint32', 'single', 'int64',
%     'uint64', and 'double'.
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  phyFolder (1,:) {mustBeA(phyFolder,'char'),mustBeVector}
  binaryFile (1,:) {mustBeA(binaryFile,'char'),mustBeVector}
  outputFolder (1,:) {mustBeA(outputFolder,'char'),mustBeVector}
  options.nChannels (1,1) {mustBePositive} = 16
  options.samplingRate (1,1) {mustBePositive} = 32000
  options.precision  (1,:) {mustBeA(options.precision,'char'),mustBeVector} = 'int16'
end

cd(phyFolder);

% Convert to CellExplorer path naming convention
basepath = phyFolder;
[~, basename] = fileparts(binaryFile);

% Create the session file and save it
session = sessionTemplate(basepath);
session.general.name = basename;
session.general.basepath = outputFolder;

session.animal.name = 'Unknown';
session.animal.species = 'human';

session.extracellular.sr = options.samplingRate;
session.extracellular.srLfp = options.samplingRate;
session.extracellular.nChannels = options.nChannels;
session.extracellular.fileName = binaryFile;
session.extracellular.precision = options.precision;

save(fullfile(session.general.basepath, ...
  [session.general.name, '.session.mat']), 'session');

% Run the Phy to CellExplorer pipeline
loadSpikes('basepath',basepath, 'basename',basename, ...
  'labelsToRead',{'good','mua'}, 'forceReload',true, 'showWaveforms',false);

% Move pipeline output files to a dedicated folder
if ~strcmpi(phyFolder, outputFolder)
  source = fullfile(phyFolder, [basename '.noiseLevel.channelInfo.mat']);
  destination = fullfile(outputFolder, [basename '.noiseLevel.channelInfo.mat']);
  movefile(source, destination);
  source = fullfile(phyFolder, [basename '.spikes.cellinfo.mat']);
  destination = fullfile(outputFolder, [basename '.spikes.cellinfo.mat']);
  movefile(source, destination);
end

% Delete any remainin output
trash = fullfile(phyFolder, 'summaryFigures');
try
  rmdir(trash, 's');
catch
  % Do nothing
end

cd(outputFolder);