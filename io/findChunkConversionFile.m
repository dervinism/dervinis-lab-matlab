function conversionMatFile = findChunkConversionFile(chunkFolder)
% conversionMatFile = findChunkConversionFile(chunkFolder)
%
% Function returns the path to a chunk folder's companion conversion .mat
% file (saved by nwb2binary.m), or '' if there is none.
%
% nwb2binary.m names the file after its output binary with unpadded
% segment numbers (e.g. '..._segments_1_20.mat'), while chunk folders are
% often renamed with zero-padded ones ('..._segments_001_020') so they
% sort correctly. Both the folder's own name and the unpadded variant are
% therefore accepted.
%
% Args:
%   chunkFolder (char, required, positional): a shape-(1, n) character
%     array containing the path to the chunk folder.
%
% Returns:
%   conversionMatFile (char): a shape-(1, m) character array containing
%     the path to the companion .mat file, or '' if not found.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  chunkFolder (1,:) {mustBeA(chunkFolder,'char'),mustBeVector}
end

[~, chunkFolderName] = fileparts(chunkFolder);
candidates = {[chunkFolderName '.mat']};
tok = regexp(chunkFolderName, '^(.*_segments_)(\d+)_(\d+)$', 'tokens', 'once');
if ~isempty(tok)
  candidates{end+1} = sprintf('%s%d_%d.mat', tok{1}, str2double(tok{2}), str2double(tok{3}));
end

conversionMatFile = '';
for iCandidate = 1:numel(candidates)
  candidate = fullfile(chunkFolder, candidates{iCandidate});
  if isfile(candidate)
    conversionMatFile = candidate;
    return
  end
end
end
