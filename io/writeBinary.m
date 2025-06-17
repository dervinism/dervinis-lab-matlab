function outputFilename = writeBinary(data, fileName, options)
% writeBinary(data, fileName, <writeBufferSize>)
%
% Function saves a row matrix data as a binary file (.dat).
%
% Args:
%   data (numeric, required, positional): a shape-(m, n) numeric array of
%     data.
%   fileName (char, required, positional): a shape-(1, n) character array
%     containing the output file name with a full path. It must end with
%     '.dat'.
%   writeBufferSize (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing the size of the binary write buffer in bytes.
%     To improve the processing speed, it should be a multiple of 4096 as
%     this is the size of most hard drive buffers. Default is 500*4096
%     bytes.
%
% Output: outputFilename is the name of the binary file when saved.

arguments
  data (:,:) {mustBeNumeric}
  fileName (1,:) {mustBeA(fileName,'char'),mustBeVector,endsWith(fileName,'.dat')}
  options.writeBufferSize (1,1) {mustBePositive} = 500*4096
end

fidOut = [];

w = whos('data');
format = w.class;
if strcmp(format, 'int8') || strcmp(format, 'uint8')
  nSampsTotal = format;
elseif strcmp(format, 'int16') || strcmp(format, 'uint16')
  nSampsTotal = w.bytes/2;
elseif strcmp(format, 'int32') || strcmp(format, 'uint32') || strcmp(format, 'single')
  nSampsTotal = w.bytes/4;
elseif strcmp(format, 'int64') || strcmp(format, 'uint64') || strcmp(format, 'double')
  nSampsTotal = w.bytes/8;
else
  error(['Unsupported data format. '...
    'Supported formats are double, single, int8, int16, int32, int64, uint8, uint16, uint32, uint64']);
end
  
nChunksTotal = ceil(nSampsTotal/options.writeBufferSize);

try
  outputFilename  = fileName;
  fidOut = fopen(outputFilename, 'w');
  
  chunkInd = 1;
  while 1
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    inds = (1:options.writeBufferSize) + (chunkInd-1)*options.writeBufferSize;
    if inds(1) > numel(data)
      break
    elseif inds(end) > numel(data)
      inds = inds(1):numel(data);
    end
    dat = data(inds);
    fwrite(fidOut, dat, format);
    chunkInd = chunkInd+1;
  end
  
  fclose(fidOut);
  
catch me
  if ~isempty(fidOut)
    fclose(fidOut);
  end
  
  rethrow(me)
  
end