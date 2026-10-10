function [spikeTimes, somethingelse] = nwbSpikes2mat(inputFile, options)
% [spikeTimes, somethingelse] = nwbSpikes2mat(inputFile, <options>)
%
% Function loads spiking data from an NWB file and outputs it as a cell
% array of spike times.
%
% Args:
%   inputFile (char, required, positional): a shape-(1, n) character
%     array containing the input NWB file name with a full path. It must
%     end with '.nwb'.
%   verbose (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling the verbosity level of the function (default=false).
%
% Returns:
%   spikeTimes (cell): a shape-(m, 1) cell array of spike times with
%     individual elements corresponding to units.
%   somethingelse (numeric): a shape-(1, h) numeric array with timestamps
%     corresponding to the columns of timeseriesData.
%
% Dependencies:
%   matnwb (https://neurodatawithoutborders.github.io/matnwb/)
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  inputFile (1,:) {mustBeA(inputFile,'char'),mustBeVector,endsWith(inputFile,'.nwb')}
  options.verbose (1,1) {islogical} = false
end

% Load the nwb file
if options.verbose
  disp('Loading NWB data')
end
nwbData = nwbRead(inputFile);

% Access the units table
units = nwbData.units;

% Get spike times for all units
spikeTimes_vector = units.spike_times.data.load(); % all times concatenated
spikeTimes_index  = units.spike_times_index.data.load(); % end indices per unit

% Reconstruct per-unit spike times
num_units = length(spikeTimes_index);
spikeTimes = cell(num_units, 1);

prev_idx = 0;
for i = 1:num_units
    end_idx = spikeTimes_index(i);
    spikeTimes{i} = spikeTimes_vector(prev_idx+1 : end_idx);
    prev_idx = end_idx;
end

% Plot raster for unit 1
figure;
plot(spikeTimes{1}, ones(size(spikeTimes{1})), '|k', 'MarkerSize', 10);
xlabel('Time (s)'); title('Unit 1 Spike Raster');