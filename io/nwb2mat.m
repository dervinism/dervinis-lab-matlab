function [timeseriesData, timestamps] = nwb2mat(inputFile, options)
% [timeseriesData, timestamps] = nwb2mat(inputFile, <options>)
%
% Function loads a file in NWB format and outputs it as a row matrix with
% rows corresponding to individual recording channels.
%
% Args:
%   inputFile (char, required, positional): a shape-(1, n) character
%     array containing the input NWB file name with a full path. It must
%     end with '.nwb'.
%   timeseriesGroup (char | cell, optional, keyword): a shape-(1, m)
%     character or cell array containing the NWB timeseries group or
%     multiple groups to be converted into the binary format. If not
%     specified, loads microelectrode data.
%   filter (logical, optional, keyword): a shape-(1, 1) logical scalar
%     indicating whether to bandpass filter the data (default=false).
%   freqRange (numeric, optional, keyword): a shape-(1, 2) numeric array
%     representing the low and high cutoff filter frequencies in Hz. It is
%     only used if the filter is set to true (default=[300 3000]).
%   removeArtifacts (logical, optional, keyword): a shape-(1, 1) logical
%     scalar indicating whether to remove high amplitude noise periods in
%     recordings (default=false).
%   artifactCutoff (numeric, optional, keyword): a shape-(1, 2) numeric
%     array of upper and lower voltage cutoff values in microvolts used
%     when removeArtifacts is set to true (default=[1500 -1500]).
%   gain (numeric, optional, keyword): a shape-(1, 1) numeric scalar used
%     to multiple data to convert it into appropriate
%     precision format (default=1).
%   car (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling the application of the common average reference (CAR).
%     Referencing is applied on a per lead basis described in
%     channelGroups (default=false).
%   channelGroups (numeric, optional, keyword): a shape-(k, l) numeric
%     array of channel groups. Rows correspond to tetrodes and columns
%     correspond to individual channels. By default assumes consequtive
%     tetrode channel pairings.
%   segmentSize (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     defining data segment size to be extracted in samples (default=3e8).
%   segmentStartSample (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar defining data segment start sample (default=1).
%   verbose (logical, optional, keyword): a shape-(1, 1) logical scalar
%     controlling the verbosity level of the function (default=false).
%
% Returns:
%   timeseriesData (numeric): a shape-(g, h) numeric row matrix containing
%     the data output with rows corresponding to individual recording
%     channels.
%   timestamps (numeric): a shape-(1, h) numeric array with timestamps
%     corresponding to the columns of timeseriesData.
%
% Dependencies:
%   matnwb (https://neurodatawithoutborders.github.io/matnwb/)
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  inputFile (1,:) {mustBeA(inputFile,'char'),mustBeVector,endsWith(inputFile,'.nwb')}
  options.timeseriesGroup (1,:) {mustBeA(options.timeseriesGroup,'char')} = 'TimeSeries_32000_Hz'
  options.filter (1,1) {islogical} = false
  options.freqRange (1,2) {mustBePositive} = [300 3000]
  options.removeArtifacts (1,1) {islogical} = false
  options.artifactCutoff (1,:) {mustBeNumeric} = [1500 -1500]
  options.gain (1,1) {mustBePositive} = 1
  options.car (1,1) {islogical} = false
  options.channelGroups (:,:) {mustBePositive} = [1:4; 5:8; 9:12; 13:16; 17:20; 21:24; 25:28; 29:32; 33:36; 37:40; 41:44; 45:48; 49:52; 53:56; 57:60; 61:64];
  options.segmentSize (1,1) {mustBeNumeric} = 3e8
  options.segmentStartSample (1,1) {mustBeNumeric} = 1
  options.verbose (1,1) {islogical} = false
end

% Parameters
options.absAmpSmoothCutoffFactor = 1.25; % baseline factor
options.artifactExpandSamples = 20; % samples
options.artifactForwardExpandSamples = 20; % samples
options.largeArtifactExpandTime = 0.05; % seconds
options.baselinePrctile = 25; % percentile

% Load the nwb file
if options.verbose
  disp('Loading NWB data')
end
nwbData = nwbRead(inputFile);
dataContainer = nwbData.acquisition.get(options.timeseriesGroup);

% Get the data
segmentInds = (1:options.segmentSize) + options.segmentStartSample - 1;
segmentInds(segmentInds > size(dataContainer.data,2)) = [];
timeseriesData = dataContainer.data(:,segmentInds);
timeseriesData = double(timeseriesData).*options.gain;
samplingRate = dataContainer.starting_time_rate;
if isempty(samplingRate)
  timestamps = dataContainer.timestamps(segmentInds)';
  timestamps = timestamps - dataContainer.timestamps(1);
  samplingRate = 1/median(diff(timestamps));
else
  timestamps = segmentInds./samplingRate;
end

% Filter data
if options.filter
  if options.verbose
    disp('Filtering data')
  end
  timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
    sampleRate=samplingRate, frequencyRange=options.freqRange);
end

% Subtract CAR
nChans = size(timeseriesData,1);
options.channelGroups(nChans/4+1:end,:) = [];
if options.car
  if options.verbose
    disp('Applying common average reference')
  end
  if size(options.channelGroups,1) == 1
    leadGroups = options.channelGroups(1,:);
  elseif size(options.channelGroups,1) == 2
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:)];
  elseif size(options.channelGroups,1) == 4
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:)];
  elseif size(options.channelGroups,1) == 6
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:)];
  elseif size(options.channelGroups,1) == 8
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:); ...
      options.channelGroups(7,:) options.channelGroups(8,:)];
  elseif size(options.channelGroups,1) == 10
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:); ...
      options.channelGroups(7,:) options.channelGroups(8,:); ...
      options.channelGroups(9,:) options.channelGroups(10,:)];
  elseif size(options.channelGroups,1) == 12
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:); ...
      options.channelGroups(7,:) options.channelGroups(8,:); ...
      options.channelGroups(9,:) options.channelGroups(10,:); ...
      options.channelGroups(11,:) options.channelGroups(12,:)];
  elseif size(options.channelGroups,1) == 14
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:); ...
      options.channelGroups(7,:) options.channelGroups(8,:); ...
      options.channelGroups(9,:) options.channelGroups(10,:); ...
      options.channelGroups(11,:) options.channelGroups(12,:); ...
      options.channelGroups(13,:) options.channelGroups(14,:)];
  elseif size(options.channelGroups,1) == 16
    leadGroups = [options.channelGroups(1,:) options.channelGroups(2,:); ...
      options.channelGroups(3,:) options.channelGroups(4,:); ...
      options.channelGroups(5,:) options.channelGroups(6,:); ...
      options.channelGroups(7,:) options.channelGroups(8,:); ...
      options.channelGroups(9,:) options.channelGroups(10,:); ...
      options.channelGroups(11,:) options.channelGroups(12,:); ...
      options.channelGroups(13,:) options.channelGroups(14,:); ...
      options.channelGroups(15,:) options.channelGroups(16,:)];
  end
  for iLead = 1:size(leadGroups,1)
    [~, chanInds] = unique(leadGroups(iLead,:));
    medianChannel = median(timeseriesData(leadGroups(iLead,chanInds),:), 2, 'omitnan');
    timeseriesData(leadGroups(iLead,chanInds),:) = bsxfun( ...
      @minus, timeseriesData(leadGroups(iLead,chanInds),:), medianChannel);
    medianTrace = median(timeseriesData(leadGroups(iLead,chanInds),:), 1, 'omitnan');
    timeseriesData(leadGroups(iLead,chanInds),:) = bsxfun( ...
      @minus, timeseriesData(leadGroups(iLead,chanInds),:), medianTrace);
  end
end

% Zero out stimulation data (or artifacts)
if options.verbose && options.removeArtifacts
  disp('Removing large amplitude noise')
end
baseline = 0;
for iChan = 1:nChans
  nanMask = isnan(timeseriesData(iChan,:));
  timeseriesData(iChan,nanMask) = 0;
  if options.removeArtifacts
    if options.verbose
      disp(['   Processing channel ' num2str(iChan) '/' num2str(nChans)])
    end
    if ismember(iChan, options.channelGroups(:,1))
      stimMask = false(1,size(timeseriesData,2));
    end
    absTimeseriesSmooth = smooth(abs(timeseriesData(iChan,:)), round(samplingRate/1));
    baselineSmooth = prctile(absTimeseriesSmooth, options.baselinePrctile);
    stimMask(absTimeseriesSmooth >= baselineSmooth*options.absAmpSmoothCutoffFactor) = true;
    if ~isempty(options.artifactCutoff)
      stimMask(timeseriesData(iChan,:) >= options.artifactCutoff(1)) = true;
      stimMask(timeseriesData(iChan,:) <= options.artifactCutoff(2)) = true;
    end
    chanGroupMask = ismember(options.channelGroups(:,end), iChan);
    if any(chanGroupMask)
      for iShift = 1:options.artifactExpandSamples
        if iShift > options.artifactForwardExpandSamples
          stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
        else
          stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
            [false(1,iShift) stimMask(1:end-iShift)]);
        end
      end
      diffStimMask = [0 diff(stimMask)];
      onsets = find(diffStimMask > 0);
      offsets = find(diffStimMask < 0);
      if isempty(onsets) && ~isempty(offsets)
        onsets = 1;
      elseif isempty(onsets) && isempty(offsets) && diffStimMask(1)
        onsets = 1;
        offsets = numel(diffStimMask);
      elseif ~isempty(onsets) && isempty(offsets)
        offsets = numel(diffStimMask);
      elseif isempty(onsets) && isempty(offsets)
        continue
      end
      if onsets(1) > offsets(1)
        onsets = [1 onsets]; %#ok<*AGROW>
      end
      if offsets(end) < onsets(end)
        offsets = [offsets numel(stimMask)];
      end
      assert(numel(onsets) == numel(offsets));
      durations = (offsets-onsets)./samplingRate;
      for iDuration = 1:numel(durations)
        if durations(iDuration) >= 1
          stimMask(max([1 onsets(iDuration) - round(options.largeArtifactExpandTime*samplingRate)]): ...
            min([offsets(iDuration) + round(options.largeArtifactExpandTime*samplingRate) numel(stimMask)])) = true;
        end
      end
      if any(stimMask)
        %figure;
        %plot(timeseriesData./max(timeseriesData(iChan,:)));
        %hold on; plot(stimMask); hold off
        timeseriesData(options.channelGroups(chanGroupMask,:),stimMask) = baseline;
      end
    end
  end
end