function nwb2binary(inputFile, outputFile, options)
% nwb2binary(inputFile, outputFile, <options>)
%
% Function loads a file in NWB format and saves it as a binary file (.dat).
% Individual timeseries groups are saved in separate files.
%
% Args:
%   inputFile (char, required, positional): a shape-(1, n) character
%     array containing the input NWB file name with a full path. It must
%     end with '.nwb'.
%   outputFile (char, required, positional): a shape-(1, n) character
%     array containing the output file name with a full path. It must
%     end with '.dat'. Individual timeseries groups will be saved with
%     filenames appended by the timeseries group name.
%   timeseriesGroup (char | cell, optional, keyword): a shape-(1, k)
%     character or cell array containing the NWB timeseries group or
%     multiple groups to be converted into the binary format. If not
%     specified, all timeseries groups will be saved in separate binary
%     files.
%   writeBufferSize (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing the size of the binary write buffer in bytes.
%     To improve the processing speed, it should be a multiple of 4096 as
%     this is the size of most hard drive buffers. Default is 500*4096
%     bytes.
%   precision (char, optional, keyword): a shape-(1, m) character array
%     describing the format that the binary should be saved in. By default,
%     the binary is saved in the same format as it was stored in the NWB
%     file. Available formats include 'int8', 'uint8', 'int16', 'uint16',
%     'int32', 'uint32', 'single', 'int64', 'uint64', and 'double'.
%   filter (char, optional, keyword): a shape-(1, l) character array
%     describing the data filter type, if any. Available options include:
%     'bandpass_prior', 'lowpass_prior', 'highpass_prior', 'bandpass_post',
%     'lowpass_post', 'highpass_post', 'none' (default).
%   lowCutoffFreq (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing the low cutoff filter frequency in Hz. It is only
%     used if the filter is specified as 'bandpass' or 'lowpass'
%     (default=300).
%   highCutoffFreq (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar representing the high cutoff filter frequency in Hz. It is
%     only used if the filter is specified as 'bandpass' or 'highpass'
%     (default=6000).
%   samplingRate (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar corresponding to data sampling frequency in Hz
%     (default=32000).
%   removeArtifacts (char, optional, keyword): a shape-(1, j) character
%     array describing the method for artifact removal. Available options
%     include:
%       'none' - no artifact removal (default);
%       'stim' - artifact removal based on stimulation channels (not
%                recommended). This option is not available when procByChan
%                is set to true;
%       'recNormThr' - removal of absolute voltage periods exceeding a
%                      threshold set relative to the maximum value.
%       'recAbsThr' - removal of voltage periods exceeding an explicitly
%                     set double threshold in positive and negative
%                     directions.
%       'amp' - removal of voltage periods exceeding with sustained high
%               amplitude fluctuations (recommended).
%   artifactCutoff (numeric, optional, keyword): a shape-(1, 2) numeric
%     array of voltage cutoff values used when 'recNormThr' or 'recAbsThr'
%     options are chosen for removeArtifacts. The second value is only used
%     in the latter case. Default is [0.01 0.01].
%   zeroPeriods (numeric, optional, keyword): a shape-(l, 3) numeric array
%     containing time periods to be zeroed out in the data. Rows of this
%     matrix corresponds to individual time periods, the first column
%     corresponds to the channel number (use 0 to indicate all channels),
%     the second column corresponds to start times, whereas the third one
%     corresponds to end times. If left empty, no data will be zeroed out
%     (default).
%   visualiseData (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for visualising data that is being saved to a binary file
%     (default=false).
%   dataConversionFactor (numeric, optional, keyword): a shape-(1, 1)
%     numeric scalar used to multiple data to convert it into appropriate
%     precision format (default=32000).
%   conversionMode (char, optional, keyword): a shape-(1, i) character
%     array describing the method of binary conversion. Available options
%     include:
%       'full' - loading full data matrix from NWB and converting it into
%                binary (default). This is the best approach in converting
%                small NWB files;
%       'channel' - loading a single data matrix row (channel) from the NWB
%                   file, processing it, and storing it in a data container
%                   in an integer format before all channels are loaded and
%                   processed, and finally saving them in the binary format.
%                   This approach works with intermediate size files;
%       'segment' - iteratively loading all channels of the data matrix but
%                   of a defined size segment and converting it into the
%                   binary format. This is the best method for dealing with
%                   arbitrarily large NWB files.
%
% Returns:
%   None
%
% Dependencies:
%   matnwb (https://neurodatawithoutborders.github.io/matnwb/)
%   dervinism/dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab).
%   buzsakilab/buzcode (https://github.com/buzsakilab/buzcode).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  inputFile (1,:) {mustBeA(inputFile,'char'),mustBeVector,endsWith(inputFile,'.nwb')}
  outputFile (1,:) {mustBeA(outputFile,'char'),mustBeVector,endsWith(outputFile,'.dat')}
  options.timeseriesGroup (:,:) {mustBeCharOrListedType(options.timeseriesGroup, 'cell')} = ''
  options.writeBufferSize (1,1) {mustBePositive} = 500*4096
  options.precision (1,:) {mustBeA(options.precision,'char'),mustBeVector} = ''
  options.filter (1,:) {mustBeMember(options.filter,{'none','bandpass_prior','lowpass_prior','highpass_prior','bandpass_post','lowpass_post','highpass_post'})} = 'none'
  options.lowCutoffFreq (1,1) {mustBePositive} = 300
  options.highCutoffFreq (1,1) {mustBePositive} = 6000
  options.removeArtifacts (1,:) {mustBeMember(options.removeArtifacts,{'none','stim','recNormThr','recAbsThr','amp'})} = 'none'
  options.artifactCutoff (1,:) {mustBeNumeric} = [0.01 0.01]
  options.zeroPeriods (:,:) {mustBeNumeric} = []
  options.visualiseData (1,1) {islogical} = false
  options.dataConversionFactor (1,1) {mustBePositive} = 1
  options.conversionMode (1,:) {mustBeMember(options.conversionMode,{'full','channel','segment'})} = 'full'
end

% Parse input
if ~isempty(options.timeseriesGroup) && ischar(options.timeseriesGroup)
  options.timeseriesGroup = {options.timeseriesGroup};
end

% Parameters
segmentSize = 1e8; % samples. 1875 seconds for 32kHz sampling frequency

% Load NWB file
nwbData = nwbRead(inputFile);

% Determine which timeseries groups to convert
timeseriesGroups = nwbData.acquisition.keys;
if isempty(options.timeseriesGroup)
  options.timeseriesGroup = timeseriesGroups;
end

if strcmpi(options.conversionMode, 'segment')
  % Convert each timeseries group
  nGroups = numel(options.timeseriesGroup);
  for iGroup = 1:nGroups
    dataContainer = nwbData.acquisition.get(options.timeseriesGroup{iGroup});
    % Convert each segment of a timeseries group
    nSegments = ceil(size(dataContainer.data,2)/segmentSize);
    for iSegment = 1:nSegments
      segmentInds = (1:segmentSize) + (iSegment-1)*segmentSize;
      segmentInds(segmentInds > size(dataContainer.data,2)) = [];
      timeseriesData = dataContainer.data(:,segmentInds);
      samplingRate = dataContainer.starting_time_rate;
      if isempty(samplingRate)
        timestamps = dataContainer.timestamps(segmentInds)';
        timestamps = timestamps - dataContainer.timestamps(1);
        samplingRate = 1/median(diff(timestamps));
      else
        timestamps = segmentInds./samplingRate;
      end

      % Convert to double
      if isfloat(timeseriesData)
        warning(['Timeseries data type is float. This is likely an error that ' ...
          'will result in incorrect binary data conversion. Or you may need ' ...
          'to provide appropriate dataConversionFactor.']);
      end
      dataType = class(timeseriesData);
      timeseriesData = double(timeseriesData);
      timeseriesData = timeseriesData.*options.dataConversionFactor;

      % Filter data
      SD_prior = mean(std(timeseriesData,[],2));
      if strcmpi(options.filter, 'bandpass_prior')
        timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
      elseif strcmpi(options.filter, 'lowpass_prior')
        timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.lowCutoffFreq);
      elseif strcmpi(options.filter, 'highpass_prior')
        timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.highCutoffFreq);
      end
      SD_post = mean(std(timeseriesData,[],2));
      timeseriesData = timeseriesData.*(SD_prior/SD_post);

      % Visualise data
      nChans = size(timeseriesData,1);
      if options.visualiseData
        for iChan = 1:nChans
          figure; plot(timestamps, timeseriesData(iChan,:));
          %figure; MTSpectrogram(timeseriesData(iChan,:)', ...
          %  'frequency',samplingRate, 'show','on');
        end
      end

      % Zero out stimulation data (or artifacts)
      stimGroupMask = contains(timeseriesGroups, 'Stim'); %#ok<*NASGU>
      if ~strcmpi(options.removeArtifacts, 'none')
        if strcmpi(options.removeArtifacts, 'stim')
          if sum(stimGroupMask) %#ok<*UNRCH>
            stimGroupName = timeseriesGroups{stimGroupMask};
            dataContainer = nwbData.acquisition.get(stimGroupName);
            stimTimeseriesData = dataContainer.data(:,:);
            thr = min(max(stimTimeseriesData,[],2))/4;
            stimMask = false(size(stimTimeseriesData));
            stimMask(stimTimeseriesData >= thr) = true;
            stimMask = sum(stimMask);
            stimMask = logical([stimMask(2:end) false] + [false stimMask(1:end-1)]);
            stimSamplingRate = dataContainer.starting_time_rate;
            if isempty(stimSamplingRate)
              stimTimestamps = dataContainer.timestamps(:)';
              stimSamplingRate = 1/median(diff(stimTimestamps));
            end
            if any(stimMask)
              for iChan = 1:nChans
                timeseriesData(iChan,stimMask) = median(timeseriesData(iChan,:));
              end
            end
          end
        elseif strcmpi(options.removeArtifacts, 'recNormThr')
          for iChan = 1:nChans
            baseline = median(timeseriesData(iChan,:));
            if strcmpi(options.filter, 'none') || ~strcmpi(options.filter, 'bandpass_prior')
              absTimeseries = abs(bandpassFilterTimeSeries(double(timeseriesData(iChan,:)), ...
                sampleRate=samplingRate, frequencyRange=[1 9500]));
            else
              absTimeseries = abs(timeseriesData(iChan,:));
            end
            stimMask = absTimeseries./max(absTimeseries);
            thr = options.artifactCutoff;
            stimMask(stimMask >= thr) = 1;
            stimMask(stimMask < thr) = 0;
            stimMask = logical(stimMask);
            for iShift = 1:20
              if iShift > 20
                stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
              else
                stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
                  [false(1,iShift) stimMask(1:end-iShift)]);
              end
            end
            if any(stimMask)
              %figure;
              %plot(timeseriesData(iChan,:)./max(timeseriesData(iChan,:)));
              %hold on; plot(stimMask); hold off
              timeseriesData(iChan,stimMask) = baseline;
            end
          end
        elseif strcmpi(options.removeArtifacts, 'recAbsThr')
          baseline = 0;
          for iChan = 1:nChans
            stimMask = false(1,size(timeseriesData,2));
            stimMask(timeseriesData(iChan,:) >= options.artifactCutoff(1)) = true;
            stimMask(timeseriesData(iChan,:) <= options.artifactCutoff(2)) = true;
            for iShift = 1:20
              if iShift > 20
                stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
              else
                stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
                  [false(1,iShift) stimMask(1:end-iShift)]);
              end
            end
            if any(stimMask)
              %figure;
              %plot(timeseriesData./max(timeseriesData(iChan,:)));
              %hold on; plot(stimMask); hold off
              timeseriesData(iChan,stimMask) = baseline;
            end
          end
        elseif strcmpi(options.removeArtifacts, 'amp')
          baseline = 0;
          for iChan = 1:nChans
            stimMask = false(1,size(timeseriesData,2));
            stimMask(timeseriesData(iChan,:) >= options.artifactCutoff(1)) = true;
            stimMask(timeseriesData(iChan,:) <= options.artifactCutoff(2)) = true;
            absTimeseriesSmooth = smooth(abs(timeseriesData(iChan,:)), round(samplingRate/10));
            baselineSmooth = median(absTimeseriesSmooth);
            stimMask(absTimeseriesSmooth >= baselineSmooth*2) = true;
            for iShift = 1:20
              if iShift > 20
                stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
              else
                stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
                  [false(1,iShift) stimMask(1:end-iShift)]);
              end
            end
            diffStimMask = [0 diff(stimMask)];
            onsets = find(diffStimMask > 0);
            offsets = find(diffStimMask < 0);
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
                stimMask(max([1 onsets(iDuration) - round(0.1*samplingRate)]): ...
                  min([offsets(iDuration) + round(0.1*samplingRate) numel(stimMask)])) = true;
              end
            end
            if any(stimMask)
              %figure;
              %plot(timeseriesData./max(timeseriesData(iChan,:)));
              %hold on; plot(stimMask); hold off
              timeseriesData(iChan,stimMask) = baseline;
            end
          end
        end
      end

      % Zero out time periods of choice
      if ~isempty(options.zeroPeriods)
        nPeriods = size(options.zeroPeriods,1);
        for iPeriod = 1:nPeriods
          if ~options.zeroPeriods(iPeriod,2)
            inds(1) = 1;
          else
            inds(1) = find(timestamps - options.zeroPeriods(iPeriod,2) > 0, 1);
          end
          if isinf(options.zeroPeriods(iPeriod,3))
            inds(2) = inf;
          else
            inds(2) = find(timestamps - options.zeroPeriods(iPeriod,3) > 0, 1);
          end
          if options.zeroPeriods(iPeriod,1)
            if isinf(inds(2))
              timeseriesData(options.zeroPeriods(iPeriod,1),inds(1):end) = ...
                median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
            else
              timeseriesData(options.zeroPeriods(iPeriod,1),inds(1):inds(2)) = ...
                median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
            end
          else
            for iChan = 1:nChans
              if isinf(inds(2))
                timeseriesData(iChan,inds(1):end) = median(timeseriesData(iChan,:));
              else
                timeseriesData(iChan,inds(1):inds(2)) = median(timeseriesData(iChan,:));
              end
            end
          end
        end
      end

      % Filter data
      SD_prior = mean(std(timeseriesData,[],2));
      if strcmpi(options.filter, 'bandpass_post')
        timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
      elseif strcmpi(options.filter, 'lowpass_post')
        timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.lowCutoffFreq);
      elseif strcmpi(options.filter, 'highpass_post')
        timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.highCutoffFreq);
      end
      SD_post = mean(std(timeseriesData,[],2));
      timeseriesData = timeseriesData.*(SD_prior/SD_post);

      % Visualise data
      if options.visualiseData
        for iChan = 1:nChans
          figure; plot(timestamps, timeseriesData(iChan,:));
          %figure; MTSpectrogram(timeseriesData(iChan,:)', ...
          %  'frequency',samplingRate, 'show','on');
        end
      end

      % Save the binary file
      [~,~,ext] = fileparts(outputFile);
      outputFile_group = [outputFile(1:end-4) '_' options.timeseriesGroup{iGroup} ext];
      if strcmpi(dataType, options.precision) || ~isempty(options.precision)
        timeseriesData = cast(timeseriesData, options.precision);
      else
        timeseriesData = cast(timeseriesData, dataType);
      end
      if iSegment == 1
        fid = fopen(outputFile_group, 'w');
      else
        fid = fopen(outputFile_group, 'ab');
      end
      fwrite(fid, timeseriesData, options.precision);
      if iSegment == nSegments
        fclose(fid);
      end
    end
  end





elseif strcmpi(options.conversionMode, 'channel')

  % Convert each timeseries group
  nGroups = numel(options.timeseriesGroup);
  for iGroup = 1:nGroups
    dataContainer = nwbData.acquisition.get(options.timeseriesGroup{iGroup});
    nChans = size(dataContainer.data,1);
    samplingRate = dataContainer.starting_time_rate;
    if isempty(samplingRate)
      timestamps = dataContainer.timestamps(:)';
      timestamps = timestamps - timestamps(1);
      samplingRate = 1/median(diff(timestamps));
    else
      timestamps = (1:size(dataContainer.data,2))./samplingRate;
    end

    % Convert channel by channel
    timeseriesData2write = zeros(size(dataContainer.data), 'int32');
    for iChan = 1:nChans
      timeseriesData = dataContainer.data(iChan,:);

      % Convert to double
      if isfloat(timeseriesData)
        warning(['Timeseries data type is float. This is likely an error that ' ...
          'will result in incorrect binary data conversion. Or you may need ' ...
          'to provide appropriate dataConversionFactor.']);
      end
      dataType = class(timeseriesData);
      timeseriesData = double(timeseriesData);
      timeseriesData = timeseriesData.*options.dataConversionFactor;

      % Filter data
      SD_prior = std(timeseriesData);
      if strcmpi(options.filter, 'bandpass_prior')
        timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
      elseif strcmpi(options.filter, 'lowpass_prior')
        timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.lowCutoffFreq);
      elseif strcmpi(options.filter, 'highpass_prior')
        timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
          sampleRate=samplingRate, ...
          frequencyCutoff=options.highCutoffFreq);
      end
      SD_post = std(timeseriesData);
      timeseriesData = timeseriesData.*(SD_prior/SD_post);

      % Visualise data
      if options.visualiseData
        %figure; plot(timestamps, timeseriesData);
        %figure; MTSpectrogram(double(timeseriesData'), ...
        %  'frequency',samplingRate, 'show','on');
      end

      % Zero out stimulation data (or artifacts)
      stimGroupMask = contains(timeseriesGroups, 'Stim'); %#ok<*NASGU>
      if strcmpi(options.removeArtifacts, 'recNormThr')
        baseline = median(timeseriesData);
        if strcmpi(options.filter, 'none') || ~strcmpi(options.filter, 'bandpass_prior')
          absTimeseries = abs(bandpassFilterTimeSeries(timeseriesData, ...
            sampleRate=samplingRate, frequencyRange=[1 9500]));
        else
          absTimeseries = abs(timeseriesData);
        end
        stimMask = absTimeseries./max(absTimeseries);
        thr = options.artifactCutoff;
        stimMask(stimMask >= thr) = 1;
        stimMask(stimMask < thr) = 0;
        stimMask = logical(stimMask);
        for iShift = 1:20
          if iShift > 20
            stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
          else
            stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
              [false(1,iShift) stimMask(1:end-iShift)]);
          end
        end
        if any(stimMask)
          %figure;
          %plot(timeseriesData./max(timeseriesData));
          %hold on; plot(stimMask); hold off
          timeseriesData(stimMask) = baseline;
        end
      elseif strcmpi(options.removeArtifacts, 'recAbsThr')
        baseline = 0;
        stimMask = false(size(timeseriesData));
        stimMask(timeseriesData >= options.artifactCutoff(1)) = true;
        stimMask(timeseriesData <= options.artifactCutoff(2)) = true;
        for iShift = 1:20
          if iShift > 20
            stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
          else
            stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
              [false(1,iShift) stimMask(1:end-iShift)]);
          end
        end
        if any(stimMask)
          %figure;
          %plot(timeseriesData./max(timeseriesData));
          %hold on; plot(stimMask); hold off
          timeseriesData(stimMask) = baseline;
        end
      elseif strcmpi(options.removeArtifacts, 'amp')
        baseline = 0;
        stimMask = false(size(timeseriesData));
        stimMask(timeseriesData >= options.artifactCutoff(1)) = true;
        stimMask(timeseriesData <= options.artifactCutoff(2)) = true;
        absTimeseriesSmooth = smooth(abs(timeseriesData), round(samplingRate/10));
        baselineSmooth = median(absTimeseriesSmooth);
        stimMask(absTimeseriesSmooth >= baselineSmooth*2) = true;
        for iShift = 1:20
          if iShift > 20
            stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
          else
            stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
              [false(1,iShift) stimMask(1:end-iShift)]);
          end
        end
        diffStimMask = [0 diff(stimMask)];
        onsets = find(diffStimMask > 0);
        offsets = find(diffStimMask < 0);
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
            stimMask(max([1 onsets(iDuration) - round(0.1*samplingRate)]): ...
              min([offsets(iDuration) + round(0.1*samplingRate) numel(stimMask)])) = true;
          end
        end
        if any(stimMask)
          %figure;
          %plot(timeseriesData./max(timeseriesData));
          %hold on; plot(stimMask); hold off
          timeseriesData(stimMask) = baseline;
        end
      end
      clear stimMask diffStimMask absTimeseries absTimeseriesSmooth onsets
      clear offsets durations SD_post SD_prior baseline baselineSmooth
      clear iDuration iShift

      % Zero out time periods of choice
      if ~isempty(options.zeroPeriods)
        nPeriods = size(options.zeroPeriods,1);
        for iPeriod = 1:nPeriods
          if options.zeroPeriods(iPeriod,1) == iChan || ~options.zeroPeriods(iPeriod,1)
            if ~options.zeroPeriods(iPeriod,2)
              inds(1) = 1;
            else
              inds(1) = find(timestamps - options.zeroPeriods(iPeriod,2) > 0, 1);
            end
            if isinf(options.zeroPeriods(iPeriod,3))
              inds(2) = inf;
            else
              inds(2) = find(timestamps - options.zeroPeriods(iPeriod,3) > 0, 1);
            end
            if isinf(inds(2))
              timeseriesData(inds(1):end) = ...
                median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
            else
              timeseriesData(inds(1):inds(2)) = ...
                median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
            end
          end
        end
      end

      % Filter data
      if strcmpi(options.filter, 'bandpass_post') || ...
          strcmpi(options.filter, 'lowpass_post') || ...
          strcmpi(options.filter, 'highpass_post')
        SD_prior = std(timeseriesData);
        if strcmpi(options.filter, 'bandpass_post')
          timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
            sampleRate=samplingRate, ...
            frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
        elseif strcmpi(options.filter, 'lowpass_post')
          timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
            sampleRate=samplingRate, ...
            frequencyCutoff=options.lowCutoffFreq);
        elseif strcmpi(options.filter, 'highpass_post')
          timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
            sampleRate=samplingRate, ...
            frequencyCutoff=options.highCutoffFreq);
        end
        SD_post = mean(std(timeseriesData,[],2));
        timeseriesData = timeseriesData.*(SD_prior/SD_post);
      end

      % Visualise data
      if options.visualiseData
        figure; plot(timestamps, timeseriesData);
        %figure; MTSpectrogram(timeseriesData', ...
        %  'frequency',samplingRate, 'show','on');
      end

      % Save the channel MAT file
      outputFile_ch = [outputFile(1:end-4) '_' options.timeseriesGroup{iGroup} '_' num2str(iChan) '.mat'];
      timeseriesData = cast(timeseriesData, dataType);
      if ~isempty(options.precision)
        timeseriesData2write(iChan,:) = cast(timeseriesData, options.precision);
      else
        timeseriesData2write(iChan,:) = timeseriesData;
      end
      %save(outputFile_ch, 'timeseriesData', '-v7.3');
    end

    % Save the binary file
    [~,~,ext] = fileparts(outputFile);
    outputFile_group = [outputFile(1:end-4) '_' options.timeseriesGroup{iGroup} ext];
    writeBinary(timeseriesData, outputFile_group, writeBufferSize=options.writeBufferSize);
  end





elseif strcmpi(options.conversionMode, 'full')

  % Convert each timeseries group
  nGroups = numel(options.timeseriesGroup);
  for iGroup = 1:nGroups
    dataContainer = nwbData.acquisition.get(options.timeseriesGroup{iGroup});
    timeseriesData = dataContainer.data(:,:);
    samplingRate = dataContainer.starting_time_rate;
    if isempty(samplingRate)
      timestamps = dataContainer.timestamps(:)';
      timestamps = timestamps - timestamps(1);
      samplingRate = 1/median(diff(timestamps));
    else
      timestamps = (1:size(timeseriesData,2))./samplingRate;
    end

    % Convert to double
    if isfloat(timeseriesData)
      warning(['Timeseries data type is float. This is likely an error that ' ...
        'will result in incorrect binary data conversion. Or you may need ' ...
        'to provide appropriate dataConversionFactor.']);
    end
    dataType = class(timeseriesData);
    timeseriesData = double(timeseriesData);
    timeseriesData = timeseriesData.*options.dataConversionFactor;

    % Filter data
    SD_prior = mean(std(timeseriesData,[],2));
    if strcmpi(options.filter, 'bandpass_prior')
      timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
    elseif strcmpi(options.filter, 'lowpass_prior')
      timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyCutoff=options.lowCutoffFreq);
    elseif strcmpi(options.filter, 'highpass_prior')
      timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyCutoff=options.highCutoffFreq);
    end
    SD_post = mean(std(timeseriesData,[],2));
    timeseriesData = timeseriesData.*(SD_prior/SD_post);

    % Visualise data
    nChans = size(timeseriesData,1);
    if options.visualiseData
      for iChan = 1:nChans
        figure; plot(timestamps, timeseriesData(iChan,:));
        %figure; MTSpectrogram(timeseriesData(iChan,:)', ...
        %  'frequency',samplingRate, 'show','on');
      end
    end

    % Zero out stimulation data (or artifacts)
    stimGroupMask = contains(timeseriesGroups, 'Stim'); %#ok<*NASGU>
    if ~strcmpi(options.removeArtifacts, 'none')
      if strcmpi(options.removeArtifacts, 'stim')
        if sum(stimGroupMask) %#ok<*UNRCH>
          stimGroupName = timeseriesGroups{stimGroupMask};
          dataContainer = nwbData.acquisition.get(stimGroupName);
          stimTimeseriesData = dataContainer.data(:,:);
          thr = min(max(stimTimeseriesData,[],2))/4;
          stimMask = false(size(stimTimeseriesData));
          stimMask(stimTimeseriesData >= thr) = true;
          stimMask = sum(stimMask);
          stimMask = logical([stimMask(2:end) false] + [false stimMask(1:end-1)]);
          stimSamplingRate = dataContainer.starting_time_rate;
          if isempty(stimSamplingRate)
            stimTimestamps = dataContainer.timestamps(:)';
            stimSamplingRate = 1/median(diff(stimTimestamps));
          end
          if any(stimMask)
            for iChan = 1:nChans
              timeseriesData(iChan,stimMask) = median(timeseriesData(iChan,:));
            end
          end
        end
      elseif strcmpi(options.removeArtifacts, 'recNormThr')
        for iChan = 1:nChans
          baseline = median(timeseriesData(iChan,:));
          if strcmpi(options.filter, 'none') || ~strcmpi(options.filter, 'bandpass_prior')
            absTimeseries = abs(bandpassFilterTimeSeries(double(timeseriesData(iChan,:)), ...
              sampleRate=samplingRate, frequencyRange=[1 9500]));
          else
            absTimeseries = abs(timeseriesData(iChan,:));
          end
          stimMask = absTimeseries./max(absTimeseries);
          thr = options.artifactCutoff;
          stimMask(stimMask >= thr) = 1;
          stimMask(stimMask < thr) = 0;
          stimMask = logical(stimMask);
          for iShift = 1:20
            if iShift > 20
              stimMask = logical(stimMask + [false(1,iShift) stimMask(1:end-iShift)]);
            else
              stimMask = logical([stimMask(1+iShift:end) false(1,iShift)] + ...
                [false(1,iShift) stimMask(1:end-iShift)]);
            end
          end
          if any(stimMask)
            %figure;
            %plot(timeseriesData(iChan,:)./max(timeseriesData(iChan,:)));
            %hold on; plot(stimMask); hold off
            timeseriesData(iChan,stimMask) = baseline;
          end
        end
      end
    end

    % Zero out time periods of choice
    if ~isempty(options.zeroPeriods)
      nPeriods = size(options.zeroPeriods,1);
      for iPeriod = 1:nPeriods
        if ~options.zeroPeriods(iPeriod,2)
          inds(1) = 1;
        else
          inds(1) = find(timestamps - options.zeroPeriods(iPeriod,2) > 0, 1);
        end
        if isinf(options.zeroPeriods(iPeriod,3))
          inds(2) = inf;
        else
          inds(2) = find(timestamps - options.zeroPeriods(iPeriod,3) > 0, 1);
        end
        if options.zeroPeriods(iPeriod,1)
          if isinf(inds(2))
            timeseriesData(options.zeroPeriods(iPeriod,1),inds(1):end) = ...
              median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
          else
            timeseriesData(options.zeroPeriods(iPeriod,1),inds(1):inds(2)) = ...
              median(timeseriesData(options.zeroPeriods(iPeriod,1),:));
          end
        else
          for iChan = 1:nChans
            if isinf(inds(2))
              timeseriesData(iChan,inds(1):end) = median(timeseriesData(iChan,:));
            else
              timeseriesData(iChan,inds(1):inds(2)) = median(timeseriesData(iChan,:));
            end
          end
        end
      end
    end

    % Filter data
    SD_prior = mean(std(timeseriesData,[],2));
    if strcmpi(options.filter, 'bandpass_post')
      timeseriesData = bandpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyRange=[options.lowCutoffFreq options.highCutoffFreq]);
    elseif strcmpi(options.filter, 'lowpass_post')
      timeseriesData = lowpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyCutoff=options.lowCutoffFreq);
    elseif strcmpi(options.filter, 'highpass_post')
      timeseriesData = highpassFilterTimeSeries(timeseriesData, ...
        sampleRate=samplingRate, ...
        frequencyCutoff=options.highCutoffFreq);
    end
    SD_post = mean(std(timeseriesData,[],2));
    timeseriesData = timeseriesData.*(SD_prior/SD_post);

    % Visualise data
    if options.visualiseData
      for iChan = 1:nChans
        figure; plot(timestamps, timeseriesData(iChan,:));
        %figure; MTSpectrogram(timeseriesData(iChan,:)', ...
        %  'frequency',samplingRate, 'show','on');
      end
    end

    % Save the binary file
    [~,~,ext] = fileparts(outputFile);
    outputFile_group = [outputFile(1:end-4) '_' options.timeseriesGroup{iGroup} ext];
    timeseriesData = cast(timeseriesData, dataType);
    if isempty(options.precision)
      writeBinary(timeseriesData, outputFile_group, writeBufferSize=options.writeBufferSize);
    else
      writeBinary(cast(timeseriesData, options.precision), ...
        outputFile_group, writeBufferSize=options.writeBufferSize);
    end
  end
end

% Save binary generation parameters (for replication)
procTime = datetime;
matFile = [outputFile_group(1:end-4) '.mat'];
save(matFile, 'inputFile', 'outputFile', 'options','procTime', '-v7.3');