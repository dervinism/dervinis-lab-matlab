function timeSeries = highpassFilterTimeSeries(timeSeries, options)
% filtTimeSeries = highpassFilterTimeSeries(timeSeries, <options>)
%
% High-pass filter a row matrix containing time series vectors.
%
% Args:
%   timeSeries (numeric, required, positional): a shape-(N, M) numeric
%     array of continuous time series vectors (rows).
%   sampleRate (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%    corresponding to sampling frequency in Hz (default=500).
%   frequencyCutoff (numeric, optional, keyword): a shape-(1, 1) numeric
%    scalar with the frequency cutoff value.
%
% Returns:
%   filtTimeSeries (numeric): a shape-(N, M) numeric array of high-pass
%     filtered time series vectors (rows). Rows containing only NaNs are
%     skipped while rows containing some NaNs elicit an error.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  timeSeries (:,:) {mustBeNumeric}
  options.sampleRate (1,1) {mustBeNumeric,mustBePositive} = 500
  options.frequencyCutoff (1,1) {mustBeNumeric,mustBePositive} = 12
end

% Filter time series data
Wn_theta = options.frequencyCutoff(1)/(options.sampleRate/2); % High-pass frequency normalised by the Nyquist frequency
[btheta,atheta] = butter(3,Wn_theta, 'high'); % Apply butterworth filter
for row = 1:size(timeSeries,1)
  if sum(isnan(timeSeries(row,:))) && sum(isnan(timeSeries(row,:))) < numel(timeSeries(row,:))
    error('Input data contains NaNs');
  elseif ~sum(isnan(timeSeries(row,:)))
    timeSeries(row,:) = filtfilt(btheta,atheta,double(timeSeries(row,:)));
  end
end