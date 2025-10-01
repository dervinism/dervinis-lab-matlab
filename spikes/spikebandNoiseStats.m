function stats = spikebandNoiseStats(x, options)
% stats = spikebandNoiseStats(x, <options>)
% Functions calculates robust noise metrics for extracellular data.
%
% Args:
%   x (numeric, required, positional): a shape-(1, N) numeric array with an
%     LFP trace in microvolts (µV).
%   spikeTimes (numeric, optional, keyword): a shape-(1, M) numeric array
%     with spike times in seconds. If spike times are not provided, then
%     they will not be excluded from noise estimate calculations (not
%     ideal).
%   sr (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing sampling rate in Hz (default is 32000).
%   robustPct (numeric, optional, keyword): a shape-(1, 2) numeric array of
%     percentiles for robust peak-to-peak noise estimate (default =
%     [0.1 99.9]%).
%   winSec (numeric, optional, keyword): a shape-(1, 1) numeric scalar
%     representing window size in seconds for a windowed peak-to-peak noise
%     estimate (default = 2).
%
% Returns:
%   stats (struct): a shape-(1, 1) scalar structure with the following
%     fields:
%     rms - root mean square estimate prior to spike removal.
%     p2p - peak-to-peak noise estimate prior to spike removal.
%     rms_robust - root mean square estimate after spike removal (if spike
%                  times were supplied).
%     p2p_robust - peak-to-peak noise estimate after spike removal (if
%                  spike times were supplied).
%     sigma_MAD - robust σ estimate from MAD.
%     pct - percentiles used.
%
% Dependencies:
%   None.
%
% Comments:
%   For peak-to-peak, we report both naive (max-min) and robust (percentile
%   span), since single artifacts can inflate naive p2p.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  x (1,:) {mustBeNumeric,mustBeVector}
  options.spikeTimes (1,:) {mustBePositive,mustBeVector} = []
  options.sr (1,1) {mustBePositive} = 32000
  options.robustPct (1,:) {mustBePositive,mustBeVector} = [0.1 99.9]
  options.winSec (1,1) {mustBePositive} = 2
end

% Parse the input
x = x(:);

% ---------- 1) Naive (spikes included) metrics
stats.rms = rms_centered(x);
stats.p2p = max(x) - min(x);

% ---------- 2) Robust (spikes rejected) metrics
% Mask out spike samples (and a little margin)
spikeIdx = round(options.spikeTimes/options.sr);
spikeMask = false(numel(x),1);
for iSpike = 1:numel(spikeIdx)
  spikeStartInd = max([1 spikeIdx(iSpike) - round(0.5*options.sr)]);
  spikeEndInd = min([spikeIdx(iSpike) + round(0.5*options.sr) numel(x)]);
  spikeMask(spikeStartInd:spikeEndInd) = true;
end
x_clean = x(~spikeMask);

% Robust RMS on spike-free data
stats.rms_robust = rms_centered(x_clean);

% Robust p2p via percentiles (resistant to outliers)
lo = prctile(x_clean, options.robustPct(1));
hi = prctile(x_clean, options.robustPct(2));
stats.p2p_robust = hi - lo;
stats.pct = options.robustPct;

% ---------- 3) Optional: windowed peak-to-peak (for stability checks)
if ~isempty(options.winSec) && options.winSec > 0
  Nw = max(1, round(options.winSec * options.sr));
  n  = numel(x_clean);
  p2p_win = [];
  for i = 1:Nw: (n - Nw + 1)
    seg = x_clean(i:i+Nw-1);
    % robust window p2p
    p2p_win(end+1,1) = prctile(seg,options.robustPct(2)) - prctile(seg,options.robustPct(1)); %#ok<AGROW>
  end
  stats.p2p_window_mean = mean(p2p_win(p2p_win > 0));
  stats.p2p_window_median = median(p2p_win(p2p_win > 0));
  stats.p2p_window_std = std(p2p_win(p2p_win > 0));
end
end

% ---------- helpers ----------
function r = rms_centered(x)
x = x - mean(x);
r = sqrt(mean(x.^2));
end
