function [dataMean, dataCI95, counts, dataStd, dataSEM] = datamean(data, meanType, options)
% [dataMean, dataCI95, counts, dataStd, dataSEM] = datamean(data, <meanType>)
%
% Function calculates the mean and the 95% confidence interval of a data
% matrix. Infinite values are treated as NaNs. NaNs are ommited from
% calculations.
%
% Args:
%   data (numeric, required, positional): a shape-(M, N) numeric array of
%     data values. If data is a matrix, the mean is calculated column-wise.
%   meanType (char, optional, positional): a shape-(1, K) character array
%     indicating the type of mean calculation. The following options are
%     available:
%       'regular' - a regular mean (Matlab's native; default).
%       'circular' - a circular parametric mean.
%       'circularNP' - a circular mean with non-parametric confidence
%                      intervals.
%   directional (logical, optional, keyword): a shape-(1, 1) logical scalar
%     indicating whether the confidence interval is unidirectional. In
%     other words, this optional variable indicates whether the test
%     hypothesis specifies the direction of the difference (e.g., is the
%     mean higher than some value). If true is specified the 90% confidence
%     interval is calculated instead of the 95% one. Default is false.
%
% Returns:
%   dataMean (numeric): a shape-(1, N) numeric array representing data
%     mean(s).
%   dataCI95 (numeric) a shape-(2, N) numeric array with 95% confidence
%     limits around the data mean(s).
%   counts (numeric): a shape-(1, N) numeric array of non-NaN counts (over
%     columns if data is a matrix).
%   dataStd (numeric) a shape-(1, N) numeric array with standard deviations.
%   dataSEM (numeric) a shape-(1, N) numeric array with standard error of
%     the mean.
%
% Dependencies:
%   Circular Statistics Toolbox (https://github.com/circstat/circstat-matlab).
%   Non-parametric circular statistics functions for Matlab (https://github.com/dervinism/circStatNP).
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com)

arguments
  data (:,:) {mustBeNumeric}
  meanType {mustBeMember(meanType,{'regular','circular','circularNP'})} = 'regular'
  options.directional (1,1) {mustBeA(options.directional,'logical')} = false
end

% Initialise output containers
dataMean = []; dataCI95 = [];
if isempty(data)
  return
end
data(isinf(data)) = NaN;
counts = sum(~isnan(data),1); % number of significant cells
F = size(data, 2); % number of frequencies or time
if counts == 1
  dataMean = data(~isnan(data));
  dataCI95 = 0; dataStd = 0; dataSEM = 0;
  return
end

% Mean calculations
if strcmp(meanType, 'regular') % regular means
  dataMean = mean(data, 1, 'omitnan');
  dataStd = std(data, 1, 'omitnan');
  dataSEM = dataStd ./ counts;
  CI95 = zeros(2,F); % regular confidence intervals
  dataCI95 = zeros(2,F);
  for f = 1:F
    if options.directional
      CI95(:,f) = (tinv([0.05 0.95], counts(f)-1))';
    else
      CI95(:,f) = (tinv([0.025 0.975], counts(f)-1))';
    end
    dataCI95(:,f) = bsxfun(@times, dataSEM(f), CI95(:,f));
  end
elseif strcmp(meanType, 'circular') || strcmp(meanType, 'circularNP') % circular means
  if options.directional
    alpha = 0.1;
  else
    alpha = 0.05;
  end
  dataMean = NaN(1,size(data,2));
  dataStd = NaN(1,size(data,2));
  dataSEM = NaN(1,size(data,2));
  dataCI95 = NaN(2,size(data,2));
  for j = 1:size(data,2)
    dataExist = ~isnan(data(:,j));
    if sum(dataExist)
      dataMean(j) = circmean(data(:,j));
      dataStd(j) = std(data(dataExist,j));
      dataSEM(j) = dataStd(j) ./ counts;
      if strcmp(meanType, 'circular')
        dataCI95(2,j) = circ_confmean(data(~isnan(data(:,j)),j), alpha); % parametric circular confidence intervals
      elseif strcmp(meanType, 'circularNP')
        dataCI95(2,j) = circ_confmeanFisher(data(~isnan(data(:,j)),j), alpha); % non-parametric circular confidence intervals
      end
    end
  end
  dataCI95(1,:) = -dataCI95(2,:);
end