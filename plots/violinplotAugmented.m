function violins = violinplotAugmented(scatterGroups, groupNamesOrIDs, options)
% [fH, violins] = violinplotAugmented(scatterGroups, groupNamesOrIDs, <options>)
%
% Function produces a violin plot. It is a wrapper function for
% violinplot() with augmented functionality to display data means and
% confidence limits.
%
% Args:
%   scatterGroups (cell | num | struct, required, positional): a
%     shape-(1, n) or -(n, 1) cell array of data vectors, a shape (m, n)
%     numeric array of data column vectors, or a shape-(1, 1) scalar
%     structure with fields corresponding to data vectors (the latter two
%     types are inherited from violinplot(); can also take in data
%     structured as a table or a dataset). Violins will be plotted for
%     each data vector (a table column or a data set).
%   groupNamesOrIDs (cell | num, required, positional): a shape-(1, n) or
%     -(n, 1) cell array of data vector names (corresponds to individual
%     data cells/vectors/columns/sets in scatterGroups) or a shape-(1, l)
%     array of category IDs (scatterGroups must be a shape-(1, l) numeric
%     array; inherited from violinplot()).
%   options.dataMeans (num, optional, keyword): a shape-(1, n) numeric
%     array of data means corresponding to individual data
%     cells/vectors/columns/sets in scatterGroups. If not supplied or left
%     empty, data mean lines will not be marked on the violins (default).
%     Otherwise overrides the default ShowMean option of violinplot().
%   options.dataCIs (num, optional, keyword): a shape-(2, n) numeric array
%     of data mean confidence limits (the lower limit in row 1 and the
%     upper one in row 2; could be 95% confidence limits) corresponding to
%     individual data cells/vectors/columns/sets in scatterGroups (columns
%     in options.dataCIs). If not supplied or left empty, confidence limit
%     dashed lines will not be marked on the violins (default).
%   yScale (char, optional, keyword): a shape-(1, k) character array
%     describing the type of y axis scale. Available options are 'log' and
%     'regular' (default).
%   edgeVisibility (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for keeping (true; default) or eliminating (false) violin
%     contours (overrides the default EdgeColor option of violinplot()).
%   medianPlot (logical, optional, keyword): a shape-(1, 1) logical
%     scalar for displaying the median plot (true; default) or not (false).
%   fillColours (cell, optional, keyword): a shape-(1, n) or -(n, 1) cell
%     array of violin fill colour triplets corresponding to individual data
%     cells/vectors/columns/sets in scatterGroups. If not specified or left
%     empty, will have no effect. This options overrides ViolinColor option
%     of violinplot().
%   edgeColours (cell, optional, keyword): a shape-(1, n) or -(n, 1) cell
%     array of violin edge colour triplets corresponding to individual data
%     cells/vectors/columns/sets in scatterGroups. If not specified or left
%     empty, will have no effect. This options overrides EdgeColor option
%     of violinplot().
%   markerFaceAlpha (num, optional, keyword): a shape-(1, 1) numeric scalar
%     describing the opacity of the data points. 1 means completely opaque
%     while 0 means completely transparent.
%   Other optional keyword arguments are the same as for the violinplot()
%     function (type help violinplot for their descriptions, types,
%     dimensions, available options, and default values) and include Width,
%     Bandwidth, ViolinColor, ViolinAlpha, EdgeColor, BoxColor,
%     MedianColor, ShowData, ShowNotches, ShowMean, GroupOrder.
%
% Returns:
%   violins (violin): an object array of
%     <a href="matlab:help('Violin')">Violin</a> objects.
%
% Dependencies:
%   dervinis-lab-matlab
%     (https://github.com/dervinism/dervinis-lab-matlab/).
%
% Author:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  scatterGroups {mustBeNonempty}
  groupNamesOrIDs {mustBeNonempty}
  options.dataMeans (1,:) {mustBeVector,mustBeNumeric} = []
  options.dataCIs (2,:) {mustBeNumeric} = []
  options.yScale (1,:) {mustBeMember(options.yScale,{'log','regular'})} = 'regular'
  options.edgeVisibility (1,:) {islogical} = true
  options.medianPlot (1,:) {islogical} = true
  options.fillColours {mustBeA(options.fillColours,'cell')} = {}
  options.edgeColours {mustBeA(options.edgeColours,'cell')} = {}
  options.markerFaceAlpha (1,1) {mustBeInRange(options.markerFaceAlpha, 0, 1)} = 1
  options.Width (1,1) {mustBeNonnegative} = 0.3
  options.Bandwidth {mustBeScalarOrEmpty} = []
  options.ViolinColor {mustBeVector(options.ViolinColor,'allow-all-empties'),mustBeNonnegative} = []
  options.ViolinAlpha (1,1) {mustBeNonnegative} = 0.3
  options.EdgeColor (1,3) {mustBeNonnegative} = [0.5 0.5 0.5]
  options.BoxColor (1,3) {mustBeNonnegative} = [0.5 0.5 0.5]
  options.MedianColor (1,3) {mustBeNonnegative} = [1 1 1]
  options.ShowData (1,:) {islogical} = true
  options.ShowNotches (1,:) {islogical} = false
  options.ShowMean (1,:) {islogical} = false
  options.GroupOrder (1,:) {mustBeA(options.GroupOrder,'cell')} = {}
end

% Parse input
if ~isempty(options.dataMeans)
  options.ShowMean = true;
end
if ~isempty(options.ViolinColor) && size(options.ViolinColor,2) ~= 3
  error('options.ViolinColor must be a numeric array with the number of columns equal to 3.');
end
if iscell(scatterGroups)
  for iViolin = 1:numel(scatterGroups)
    scatterGroupsStruct.(['Violin' num2str(iViolin)]) = scatterGroups{iViolin};
  end
  scatterGroups = scatterGroupsStruct;
end

% Plot violins
if isempty(options.Bandwidth) && isempty(options.ViolinColor)
  violins = violinplot(scatterGroups, groupNamesOrIDs, ...
    'Width',options.Width, ...
    'ViolinAlpha',options.ViolinAlpha, ...
    'EdgeColor',options.EdgeColor, ...
    'BoxColor',options.BoxColor, ...
    'MedianColor',options.MedianColor, ...
    'ShowData',options.ShowData, ...
    'ShowNotches',options.ShowNotches, ...
    'ShowMean',options.ShowMean, ...
    'GroupOrder',options.GroupOrder);
elseif isempty(options.Bandwidth)
  violins = violinplot(scatterGroups, groupNamesOrIDs, ...
    'Width',options.Width, ...
    'ViolinColor',options.ViolinColor, ...
    'ViolinAlpha',options.ViolinAlpha, ...
    'EdgeColor',options.EdgeColor, ...
    'BoxColor',options.BoxColor, ...
    'MedianColor',options.MedianColor, ...
    'ShowData',options.ShowData, ...
    'ShowNotches',options.ShowNotches, ...
    'ShowMean',options.ShowMean, ...
    'GroupOrder',options.GroupOrder);
elseif isempty(options.ViolinColor)
  violins = violinplot(scatterGroups, groupNamesOrIDs, ...
    'Width',options.Width, ...
    'Bandwidth',options.Bandwidth, ...
    'ViolinAlpha',options.ViolinAlpha, ...
    'EdgeColor',options.EdgeColor, ...
    'BoxColor',options.BoxColor, ...
    'MedianColor',options.MedianColor, ...
    'ShowData',options.ShowData, ...
    'ShowNotches',options.ShowNotches, ...
    'ShowMean',options.ShowMean, ...
    'GroupOrder',options.GroupOrder);
else
  violins = violinplot(scatterGroups, groupNamesOrIDs, ...
    'Width',options.Width, ...
    'Bandwidth',options.Bandwidth, ...
    'ViolinColor',options.ViolinColor, ...
    'ViolinAlpha',options.ViolinAlpha, ...
    'EdgeColor',options.EdgeColor, ...
    'BoxColor',options.BoxColor, ...
    'MedianColor',options.MedianColor, ...
    'ShowData',options.ShowData, ...
    'ShowNotches',options.ShowNotches, ...
    'ShowMean',options.ShowMean, ...
    'GroupOrder',options.GroupOrder);
end

% Adjust violin properties
nViolins = numel(violins);
for iViolin = 1:nViolins
  if ~isempty(options.fillColours)
    violins(iViolin).ViolinColor = options.fillColours{iViolin};
  end
  if ~isempty(options.edgeColours)
    violins(iViolin).EdgeColor = options.edgeColours{iViolin};
  end
  if options.medianPlot
    violins(iViolin).MedianColor = [1 1 1];
  else
    violins(iViolin).MedianColor = [0.5 0.5 0.5];
    violins(iViolin).MedianPlot.SizeData = 1;
  end
  violins(iViolin).ViolinAlpha = options.ViolinAlpha;
  violins(iViolin).ShowNotches = options.ShowNotches;
  violins(iViolin).ShowMean = options.ShowMean;
  violins(iViolin).ScatterPlot.MarkerFaceAlpha = options.markerFaceAlpha;
  if ~options.edgeVisibility
    violins(iViolin).ViolinPlot.Visible  = 'on';
    violins(iViolin).BoxPlot.Visible  = 'off';
    violins(iViolin).WhiskerPlot.Visible  = 'off';
  end

  % Produce 95% confidene intervals on means
  if ~isempty(options.dataCIs)
    sideSizeY = round(size(violins(iViolin).ViolinPlot.Vertices,1)/2);
    if sideSizeY < 2
      continue
    end
    rightSideY = violins(iViolin).ViolinPlot.Vertices(1:sideSizeY,2);
    interpRightSideY = sort([rightSideY; options.dataMeans(iViolin)+options.dataCIs(:,iViolin)], 'ascend');
    leftSideY = violins(iViolin).ViolinPlot.Vertices(sideSizeY+1:end,2);
    interpLeftSideY = sort([leftSideY; options.dataMeans(iViolin)+options.dataCIs(:,iViolin)], 'descend');

    sideSizeX = round(size(violins(iViolin).ViolinPlot.Vertices,1)/2);
    if sideSizeX < 2
      continue
    end
    rightSideX = violins(iViolin).ViolinPlot.Vertices(1:sideSizeX,1);
    interpRightSideX = interp1(rightSideY, rightSideX, interpRightSideY);
    leftSideX = violins(iViolin).ViolinPlot.Vertices(sideSizeX+1:end,1);
    interpLeftSideX = interp1(leftSideY, leftSideX, interpLeftSideY);

    lowerLimit1 = interpLeftSideX(interpLeftSideY == options.dataMeans(iViolin)+options.dataCIs(1,iViolin));
    lowerLimit2 = interpRightSideX(interpRightSideY == options.dataMeans(iViolin)+options.dataCIs(1,iViolin));
    upperLimit1 = interpLeftSideX(interpLeftSideY == options.dataMeans(iViolin)+options.dataCIs(2,iViolin));
    upperLimit2 = interpRightSideX(interpRightSideY == options.dataMeans(iViolin)+options.dataCIs(2,iViolin));

    if ~isempty(options.fillColours)
      plot([lowerLimit1 lowerLimit2], [options.dataMeans(iViolin)+options.dataCIs(1,iViolin) options.dataMeans(iViolin)+options.dataCIs(1,iViolin)], '--', 'Color',options.fillColours{iViolin});
      plot([upperLimit1 upperLimit2], [options.dataMeans(iViolin)+options.dataCIs(2,iViolin) options.dataMeans(iViolin)+options.dataCIs(2,iViolin)], '--', 'Color',options.fillColours{iViolin});
    else
      plot([lowerLimit1 lowerLimit2], [options.dataMeans(iViolin)+options.dataCIs(1,iViolin) options.dataMeans(iViolin)+options.dataCIs(1,iViolin)], '--', 'Color',violins(iViolin).ViolinColor);
      plot([upperLimit1 upperLimit2], [options.dataMeans(iViolin)+options.dataCIs(2,iViolin) options.dataMeans(iViolin)+options.dataCIs(2,iViolin)], '--', 'Color',violins(iViolin).ViolinColor);
    end
  end
end

% Adjust figure properties
if strcmp(options.yScale,'log')
  set(gca,'Yscale','log')
end

xticklabels(groupNamesOrIDs);