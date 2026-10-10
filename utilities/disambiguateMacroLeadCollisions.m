function nameMap = disambiguateMacroLeadCollisions(chNames)
% nameMap = disambiguateMacroLeadCollisions(chNames)
%
% Function corrects '<stem>Macro' channel names that collide with a
% separate '<stem>' depth lead by reinstating a dropped 'A'.
%
% Behnke-Fried depth electrodes carry both macro ring contacts and
% microwire bundles, and their channels should share a bare lead name,
% e.g. 'LTAMacro_0001' (bare lead 'LTA') alongside 'LTAMicro_0001' (bare
% lead 'LTA'). A recording-time naming mistake sometimes drops the 'A'
% from the macro channels only, producing 'LTMacro_0001' - which
% bareLeadName reduces to 'LT', colliding with a genuine, separate 'LT'
% depth lead and leaving the 'LTA' microwire bundle with no macro sibling.
%
% This function detects that specific collision and returns a corrected
% name for the affected macro channels ('LTMacro_0001' ->
% 'LTAMacro_0001'). It is deliberately conservative: a '<stem>Macro'
% channel is rewritten only when BOTH
%   (1) a separate lead with bare name '<stem>' and no electrode-type
%       token (i.e. a plain depth lead) is present, and
%   (2) a micro/shaft/bundle channel with bare name '<stem>A' is present.
% '<stem>Macro' channels that do not collide - e.g. 'RCMacro_0001' when
% 'RC' is only ever the Behnke-Fried lead and 'RCMicro' already resolves
% to 'RC' - are left untouched.
%
% Only a single 'A' immediately before the 'macro' token is ever inserted;
% nothing else about the name changes, and the stored ChName should still
% be kept verbatim - use the returned mapping only to drive lead-name
% parsing.
%
% Args:
%   chNames (cellstr | string, required, positional): recording channel
%     names, verbatim.
%
% Returns:
%   nameMap (containers.Map): maps each rewritten original name to its
%     corrected name ('LTMacro_0001' -> 'LTAMacro_0001'). Empty
%     (containers.Map with no keys) when no collision is found.
%
% Comments:
%   One warning is emitted per corrected lead, naming the colliding depth
%   lead and the microwire bundle that motivated the correction.
%
% Dependencies:
%   bareLeadName.m (this repo's utilities folder).
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  chNames {mustBeVector}
end

chNames = cellstr(chNames);
chNames = chNames(:);
n = numel(chNames);

bareNames = cell(n, 1);
elecTypes = cell(n, 1);
for i = 1:n
  [bareNames{i}, elecTypes{i}] = bareLeadName(chNames{i});
end

isMacro  = strcmp(elecTypes, 'macro');
isPlain  = cellfun(@isempty, elecTypes);
isMicroLike = ismember(elecTypes, {'micro', 'shaft', 'bundle'});

plainLeadStems = unique(bareNames(isPlain & ~cellfun(@isempty, bareNames)));
microLeadStems = unique(bareNames(isMicroLike));
macroStems     = unique(bareNames(isMacro));

nameMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

for k = 1:numel(macroStems)
  stem = macroStems{k};
  collidesWithDepthLead = ismember(stem, plainLeadStems);
  hasOrphanMicroSibling = ismember([stem 'A'], microLeadStems);
  if collidesWithDepthLead && hasOrphanMicroSibling
    corrected = strings(0);
    for i = 1:n
      if isMacro(i) && strcmp(bareNames{i}, stem)
        fixed = local_insert_A_before_macro(chNames{i});
        nameMap(chNames{i}) = fixed;
        corrected(end+1) = string(chNames{i}); %#ok<AGROW>
      end
    end
    warning('disambiguateMacroLeadCollisions:droppedA', ...
      ['"%sMacro*" channels collide with depth lead "%s" and the "%sAMicro" ' ...
       'bundle has no macro sibling; reinterpreting them as "%sA" ' ...
       '(dropped-A naming mistake). Affected: %s'], ...
      stem, stem, stem, stem, strjoin(cellstr(corrected), ', '));
  end
end
end


function out = local_insert_A_before_macro(name)
% Insert a single 'A' immediately before the 'macro' token (and before any
% underscore that precedes it), leaving the rest of the name untouched:
%   'LTMacro_0001'  -> 'LTAMacro_0001'
%   'LTmacro3'      -> 'LTAmacro3'
%   'LT_macro_0001' -> 'LTA_macro_0001'
tokStart = regexp(name, '_?[Mm]acro', 'once');
if isempty(tokStart)
  out = name;   % nothing to do (should not happen for an etype 'macro' name)
else
  out = [name(1:tokStart-1) 'A' name(tokStart:end)];
end
end
