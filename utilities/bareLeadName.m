function [bareName, electrodeType] = bareLeadName(channelName)
% [bareName, electrodeType] = bareLeadName(channelName)
%
% Function reduces a recording channel or lead name to its bare lead name
% by removing the trailing contact number and the trailing macro/micro
% electrode-type token. The result is meant to be used as a comparison key
% only (e.g. to work out which macroelectrode lead a microelectrode wire
% belongs to). It must never be stored as a channel or lead identity
% because on its own it does not preserve the macro/micro distinction (that
% is what the second output is for).
%
% Args:
%   channelName (char, required, positional): a shape-(1, n) character
%     array containing the channel or lead name (e.g. 'LTAMicro_0001',
%     'LTmacro3', 'LT_0012', 'LTA').
%
% Returns:
%   bareName (char): a shape-(1, m) character array containing the bare
%     lead name with the trailing contact number and the trailing
%     electrode-type token removed (e.g. 'LTAMicro_0001' -> 'LTA',
%     'LTmacro3' -> 'LT', 'LT_0012' -> 'LT', 'LTA' -> 'LTA').
%   electrodeType (char): a shape-(1, p) character array naming the
%     electrode-type token that was removed from the end of the name, in
%     lower case ('macro', 'micro', 'shaft' or 'bundle'), or an empty
%     character array if the name carried no such token.
%
% Comments:
%   Electrode-type tokens are matched as whole words anchored at the end of
%   the string (case-insensitively), so a partial match that would confuse,
%   for instance, lead 'LT' with lead 'LTA' cannot occur. A single trailing
%   'extra' token (used to denote an additional microelectrode group
%   belonging to the same lead) is also removed. Only a single
%   electrode-type token is removed.
%
% Dependencies:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  channelName (1,:) {mustBeA(channelName,'char'), mustBeVector}
end

bareName = channelName;

% Remove the trailing contact number, with or without a preceding
% underscore (e.g. '_0001' or a bare '3').
bareName = regexprep(bareName, '_?\d+$', '', 'once');

% Remove a single trailing 'extra' token (an additional microelectrode
% group on the same lead), with or without a preceding underscore.
bareName = regexprep(bareName, '_?extra$', '', 'once', 'ignorecase');

% Remove a second trailing number that may sit in front of the
% electrode-type token (e.g. 'LTAmicro2').
bareName = regexprep(bareName, '\d+$', '', 'once');

% Remove a single trailing electrode-type token, matched as a whole word
% anchored at the end of the string, with or without a preceding
% underscore.
electrodeTypeToken = regexpi(bareName, '_?(macro|micro|shaft|bundle)$', 'tokens', 'once');
if isempty(electrodeTypeToken)
  electrodeType = '';
else
  electrodeType = lower(electrodeTypeToken{1});
  bareName = regexprep(bareName, '_?(macro|micro|shaft|bundle)$', '', 'once', 'ignorecase');
end

% Remove a dangling trailing underscore.
bareName = regexprep(bareName, '_$', '', 'once');
