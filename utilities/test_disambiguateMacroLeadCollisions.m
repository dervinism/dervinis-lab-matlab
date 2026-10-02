function test_disambiguateMacroLeadCollisions()
% test_disambiguateMacroLeadCollisions()
%
% Unit tests for disambiguateMacroLeadCollisions. Run from the MATLAB
% command window:
%   >> test_disambiguateMacroLeadCollisions
%
% Prerequisites (must be on the MATLAB path):
%   - bareLeadName.m (this repo's utilities folder)
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

w = warning('off', 'disambiguateMacroLeadCollisions:droppedA');
cleanup = onCleanup(@() warning(w));

% ---- 1. the MSEL02799 case: LT depth lead + LTMacro (=LTA) + LTAMicro -
names = [ ...
  compose("LT_%04d", 1:15), ...
  compose("LTMacro_%04d", 1:9), ...
  compose("LTAMicro_%04d", 1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 9, 'expected 9 rewritten channels');
assert(strcmp(m('LTMacro_0001'), 'LTAMacro_0001'), 'wrong correction');
assert(strcmp(m('LTMacro_0009'), 'LTAMacro_0009'), 'wrong correction');
% and the corrected names now reduce to LTA
assert(strcmp(bareLeadName(m('LTMacro_0005')), 'LTA'), 'corrected name must be bare LTA');
% LT_* and LTAMicro_* are untouched
assert(~isKey(m, 'LT_0001') && ~isKey(m, 'LTAMicro_0001'), 'only *Macro* names change');

% ---- 2. both sides, RT as well ----------------------------------
names = [compose("LT_%04d",1:15), compose("LTMacro_%04d",1:9), compose("LTAMicro_%04d",1:8), ...
         compose("RT_%04d",1:15), compose("RTMacro_%04d",1:9), compose("RTAMicro_%04d",1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 18, 'expected 18 rewritten channels (LT + RT)');
assert(strcmp(m('RTMacro_0003'), 'RTAMacro_0003'), 'RT side not corrected');

% ---- 3. RC: no separate depth lead, RCMicro already resolves -> no-op --
names = [compose("RCMacro_%04d",1:9), compose("RCMicro_%04d",1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 0, 'RCMacro must NOT be rewritten (no collision)');

% ---- 4. XMacro with a depth lead but NO XAMicro -> no-op ----------
names = [compose("LB_%04d",1:15), compose("LBMacro_%04d",1:6)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 0, 'no XAMicro sibling -> do not rewrite');

% ---- 5. XMacro with XAMicro but NO plain X depth lead -> no-op ------
names = [compose("PQMacro_%04d",1:9), compose("PQAMicro_%04d",1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 0, 'no plain X depth lead -> do not rewrite');

% ---- 6. lowercase / no-underscore token variants -----------------
names = [compose("LT%d",1:15), compose("LTmacro%d",1:9), compose("LTAmicro_%04d",1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 9 && strcmp(m('LTmacro3'), 'LTAmacro3'), 'lowercase token not handled');

% ---- 7. idempotency: already-correct names are a no-op ------------
names = [compose("LT_%04d",1:15), compose("LTAMacro_%04d",1:9), compose("LTAMicro_%04d",1:8)];
m = disambiguateMacroLeadCollisions(names);
assert(m.Count == 0, 'already-correct LTAMacro must not be touched');

% ---- 8. empty / cellstr input -----------------------------------
m = disambiguateMacroLeadCollisions({'EKG','LB_0001','LB_0002'});
assert(isa(m, 'containers.Map') && m.Count == 0, 'benign input must give empty map');

fprintf('test_disambiguateMacroLeadCollisions: all checks passed\n');
end
