function test_bareLeadName()
% test_bareLeadName()
%
% Unit tests for bareLeadName. Run from the MATLAB command window:
%   >> test_bareLeadName
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

cases = {
% channelName            expected bareName   expected electrodeType
  'LT_0001',             'LT',               ''
  'LT_0012',             'LT',               ''
  'LT1',                 'LT',               ''
  'LTA',                 'LTA',              ''
  'LTA_0009',            'LTA',              ''
  'LTMacro_0001',        'LT',               'macro'
  'LTmacro_0003',        'LT',               'macro'
  'LTMacro3',            'LT',               'macro'
  'LTAMacro_0009',       'LTA',              'macro'
  'LTAmacro9',           'LTA',              'macro'
  'LTMicro_0001',        'LT',               'micro'
  'LTmicro_0002',        'LT',               'micro'
  'LTAMicro_0001',       'LTA',              'micro'
  'LTAmicro_0008',       'LTA',              'micro'
  'LTA_micro_0001',      'LTA',              'micro'
  'RCMicro_0004',        'RC',               'micro'
  'RCMicroExtra_0001',   'RC',               'micro'
  'LBmicroExtra_0005',   'LB',               'micro'
  'LBshaft_0001',        'LB',               'shaft'
  'LBbundle_0001',       'LB',               'bundle'
  'LTAmicro2_0001',      'LTA',              'micro'
  'EKG',                 'EKG',              ''
  };

nPass = 0; nFail = 0;
for i = 1:size(cases,1)
  [bareName, electrodeType] = bareLeadName(cases{i,1});
  ok = strcmp(bareName, cases{i,2}) && strcmp(electrodeType, cases{i,3});
  if ok
    nPass = nPass + 1;
  else
    nFail = nFail + 1;
    fprintf('  FAIL: %-20s -> ("%s","%s"), expected ("%s","%s")\n', ...
      cases{i,1}, bareName, electrodeType, cases{i,2}, cases{i,3});
  end
end

% Critical property: LT and LTA must never collapse into each other.
assert(~strcmp(bareLeadName('LTMicro_0001'), bareLeadName('LTAMicro_0001')), ...
  'LT and LTA bare names collided');

fprintf('test_bareLeadName: %d passed, %d failed\n', nPass, nFail);
assert(nFail == 0, 'test_bareLeadName had failures');
