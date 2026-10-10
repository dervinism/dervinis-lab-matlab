function test_microwireContactOffset()
% test_microwireContactOffset()
%
% Unit tests for microwireContactOffset. Run from the MATLAB command
% window:
%   >> test_microwireContactOffset
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

nTetrodes = 2;
wireSpacing = 0.035;
interTetrodeSpacing = 2;

% Collect the 8 offsets of a 2-tetrode bundle.
offsets = zeros(nTetrodes*4, 3);
row = 0;
for iTetrode = 1:nTetrodes
  for iContact = 1:4
    row = row + 1;
    offsets(row,:) = microwireContactOffset(iTetrode, nTetrodes, iContact);
  end
end

% 1. All 8 contacts are distinct.
assert(size(unique(offsets, 'rows'), 1) == 8, 'contact offsets are not all distinct');

% 2. Every tetrode's 4 wires trace a square of side wireSpacing.
for iTetrode = 1:nTetrodes
  tet = offsets((iTetrode-1)*4 + (1:4), :);
  sideLengths = sqrt(sum(diff(tet([1 2 3 4 1],:)).^2, 2));
  assert(all(abs(sideLengths - wireSpacing) < 1e-12), ...
    sprintf('tetrode %d is not a %g mm square', iTetrode, wireSpacing));
  assert(all(tet(:,3) == 0), 'tetrode square is not in the x-y plane');
end

% 3. Tetrode centres are interTetrodeSpacing apart, along y only.
c1 = mean(offsets(1:4,:), 1);
c2 = mean(offsets(5:8,:), 1);
assert(abs((c2(2) - c1(2)) - interTetrodeSpacing) < 1e-12, 'inter-tetrode spacing wrong');
assert(abs(c2(1) - c1(1)) < 1e-12 && abs(c2(3) - c1(3)) < 1e-12, ...
  'tetrodes offset outside the y axis');

% 4. First tetrode sits at the lead tip (square corner 1 has zero offset).
assert(isequal(microwireContactOffset(1, nTetrodes, 1), [0 0 0]), ...
  'first tetrode is not anchored at the lead tip');

% 5. Options are honoured.
o = microwireContactOffset(2, 2, 3, 'interTetrodeSpacing', 3, 'wireSpacing', 0.05);
assert(abs(o(1) - 0.05) < 1e-12, 'wireSpacing option ignored');
assert(abs(o(2) - (3 + 0.05)) < 1e-12, 'interTetrodeSpacing option ignored');

fprintf('test_microwireContactOffset: all checks passed\n');
