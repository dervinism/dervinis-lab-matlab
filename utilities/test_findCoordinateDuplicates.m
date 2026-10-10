function test_findCoordinateDuplicates()
% test_findCoordinateDuplicates()
%
% Unit tests for findCoordinateDuplicates. Run from the MATLAB command
% window:
%   >> test_findCoordinateDuplicates
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

% ---- 1. no duplicates -> empty 0x1 cell ------------------------------
d = findCoordinateDuplicates(["A";"B";"C"], [0 0 0; 1 1 1; 2 2 2]);
assert(iscell(d) && isempty(d) && isequal(size(d), [0 1]), ...
  'expected empty 0x1 cell when there are no duplicates');

% ---- 2. one exact duplicate pair ------------------------------------
d = findCoordinateDuplicates(["RB1";"RB2";"RB3"], ...
  [16.15 -10.44 -13.52; 16.15 -10.44 -13.52; 20 -10 -13]);
assert(isscalar(d), 'expected exactly one duplicate group');
assert(isequal(sort(d{1}.names), ["RB1";"RB2"]), 'wrong members');
assert(isequal(d{1}.xyz, [16.15 -10.44 -13.52]), 'wrong shared coordinate');

% ---- 3. a triple ---------------------------------------------------
d = findCoordinateDuplicates(["M1";"M2";"M3";"M4"], ...
  [0 0 0; 0 0 0; 0 0 0; 5 5 5]);
assert(isscalar(d) && numel(d{1}.names) == 3, 'expected one group of three');

% ---- 4. non-finite rows ignored ----------------------------------
d = findCoordinateDuplicates(["A";"B";"C"], [NaN 0 0; NaN 0 0; 1 1 1]);
assert(isempty(d), 'NaN rows must not be reported as duplicates');

% ---- 5. tolerance groups near-coincident points -----------------
p = [10 10 10; 10.02 9.99 10.01; 30 30 30];
assert(isempty(findCoordinateDuplicates(["A";"B";"C"], p)), ...
  'exact comparison must not group near points');
d = findCoordinateDuplicates(["A";"B";"C"], p, 'Tolerance', 0.1);
assert(isscalar(d) && isequal(sort(d{1}.names), ["A";"B"]), ...
  'tolerance must group points within it');

% ---- 6. size mismatch -> friendly error ------------------------
threw = false;
try
  findCoordinateDuplicates(["A";"B"], [0 0 0]);
catch e
  threw = true;
  assert(strcmp(e.identifier, 'findCoordinateDuplicates:sizeMismatch'), 'wrong id');
end
assert(threw, 'size mismatch must error');

fprintf('test_findCoordinateDuplicates: all checks passed\n');
end
