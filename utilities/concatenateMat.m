function concatenatedMat = concatenateMat(mat1, mat2, method)
% concatenatedMat = concatenateMat(mat1, mat2, method)
%
% Function concatenates two matrices of any dimensions by appending the
% smaller matrix with trailing zeros or NaNs.
%
% Args:
%   mat1 (numeric, required, positional): a shape-(M, N) numeric array.
%   mat2 (numeric, required, positional): a shape-(K, L) numeric array.
%   method (numeric, optional, positional): a shape-(1, J) charachter array
%     defining the method of matrix concatentation. Available options are
%     'vertical' - vertical concatenation using trailing zeros (default).
%     'horizontal' - horizontal concatenation using trailing zeros.
%     'verticalnan' - vertical concatenation using trailing NaNs.
%     'horizontalnan' - horizontal concatenation using trailing NaNs.
%
% Returns:
%   concatenatedMat (numeric): a shape-(M+K, max([N L])) or
%     (max([M K]), N+L) numeric array of concatenated input matrices mat1
%     and mat2 with trailing zeros or NaNs, if needed.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

if nargin < 3
  method = 'vertical';
end

if ~isempty(mat1) && ~isempty(mat2)
  if strcmp(method, 'vertical')
    diff = size(mat1,2) - size(mat2,2);
    if diff > 0
      trailingMat = zeros(size(mat2,1), abs(diff));
      mat2 = [mat2 trailingMat];
    elseif diff < 0
      trailingMat = zeros(size(mat1,1), abs(diff));
      mat1 = [mat1 trailingMat];
    end
    concatenatedMat = [mat1; mat2];
  elseif strcmp(method, 'horizontal')
    diff = size(mat1,1) - size(mat2,1);
    if diff > 0
      trailingMat = zeros(abs(diff), size(mat2,2));
      mat2 = [mat2; trailingMat];
    elseif diff < 0
      trailingMat = zeros(abs(diff), size(mat1,2));
      mat1 = [mat1; trailingMat];
    end
    concatenatedMat = [mat1 mat2];
  elseif strcmp(method, 'verticalnan')
    diff = size(mat1,2) - size(mat2,2);
    if diff > 0
      trailingMat = nan(size(mat2,1), abs(diff));
      mat2 = [mat2 trailingMat];
    elseif diff < 0
      trailingMat = nan(size(mat1,1), abs(diff));
      mat1 = [mat1 trailingMat];
    end
    concatenatedMat = [mat1; mat2];
  elseif strcmp(method, 'horizontalnan')
    diff = size(mat1,1) - size(mat2,1);
    if diff > 0
      trailingMat = nan(abs(diff), size(mat2,2));
      mat2 = [mat2; trailingMat];
    elseif diff < 0
      trailingMat = nan(abs(diff), size(mat1,2));
      mat1 = [mat1; trailingMat];
    end
    concatenatedMat = [mat1 mat2];
  end
elseif ~isempty(mat1) && isempty(mat2)
  concatenatedMat = mat1;
elseif isempty(mat1) && ~isempty(mat2)
  concatenatedMat = mat2;
else
  concatenatedMat = [];
end