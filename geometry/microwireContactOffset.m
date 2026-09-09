function contactOffset = microwireContactOffset(tetrodeNumber, nTetrodes, contactNumber, options)
% contactOffset = microwireContactOffset(tetrodeNumber, nTetrodes, contactNumber, <options>)
%
% Function calculates the 3D coordinate offset of an individual
% microelectrode wire within a tetrode bundle relative to the tip of the
% associated macroelectrode lead. Adding this offset to the lead tip
% coordinate yields an approximate and, crucially, unique coordinate for
% every wire in the bundle.
%
% Args:
%   tetrodeNumber (numeric, required, positional): a shape-(1, 1) numeric
%     scalar; the 1-based index of the tetrode within the bundle.
%   nTetrodes (numeric, required, positional): a shape-(1, 1) numeric
%     scalar; the total number of tetrodes in the bundle.
%   contactNumber (numeric, required, positional): a shape-(1, 1) numeric
%     scalar; the 1-based index of the wire within the tetrode (1-4).
%   interTetrodeSpacing (numeric, optional, keyword): a shape-(1, 1)
%     numeric scalar; the centre-to-centre distance between consecutive
%     tetrodes along the bundle axis, in millimetres (default=2).
%   wireSpacing (numeric, optional, keyword): a shape-(1, 1) numeric
%     scalar; the side length of the square traced by the four wires of a
%     single tetrode, in millimetres (default=0.035).
%
% Returns:
%   contactOffset (numeric): a shape-(1, 3) numeric array with the
%     [dx, dy, dz] coordinate offset in millimetres to be added to the
%     associated macroelectrode lead tip coordinate.
%
% Comments:
%   The true positions of individual microwires within a bundle are not
%   recoverable, so this is deliberately simple placeholder geometry: the
%   four wires of a tetrode are placed on the corners of a square in the
%   x-y plane and successive tetrodes are separated along the y axis,
%   centred on the associated lead tip. The only property it guarantees is
%   that every wire of a bundle receives a distinct, sub-millimetre-scale
%   coordinate clustered at the lead tip.
%
% Dependencies:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  tetrodeNumber (1,1) {mustBeInteger, mustBePositive}
  nTetrodes (1,1) {mustBeInteger, mustBePositive}
  contactNumber (1,1) {mustBeInteger, mustBeInRange(contactNumber, 1, 4)}
  options.interTetrodeSpacing (1,1) {mustBeNumeric, mustBeReal} = 2
  options.wireSpacing (1,1) {mustBeNumeric, mustBeReal} = 0.035
end

% The four wires of a tetrode sit on the corners of a square in the x-y
% plane (column 1 -> dx, column 2 -> dy).
squareCorners = options.wireSpacing.*[0 0; 1 0; 1 1; 0 1];

% Successive tetrodes are offset along the y axis and centred on the
% associated lead tip.
interTetrodeOffset = options.interTetrodeSpacing.*(tetrodeNumber - 0.5.*nTetrodes);

contactOffset = [squareCorners(contactNumber,1), ...
  interTetrodeOffset + squareCorners(contactNumber,2), 0];
