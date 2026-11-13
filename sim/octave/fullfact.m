function d = fullfact(levels)
% Minimal stand-in for the Statistics Toolbox fullfact, for Octave runs only.
n = prod(levels); k = numel(levels); d = zeros(n, k); rep = 1;
for j = 1:k
  col = repmat(kron((1:levels(j))', ones(rep, 1)), n / (rep * levels(j)), 1);
  d(:, j) = col; rep = rep * levels(j);
end
end
