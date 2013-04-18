function [val] = evalLegendre(X, pOrder)

[m, n] = size(X);
val = zeros(m, n, pOrder + 1);

for i = 1:n
  val(:, i, 1) = ones(m, 1);
  val(:, i, 2) = X(:, i);
  for j = 3:pOrder + 1
    val(:, i, j) = ((2 * j - 1) .* X(:, i) .* val(:, i, j - 1) - (j - 1) .* val(:, i, j ...
                                                      - 2)) / j;
  end
end
