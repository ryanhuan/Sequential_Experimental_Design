clear all
close all

set(0,'DefaultAxesLineWidth',2);
set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times');

% Specialized to 2D states only right now.
nStages = 3;
pOrder = 3;

gaussianPCWeightParams = [7.0 9.0]';
uniformPCWeightParams = [1e-5 9.0]';

coefsAll = {};
for i = 0:nStages
  fNameNoExt = strcat('coefs', num2str(i));
  fNameYesExt = strcat(fNameNoExt, '.dat');
  load(fNameYesExt);
  coefsAll{end + 1} = eval(fNameNoExt);
end
load refTable.dat

x = -5:0.01:5;
y = -1:0.01:1;

[X Y] = meshgrid(x,y);
mu = gaussianPCWeightParams(1) + sqrt(gaussianPCWeightParams(2)) .* X;
sigmaSq = Y .* 0.5 .* (uniformPCWeightParams(2)-uniformPCWeightParams(1)) + ...
          (uniformPCWeightParams(2)+uniformPCWeightParams(1))/2;

% Evaluate Hermite and Legendre polynomials on all grid points up to the maximum order.
XHermite = evalHermite(X, pOrder);
YLegendre = evalLegendre(Y, pOrder);

for i = 0:nStages
  Z = 0.* X;
  
  for j = 1:length(coefsAll{1})
    % Note 0th order is at index 1, 1st order at index 2, etc, since MATLAB starts from
    % 1 rather than 0. So we do refTable(j, ) + 1.
    Z = Z + coefsAll{i + 1}(j) * XHermite(:, :, refTable(j, 1) + 1) .* YLegendre(:, :, refTable(j, 2) + 1);
  end
  
  figure(i + 1)
  surf(mu, sigmaSq, Z, 'edgecolor', 'none')
  xlabelString = strcat('\mu_', num2str(i));
  ylabelString = strcat('\sigma_', num2str(i), '^2');
  zlabelString = strcat('Tilde J_', num2str(i));
  xlabel(xlabelString)
  ylabel(ylabelString)
  zlabel(zlabelString)
  axis([-8 22 0 10 0 20])

  figureName = strcat('tildeJ', num2str(i));
  fID = strcat('-f', num2str(i + 1));
%  print(fID, '-depsc', figureName)

end
