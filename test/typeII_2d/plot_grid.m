%% Basic Settings.
clear;
close all;
warning off;

p = 5;

alpha = 2;
beta = 0.7;

%% Generate points x and omega.
n = 2^p;
nx = n;
ny = n;
N = nx * ny;

fprintf("Basic info:\n")
fprintf("nx: %d, ny: %d, N: %d\n", nx, ny, N);

fprintf("Rand grid:\n");
M = ceil(alpha * N);
x = rand(M, 2);
M = size(x, 1);
fprintf("M: %d, N: %d\n", M, N);
fprintf("M/N: %.1e\n", M / N);

% figure();
% plot(x(:, 1), x(:, 2), "Marker", ".", "LineStyle", "none");
% axis equal;
% xlim([0, 1]);
% ylim([0, 1]);
% saveas(gcf, './figure/rand_grid.png');
% saveas(gcf, './figure/rand_grid.eps', 'epsc');

fprintf("Polar grid:\n")
x = PolarGrid(n, 1, beta * p);
M = size(x, 1);
fprintf("M: %d, N: %d\n", M, N);
fprintf("M/N: %.1e\n", M / N);
fprintf("M/N/logN: %.1e\n", M / N / p);
M / N / p / beta

% figure();
% plot(x(:, 1), x(:, 2), "Marker", ".", "LineStyle", "none");
% axis equal;
% xlim([0, 1]);
% ylim([0, 1]);
% saveas(gcf, './figure/polar_grid.png');
% saveas(gcf, './figure/polar_grid.eps', 'epsc');

