%% MATLAB code: Red square positioned at upper-right, smaller size, 
% other regions partitioned accordingly
% Save as block_matrix_offset.m

clear; clc; close all;

%% Setup figure
figure('Color', 'white', 'Position', [100 100 800 800]);
hold on;
axis equal;
axis off;

%% Define colors
colors.yellow = [0.85, 0.65, 0.13];  % Golden yellow
colors.green = [0.18, 0.55, 0.34];   % Sea green
colors.blue = [0.27, 0.51, 0.71];    % Steel blue
colors.red = [0.80, 0.36, 0.36];     % Indian red

%% Define red square parameters (offset to upper-right, smaller size)
red_w = 0.15;   % Red square width (smaller)
red_h = 0.10;   % Red square height
red_x = 0.65;   % Positioned to the right
red_y = 0.55;   % Positioned upward

%% Define partition lines based on red square position
% Vertical partition lines: x = red_x and x = red_x + red_w
% Horizontal partition lines: y = red_y and y = red_y + red_h

%% Define 9 regions
regions = struct();

% Row 1 (Upper)
regions(1).name = 'A';  % Upper-left
regions(1).x = 0; regions(1).y = red_y + red_h; 
regions(1).w = red_x; regions(1).h = 1 - (red_y + red_h);
regions(1).color = colors.yellow;

regions(2).name = 'B';  % Upper-middle (above red square)
regions(2).x = red_x; regions(2).y = red_y + red_h;
regions(2).w = red_w; regions(2).h = 1 - (red_y + red_h);
regions(2).color = colors.green;

regions(3).name = 'C';  % Upper-right
regions(3).x = red_x + red_w; regions(3).y = red_y + red_h;
regions(3).w = 1 - (red_x + red_w); regions(3).h = 1 - (red_y + red_h);
regions(3).color = colors.yellow;

% Row 2 (Middle, containing red square)
regions(4).name = 'D';  % Middle-left (left of red square)
regions(4).x = 0; regions(4).y = red_y;
regions(4).w = red_x; regions(4).h = red_h;
regions(4).color = colors.blue;

regions(5).name = 'E';  % Center - Red square
regions(5).x = red_x; regions(5).y = red_y;
regions(5).w = red_w; regions(5).h = red_h;
regions(5).color = colors.red;

regions(6).name = 'F';  % Middle-right (right of red square)
regions(6).x = red_x + red_w; regions(6).y = red_y;
regions(6).w = 1 - (red_x + red_w); regions(6).h = red_h;
regions(6).color = colors.blue;

% Row 3 (Lower)
regions(7).name = 'G';  % Lower-left
regions(7).x = 0; regions(7).y = 0;
regions(7).w = red_x; regions(7).h = red_y;
regions(7).color = colors.yellow;

regions(8).name = 'H';  % Lower-middle (below red square)
regions(8).x = red_x; regions(8).y = 0;
regions(8).w = red_w; regions(8).h = red_y;
regions(8).color = colors.green;

regions(9).name = 'I';  % Lower-right
regions(9).x = red_x + red_w; regions(9).y = 0;
regions(9).w = 1 - (red_x + red_w); regions(9).h = red_y;
regions(9).color = colors.yellow;

%% Draw each region (solid fill)
for i = 1:length(regions)
    r = regions(i);
    % Draw filled rectangle
    fill([r.x, r.x+r.w, r.x+r.w, r.x], [r.y, r.y, r.y+r.h, r.y+r.h], ...
         r.color, 'EdgeColor', 'k', 'LineWidth', 1.5);
end

%% Set coordinate limits
xlim([-0.02, 1.02]);
ylim([-0.02, 1.02]);

%% Save figure
saveas(gcf, './figure/low_rank_intuition.png');
saveas(gcf, './figure/low_rank_intuition.eps', 'epsc');

hold off;