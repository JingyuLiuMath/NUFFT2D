%% MATLAB code: Red square positioned at upper-right, smaller size, 
% other regions partitioned accordingly with axis labels a,b,c,d
% Save as block_matrix_offset.m

clear; clc; close all;

%% Setup figure
figure('Color', 'white', 'Position', [100 100 800 800]);
hold on;
axis equal;
axis off;

%% Define colors
colors.yellow = [0.85, 0.65, 0.13];  %% Golden yellow
colors.green = [0.18, 0.55, 0.34];   %% Sea green
colors.blue = [0.27, 0.51, 0.71];    %% Steel blue
colors.red = [0.80, 0.36, 0.36];     %% Indian red

%% Define red square parameters (offset to upper-right, smaller size)
red_w = 0.15;   %% Red square width (smaller)
red_h = 0.10;   %% Red square height
red_x = 0.65;   %% Positioned to the right
red_y = 0.55;   %% Positioned upward

%% Map to a,b,c,d notation
% x-axis: [0, red_x, red_x+red_w, 1]  ->  [0, a, b, 1]
% y-axis: [0, red_y, red_y+red_h, 1]  ->  [0, c, d, 1]
a = red_x;
b = red_x + red_w;
c = red_y;
d = red_y + red_h;

%% Define partition lines based on red square position
% Vertical partition lines: x = a and x = b
% Horizontal partition lines: y = c and y = d

%% Define 9 regions
regions = struct();

% Row 1 (Upper)
regions(1).name = 'A';  % Upper-left
regions(1).x = 0; regions(1).y = d; 
regions(1).w = a; regions(1).h = 1 - d;
regions(1).color = colors.yellow;

regions(2).name = 'B';  % Upper-middle (above red square)
regions(2).x = a; regions(2).y = d;
regions(2).w = b - a; regions(2).h = 1 - d;
regions(2).color = colors.green;

regions(3).name = 'C';  % Upper-right
regions(3).x = b; regions(3).y = d;
regions(3).w = 1 - b; regions(3).h = 1 - d;
regions(3).color = colors.yellow;

% Row 2 (Middle, containing red square)
regions(4).name = 'D';  % Middle-left (left of red square)
regions(4).x = 0; regions(4).y = c;
regions(4).w = a; regions(4).h = d - c;
regions(4).color = colors.blue;

regions(5).name = 'E';  % Center - Red square
regions(5).x = a; regions(5).y = c;
regions(5).w = b - a; regions(5).h = d - c;
regions(5).color = colors.red;

regions(6).name = 'F';  % Middle-right (right of red square)
regions(6).x = b; regions(6).y = c;
regions(6).w = 1 - b; regions(6).h = d - c;
regions(6).color = colors.blue;

% Row 3 (Lower)
regions(7).name = 'G';  % Lower-left
regions(7).x = 0; regions(7).y = 0;
regions(7).w = a; regions(7).h = c;
regions(7).color = colors.yellow;

regions(8).name = 'H';  % Lower-middle (below red square)
regions(8).x = a; regions(8).y = 0;
regions(8).w = b - a; regions(8).h = c;
regions(8).color = colors.green;

regions(9).name = 'I';  % Lower-right
regions(9).x = b; regions(9).y = 0;
regions(9).w = 1 - b; regions(9).h = c;
regions(9).color = colors.yellow;

%% Draw each region (solid fill)
for i = 1:length(regions)
    r = regions(i);
    % Draw filled rectangle
    fill([r.x, r.x+r.w, r.x+r.w, r.x], [r.y, r.y, r.y+r.h, r.y+r.h], ...
         r.color, 'EdgeColor', 'k', 'LineWidth', 1.5);
end

%% Add axis labels a, b, c, d
label_fontsize = 14;
label_offset = 0.03;  % Distance from axis

% X-axis labels (below the figure)
text(a, -label_offset, '$a$', 'Interpreter', 'latex', ...
     'FontSize', label_fontsize, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top');
text(b, -label_offset, '$b$', 'Interpreter', 'latex', ...
     'FontSize', label_fontsize, 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top');

% Y-axis labels (left of the figure)
text(-label_offset, c, '$c$', 'Interpreter', 'latex', ...
     'FontSize', label_fontsize, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'middle');
text(-label_offset, d, '$d$', 'Interpreter', 'latex', ...
     'FontSize', label_fontsize, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'middle');

% Optional: Add 0 and 1 labels for completeness
% text(0, -label_offset, '$0$', 'Interpreter', 'latex', ...
%      'FontSize', label_fontsize, 'HorizontalAlignment', 'center', ...
%      'VerticalAlignment', 'top');
% text(1, -label_offset, '$1$', 'Interpreter', 'latex', ...
%      'FontSize', label_fontsize, 'HorizontalAlignment', 'center', ...
%      'VerticalAlignment', 'top');
% text(-label_offset, 0, '$0$', 'Interpreter', 'latex', ...
%      'FontSize', label_fontsize, 'HorizontalAlignment', 'right', ...
%      'VerticalAlignment', 'middle');
% text(-label_offset, 1, '$1$', 'Interpreter', 'latex', ...
%      'FontSize', label_fontsize, 'HorizontalAlignment', 'right', ...
%      'VerticalAlignment', 'middle');

%% Set coordinate limits (expanded to accommodate labels)
xlim([-0.08, 1.05]);
ylim([-0.08, 1.05]);

%% Save figure
saveas(gcf, './figure/low_rank_intuition.png');
saveas(gcf, './figure/low_rank_intuition.eps', 'epsc');

hold off;