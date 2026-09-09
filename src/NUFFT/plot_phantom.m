function plot_phantom(n, c_ex, selected_row, figure_name)
% Plot a reconstructed phantom and a detailed view of one row.

c_ex = real(reshape(c_ex, n, n));
index = (-n / 2) : (n / 2 - 1);
save_figure = ~isempty(figure_name) && strlength(string(figure_name)) > 0;

if save_figure
    figure_name = string(figure_name);
    figure_directory = fileparts(figure_name);
    if strlength(figure_directory) > 0 && ~isfolder(figure_directory)
        mkdir(figure_directory);
    end
end

figure("Position", [100, 100, 500, 500]);
imagesc(index, index, c_ex);
axis image;
axis off;
clim([0, 1]);
colormap(gray);
if save_figure
    exportgraphics(gcf, figure_name + "_image.pdf", "ContentType", "vector");
    exportgraphics(gcf, figure_name + "_image.png", "Resolution", 300);
end

figure("Position", [100, 100, 600, 350]);
plot(index, c_ex(selected_row, :), ...
    "Color", "b", ...
    "LineWidth", 1.5);
xlim([-n / 2, n / 2]);
xticks([-n / 2, 0, n / 2]);
yticks([0, 1]);
box on;
set(gca, 'FontSize', 20);

if save_figure
    exportgraphics(gcf, figure_name + "_row.pdf", "ContentType", "vector");
    exportgraphics(gcf, figure_name + "_row.png", "Resolution", 300);
end

end
