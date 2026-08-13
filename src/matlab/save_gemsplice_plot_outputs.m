function save_gemsplice_plot_outputs(outputDir, patients, patients_results, final_table, pathway_table, names, biomass, pyr_kin, lac_dehy, index_comparison, index_reference, comparisonLabel, referenceLabel, pathwayPlotNames, pathwayX, pathwayY, pathwaySize, pathwayColor)
%SAVE_GEMSPLICE_PLOT_OUTPUTS Save polished GEMsplice figures and English notes.

if nargin < 1 || isempty(outputDir)
    outputDir = 'plot_outputs';
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

oldFiles = [dir(fullfile(outputDir,'plot*.png')); ...
            dir(fullfile(outputDir,'plot*.pdf')); ...
            dir(fullfile(outputDir,'plot*.fig')); ...
            dir(fullfile(outputDir,'plot*text.txt')); ...
            dir(fullfile(outputDir,'summary_table.csv')); ...
            dir(fullfile(outputDir,'pathway_table.csv'))];
for i = 1:numel(oldFiles)
    delete(fullfile(oldFiles(i).folder, oldFiles(i).name));
end

oldFigs = findall(0, 'Type', 'figure');
if ~isempty(oldFigs)
    close(oldFigs);
end

if nargin < 14
    pathwayPlotNames = {};
    pathwayX = [];
    pathwayY = [];
    pathwaySize = [];
    pathwayColor = [];
end

summary = table(names(:), biomass(:), pyr_kin(:), lac_dehy(:), ...
    'VariableNames', {'Group','MitochondrialComplexI','LactateDehydrogenase','GlutamateDehydrogenase'});
writetable(summary, fullfile(outputDir, 'summary_table.csv'));

if exist('pathway_table','var') && ~isempty(pathway_table)
    pathwaySummary = cell2table(pathway_table, ...
        'VariableNames', {'Pathway','ReferenceRelativeFluxChange','ComparisonRelativeFluxChange'});
    writetable(pathwaySummary, fullfile(outputDir, 'pathway_table.csv'));
end

figs = create_polished_figures(summary, pathwayPlotNames, pathwayX, pathwayY, pathwaySize, pathwayColor, index_reference, comparisonLabel, referenceLabel);

for i = 1:numel(figs)
    fig = figs(i);
    base = fullfile(outputDir, sprintf('plot%d', i));
    savefig(fig, [base '.fig']);
    exportgraphics(fig, [base '.png'], 'Resolution', 300);
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    write_plot_text([base 'text.txt'], i, patients, patients_results, summary, pathway_table, index_comparison, index_reference, comparisonLabel, referenceLabel, pathwayPlotNames, pathwayX, pathwayY, pathwaySize, pathwayColor);
end

end

function figs = create_polished_figures(summary, pathwayNames, pathwayX, pathwayY, pathwaySize, pathwayColor, index_reference, comparisonLabel, referenceLabel)
comparisonColor = [0.00 0.45 0.48];
referenceColor = [0.90 0.10 0.55];
figs = gobjects(1,3);

figs(1) = figure('Color','w','Position',[100 100 1050 760]);
ax = axes(figs(1));
hold(ax,'on');
legendHandles = gobjects(1,2);
for i = 1:height(summary)
    [color, marker, groupIx] = point_style(i, index_reference, referenceColor, comparisonColor);
    h = scatter3(ax, summary.MitochondrialComplexI(i), summary.LactateDehydrogenase(i), summary.GlutamateDehydrogenase(i), ...
        120, marker, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.9);
    if ~isgraphics(legendHandles(groupIx))
        legendHandles(groupIx) = h;
    end
    text(ax, summary.MitochondrialComplexI(i), summary.LactateDehydrogenase(i), summary.GlutamateDehydrogenase(i), ...
        ['  ' summary.Group{i}], 'FontSize', 10, 'FontName', 'Arial', 'Color', [0.05 0.05 0.05]);
end
xlabel(ax,'Mitochondrial Complex I objective');
ylabel(ax,'Lactate dehydrogenase objective');
zlabel(ax,'Glutamate dehydrogenase objective');
title(ax,'GEMsplice Objective Space: Alzheimer vs Control', 'Color', [0.05 0.05 0.05]);
legend(ax, legendHandles(isgraphics(legendHandles)), {comparisonLabel, referenceLabel}, 'Location','best');
view(ax, 38, 24);
grid(ax,'on');
box(ax,'on');
style_axes(ax);

figs(2) = figure('Color','w','Position',[100 100 1250 430]);
tiledlayout(figs(2),1,3,'Padding','compact','TileSpacing','compact');
plot_pair(nexttile, summary.MitochondrialComplexI, summary.LactateDehydrogenase, summary, index_reference, comparisonColor, referenceColor, 'Mitochondrial Complex I objective', 'Lactate dehydrogenase objective', comparisonLabel, referenceLabel, -1, 1);
plot_pair(nexttile, summary.MitochondrialComplexI, summary.GlutamateDehydrogenase, summary, index_reference, comparisonColor, referenceColor, 'Mitochondrial Complex I objective', 'Glutamate dehydrogenase objective', comparisonLabel, referenceLabel, -1, 1);
plot_pair(nexttile, summary.LactateDehydrogenase, summary.GlutamateDehydrogenase, summary, index_reference, comparisonColor, referenceColor, 'Lactate dehydrogenase objective', 'Glutamate dehydrogenase objective', comparisonLabel, referenceLabel, 1, 1);
sgtitle('Pairwise Objective Projections: Alzheimer vs Control', 'FontName','Arial', 'FontSize', 14, 'FontWeight','bold', 'Color', [0.05 0.05 0.05]);

figs(3) = figure('Color','w','Position',[100 100 1050 760]);
ax = axes(figs(3));
if ~isempty(pathwayX)
    finiteMask = isfinite(pathwayX(:)) & isfinite(pathwayY(:)) & isfinite(pathwayColor(:));
    pxRaw = pathwayX(finiteMask);
    pyRaw = pathwayY(finiteMask);
    px = signed_log10(pxRaw);
    py = signed_log10(pyRaw);
    pc = pathwayColor(finiteMask);
    ps = pathwaySize(finiteMask);
    pn = pathwayNames(finiteMask);
    if isempty(ps)
        scaledSize = 45 * ones(size(px));
    else
        ps = double(ps(:));
        if max(ps) == min(ps)
            scaledSize = 70 * ones(size(ps));
        else
            scaledSize = 35 + 165 * (ps - min(ps)) ./ (max(ps) - min(ps));
        end
    end
    scatter(ax, px, py, scaledSize, pc, 'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.4);
    colormap(ax, parula(256));
    cb = colorbar(ax);
    ylabel(cb, 'log normalized error');
    hold(ax,'on');
    refline(ax, 1, 0);
    [~, order] = sort(abs(py - px), 'descend');
    for i = 1:min(8,numel(order))
        ii = order(i);
        text(ax, px(ii), py(ii), ['  ' char(pn{ii})], 'FontSize', 8, 'FontName', 'Arial', 'Color', [0.05 0.05 0.05]);
    end
else
    text(ax, 0.5, 0.5, 'No finite pathway points available for plotting.', 'HorizontalAlignment','center');
end
xlabel(ax, ['signed log10(1 + abs(' referenceLabel ' relative pathway flux change))']);
ylabel(ax, ['signed log10(1 + abs(' comparisonLabel ' relative pathway flux change))']);
title(ax,'Pathway-Level Relative Flux Change', 'Color', [0.05 0.05 0.05]);
grid(ax,'on');
box(ax,'on');
style_axes(ax);

end

function [color, marker, groupIx] = point_style(i, index_reference, referenceColor, comparisonColor)
if ismember(i, index_reference)
    color = referenceColor;
    marker = 'o';
    groupIx = 2;
else
    color = comparisonColor;
    marker = '^';
    groupIx = 1;
end
end

function plot_pair(ax, x, y, summary, index_reference, comparisonColor, referenceColor, xlab, ylab, comparisonLabel, referenceLabel, xParetoSign, yParetoSign)
hold(ax,'on');
plot_pareto_fronts(ax, x, y, summary, index_reference, comparisonColor, referenceColor, comparisonLabel, referenceLabel, xParetoSign, yParetoSign);
for i = 1:height(summary)
    [color, marker] = point_style(i, index_reference, referenceColor, comparisonColor);
    scatter(ax, x(i), y(i), 90, marker, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    text(ax, x(i), y(i), ['  ' summary.Group{i}], 'FontSize', 9, 'FontName', 'Arial');
end
xlabel(ax, xlab);
ylabel(ax, ylab);
grid(ax,'on');
box(ax,'on');
pad_axis(ax, x, y);
style_axes(ax);
end

function plot_pareto_fronts(ax, x, y, summary, index_reference, comparisonColor, referenceColor, comparisonLabel, referenceLabel, xParetoSign, yParetoSign)
comparisonIx = setdiff(1:height(summary), index_reference);
lineHandles = gobjects(1, 2);
lineHandles(1) = plot_front(ax, x(comparisonIx), y(comparisonIx), comparisonColor, comparisonLabel, xParetoSign, yParetoSign);
lineHandles(2) = plot_front(ax, x(index_reference), y(index_reference), referenceColor, referenceLabel, xParetoSign, yParetoSign);
valid = isgraphics(lineHandles);
labels = {comparisonLabel, referenceLabel};
if any(valid)
    legend(ax, lineHandles(valid), labels(valid), 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
end
end

function h = plot_front(ax, x, y, color, label, xParetoSign, yParetoSign)
h = gobjects(1);
x = double(x(:));
y = double(y(:));
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if isempty(x)
    return;
end

if exist('non_domination_sort_mod', 'file') == 2
    front = non_domination_sort_mod([xParetoSign .* x, yParetoSign .* y], 2, 0);
    front = front(front(:,3) == 1, 1:2);
    front(:,1) = front(:,1) ./ xParetoSign;
    front(:,2) = front(:,2) ./ yParetoSign;
else
    front = [x y];
end
front = sortrows(front, 1);
h = plot(ax, front(:,1), front(:,2), '-', 'Color', color, 'LineWidth', 1.5, 'DisplayName', label);
end

function pad_axis(ax, x, y)
x = x(:);
y = y(:);
dx = max(x) - min(x);
dy = max(y) - min(y);
if dx == 0
    dx = max(abs(x(1))*0.05, 1);
end
if dy == 0
    dy = max(abs(y(1))*0.05, 1);
end
xlim(ax, [min(x)-0.15*dx, max(x)+0.30*dx]);
ylim(ax, [min(y)-0.15*dy, max(y)+0.30*dy]);
end

function y = signed_log10(x)
y = sign(x) .* log10(1 + abs(x));
end

function style_axes(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1, ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], ...
    'ZColor', [0.15 0.15 0.15], 'Color', 'w');
end

function write_plot_text(filePath, plotNumber, patients, patients_results, summary, pathway_table, index_comparison, index_reference, comparisonLabel, referenceLabel, pathwayPlotNames, pathwayX, pathwayY, pathwaySize, pathwayColor)
fid = fopen(filePath, 'w');
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, 'GEMsplice-iMAT Human-GEM Alzheimer vs Control - Plot %d\n', plotNumber);
fprintf(fid, 'Generated by the plotting script after the GEMsplice/FBA results were computed.\n\n');

if plotNumber == 1
    fprintf(fid, 'Plot type: 3D objective scatter.\n');
    fprintf(fid, 'Each point represents one patient/group profile.\n');
    fprintf(fid, 'X axis: Mitochondrial Complex I objective.\n');
    fprintf(fid, 'Y axis: Lactate dehydrogenase objective.\n');
    fprintf(fid, 'Z axis: Glutamate dehydrogenase objective.\n\n');
elseif plotNumber == 2
    fprintf(fid, 'Plot type: Multi-panel pairwise objective projection view.\n');
    fprintf(fid, 'Each point represents one patient/group profile.\n');
    fprintf(fid, 'The panels show the same objective values from different 2D projections.\n');
    fprintf(fid, 'Pareto/front lines are drawn with the same non-domination sign convention used in the original GEMsplice plot2.\n\n');
elseif plotNumber == 3
    fprintf(fid, 'Plot type: Pathway-level relative flux-change scatter.\n');
    fprintf(fid, 'Each point represents one metabolic pathway.\n');
    fprintf(fid, 'X axis: signed log10(1 + abs(%s relative pathway flux change)); sign is preserved.\n', referenceLabel);
    fprintf(fid, 'Y axis: signed log10(1 + abs(%s relative pathway flux change)); sign is preserved.\n', comparisonLabel);
    fprintf(fid, 'Point size: scaled number of reactions assigned to the pathway.\n');
    fprintf(fid, 'Point color: log normalized error as computed by the original plotting logic.\n\n');
else
    fprintf(fid, 'Plot type: Additional figure generated by the original plotting script.\n\n');
end

fprintf(fid, 'Group mapping:\n');
fprintf(fid, 'Reference group: %s (indices: %s)\n', referenceLabel, mat2str(index_reference));
fprintf(fid, 'Comparison group: %s (indices: %s)\n\n', comparisonLabel, mat2str(index_comparison));

fprintf(fid, 'Objective values by plotted point:\n');
for i = 1:height(summary)
    fprintf(fid, '%d. %s | Mitochondrial Complex I = %.10g | Lactate dehydrogenase = %.10g | Glutamate dehydrogenase = %.10g\n', ...
        i, summary.Group{i}, summary.MitochondrialComplexI(i), summary.LactateDehydrogenase(i), summary.GlutamateDehydrogenase(i));
end

fprintf(fid, '\nPatient result vector f1 values:\n');
for i = 1:numel(patients_results)
    fprintf(fid, '%s: ', patients(i).CancerType);
    fprintf(fid, '%.10g ', patients_results(i).f1);
    fprintf(fid, '\n');
end

if plotNumber == 3 && exist('pathwayPlotNames','var') && ~isempty(pathwayPlotNames)
    finiteMask = isfinite(pathwayX(:)) & isfinite(pathwayY(:)) & isfinite(pathwayColor(:));
    finiteNames = pathwayPlotNames(finiteMask);
    finiteX = pathwayX(finiteMask);
    finiteY = pathwayY(finiteMask);
    finiteSize = pathwaySize(finiteMask);
    finiteColor = pathwayColor(finiteMask);
    fprintf(fid, '\nPlotted pathway point mapping:\n');
    fprintf(fid, 'Pathway | Reference relative flux change | Comparison relative flux change | Reaction count | Color value\n');
    for i = 1:numel(finiteNames)
        fprintf(fid, '%s | %.10g | %.10g | %.10g | %.10g\n', finiteNames{i}, finiteX(i), finiteY(i), finiteSize(i), finiteColor(i));
    end
    fprintf(fid, '\nNon-finite pathway values are omitted from the plotted scatter but preserved in pathway_table.csv.\n');
elseif plotNumber == 3 && exist('pathway_table','var') && ~isempty(pathway_table)
    fprintf(fid, '\nFull pathway table is saved as pathway_table.csv.\n');
end
end
