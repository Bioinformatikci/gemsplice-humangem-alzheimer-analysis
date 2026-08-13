%% Plot GEMsplice-iMAT Alzheimer vs Control outputs
% Uses the same polished plot-output logic as the Human-GEM Alzheimer/Control
% tool, adapted for group-specific iMAT models.

close all;

thisDir = fileparts(mfilename('fullpath'));
resultsFile = fullfile(thisDir, 'patients_results_imat_humanGEM_salmon.mat');
if ~exist(resultsFile, 'file')
    error('patients_results_imat_humanGEM_salmon.mat not found. Run GEMsplice_RNAseq first.');
end

loaded = load(resultsFile, 'patients_results', 'patients');
patients_results = loaded.patients_results;
patients = loaded.patients;

nPatients = numel(patients_results);
mito_complex_I = zeros(1, nPatients);
pyr_kin = zeros(1, nPatients);
lac_dehy = zeros(1, nPatients);
names = cell(1, nPatients);
groups = strings(1, nPatients);

for i = 1:nPatients
    mito_complex_I(i) = patients_results(i).f1(1);
    pyr_kin(i) = patients_results(i).f1(2);
    lac_dehy(i) = patients_results(i).f1(5);
    names{i} = patients(i).CancerType;
    groups(i) = string(patients(i).Group);
end

final_table = [names' num2cell(mito_complex_I') num2cell(pyr_kin') num2cell(lac_dehy')];

index_normal = find(strcmpi(groups, 'Control'));
index_cancer = find(strcmpi(groups, 'Alzheimer'));
if isempty(index_normal)
    error('Reference group not found. Expected patients.Group == Control.');
end
if isempty(index_cancer)
    error('Comparison group not found. Expected patients.Group == Alzheimer.');
end

comparisonLabel = 'Alzheimer';
referenceLabel = 'Control';

[pathway_table, names_nonzero_reduced, x, y, r, c] = build_imat_pathway_plot_data( ...
    patients, patients_results, index_cancer, index_normal);

save_gemsplice_plot_outputs(fullfile(thisDir, 'plot_outputs'), patients, patients_results, ...
    final_table, pathway_table, names, mito_complex_I, pyr_kin, lac_dehy, ...
    index_cancer, index_normal, comparisonLabel, referenceLabel, ...
    names_nonzero_reduced, x, y, r, c);

function [pathway_table, plotNames, x, y, r, c] = build_imat_pathway_plot_data(patients, patients_results, index_cancer, index_normal)
groupNames = ["Alzheimer", "Control"];
groupData = struct();
allSubsystems = strings(0, 1);

for g = 1:numel(groupNames)
    groupName = groupNames(g);
    [fbabreast, geni, reaction_expression, pos_genes_in_react_expr, geni_transcripts] = ...
        load_imat_gemsplice_model(groupName);
    subs = normalize_subsystems(fbabreast.subSystems);
    uniqueSubs = unique(subs);
    allSubsystems = unique([allSubsystems; uniqueSubs(:)]);

    neutralExpression = ones(1, numel(geni));
    [neutralTable, ~] = evaluate_objective(neutralExpression, 3, numel(geni), ...
        fbabreast, geni, reaction_expression, pos_genes_in_react_expr);

    baseline = pathway_means_from_table(neutralTable, subs, uniqueSubs);
    sampleIxs = find(strcmpi(string({patients.Group}), groupName));
    sampleFlux = nan(numel(uniqueSubs), numel(sampleIxs));
    for j = 1:numel(sampleIxs)
        sampleFlux(:, j) = pathway_means_from_table(patients_results(sampleIxs(j)).v1_table, subs, uniqueSubs);
    end

    groupData(g).name = groupName;
    groupData(g).subsystems = uniqueSubs;
    groupData(g).baseline = baseline;
    groupData(g).sampleFlux = sampleFlux;
    groupData(g).occurrences = histcounts(categorical(subs), categorical(uniqueSubs));
end

allSubsystems = sort(allSubsystems);
alzDerivative = nan(numel(allSubsystems), 1);
ctrlDerivative = nan(numel(allSubsystems), 1);
occurrences = zeros(numel(allSubsystems), 1);
errColor = nan(numel(allSubsystems), 1);

for g = 1:numel(groupData)
    [tf, loc] = ismember(groupData(g).subsystems, allSubsystems);
    loc = loc(tf);
    avgFlux = mean(groupData(g).sampleFlux(tf, :), 2, 'omitnan');
    baseline = groupData(g).baseline(tf);
    derivative = nan(size(avgFlux));
    valid = isfinite(avgFlux) & avgFlux ~= 0 & isfinite(baseline);
    derivative(valid) = (avgFlux(valid) - baseline(valid)) ./ avgFlux(valid);

    normalizedFlux = groupData(g).sampleFlux(tf, :) ./ baseline;
    stderr = std(normalizedFlux, 0, 2, 'omitnan') ./ sqrt(max(groupData(g).occurrences(tf)', 1));
    groupColor = log(stderr);

    if strcmpi(groupData(g).name, 'Alzheimer')
        alzDerivative(loc) = derivative;
    else
        ctrlDerivative(loc) = derivative;
    end
    occurrences(loc) = max(occurrences(loc), groupData(g).occurrences(tf)');
    errColor(loc) = max_or_fill(errColor(loc), groupColor);
end

pathway_table = [cellstr(allSubsystems(:)) num2cell(ctrlDerivative(:)) num2cell(alzDerivative(:))];

excludeNames = ["", "exchange", "Transport, Nuclear", "Miscellaneous"];
finiteMask = isfinite(ctrlDerivative) & isfinite(alzDerivative) & ...
    ~ismember(lower(allSubsystems), lower(excludeNames));
plotNames = cellstr(allSubsystems(finiteMask));
x = ctrlDerivative(finiteMask);
y = alzDerivative(finiteMask);
r = occurrences(finiteMask);
c = errColor(finiteMask);
c(~isfinite(c)) = 0;
end

function means = pathway_means_from_table(v1Table, subs, uniqueSubs)
fluxColumn = cell2mat(v1Table(:, 5));
means = nan(numel(uniqueSubs), 1);
for i = 1:numel(uniqueSubs)
    means(i) = mean(fluxColumn(subs == uniqueSubs(i)), 'omitnan');
end
end

function subs = normalize_subsystems(subSystems)
subs = strings(numel(subSystems), 1);
for i = 1:numel(subSystems)
    value = subSystems{i};
    if iscell(value)
        value = value{1};
    end
    if isstring(value)
        value = char(value);
    end
    if isempty(value)
        value = "";
    end
    subs(i) = string(value);
end
end

function out = max_or_fill(existing, candidate)
out = existing;
for i = 1:numel(candidate)
    if ~isfinite(out(i))
        out(i) = candidate(i);
    elseif isfinite(candidate(i))
        out(i) = max(out(i), candidate(i));
    end
end
end
