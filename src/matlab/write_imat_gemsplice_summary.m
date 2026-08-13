function write_imat_gemsplice_summary(patients, patients_results)
%WRITE_IMAT_GEMSPLICE_SUMMARY Save a compact text summary for the iMAT-GEMsplice run.

thisDir = fileparts(mfilename('fullpath'));
outdir = fullfile(thisDir, 'outputs');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

n = numel(patients_results);
sample = strings(n, 1);
group = strings(n, 1);
f = nan(n, 7);
for i = 1:n
    sample(i) = string(patients(i).SampleID);
    group(i) = string(patients(i).Group);
    vals = patients_results(i).f1;
    f(i, 1:numel(vals)) = vals;
end

T = table(sample, group, f(:,1), f(:,2), f(:,3), f(:,4), f(:,5), f(:,6), f(:,7), ...
    'VariableNames', {'sample','group','mitochondrial_complex_I','lactate_max','lactate_min', ...
    'glutamate_lactate_max_max','glutamate_lactate_max_min','glutamate_lactate_min_max','glutamate_lactate_min_min'});
writetable(T, fullfile(outdir, 'sample_objectives.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

[groups, ~, gix] = unique(group);
G = table();
for i = 1:numel(groups)
    rows = gix == i;
    G = [G; table(groups(i), sum(rows), ...
        mean(f(rows,1)), mean(f(rows,2)), mean(f(rows,5)), ...
        std(f(rows,1)), std(f(rows,2)), std(f(rows,5)), ...
        'VariableNames', {'group','n','mean_mitochondrial_complex_I','mean_lactate_max','mean_glutamate_lactate_max_min', ...
        'std_mitochondrial_complex_I','std_lactate_max','std_glutamate_lactate_max_min'})]; %#ok<AGROW>
end
writetable(G, fullfile(outdir, 'group_objective_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

disp(T);
disp(G);
end
