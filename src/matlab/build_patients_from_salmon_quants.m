function build_patients_from_salmon_quants()
%BUILD_PATIENTS_FROM_SALMON_QUANTS Create patients.mat from the 3 AD + 3 Control Salmon quant.sf files.

thisDir = fileparts(mfilename('fullpath'));
workspaceRoot = fileparts(thisDir);
calismaRoot = fileparts(workspaceRoot);
metadataFile = fullfile(workspaceRoot, 'imat', 'prep_inputs', 'sample_metadata.csv');
metadata = readtable(metadataFile, 'TextType', 'string');

patients = struct('CancerType', {}, 'Group', {}, 'SampleID', {}, ...
    'GeneName', {}, 'TranscriptID', {}, 'Exp', {}, 'posizione_gene', {});

requiredTranscripts = containers.Map('KeyType','char','ValueType','any');
for groupName = ["Alzheimer", "Control"]
    [~, ~, ~, ~, groupTranscripts] = load_imat_gemsplice_model(groupName);
    groupKeys = strip_versions_local(string(groupTranscripts(:)));
    groupKeys(groupKeys == "") = [];
    requiredTranscripts(char(lower(groupName))) = unique(groupKeys);
end

for i = 1:height(metadata)
    quantFile = fullfile(calismaRoot, metadata.quant_file(i));
    if ~exist(quantFile, 'file')
        quantFile = fullfile(calismaRoot, 'data', metadata.quant_file(i));
    end
    if ~exist(quantFile, 'file')
        error('Quant file not found: %s', quantFile);
    end
    q = readtable(quantFile, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
    groupKey = char(lower(metadata.group(i)));
    keepKeys = requiredTranscripts(groupKey);
    quantKeys = strip_versions_local(q.Name);
    q = q(ismember(quantKeys, keepKeys), :);

    patients(i).CancerType = char(metadata.group(i) + "_" + metadata.sample_id(i));
    patients(i).Group = char(metadata.group(i));
    patients(i).SampleID = char(metadata.sample_id(i));
    patients(i).GeneName = cellstr(strip_versions_local(q.Name));
    patients(i).TranscriptID = cellstr(q.Name);
    patients(i).Exp = q.TPM;
    patients(i).posizione_gene = [];
    fprintf('%s: kept %d model-relevant transcripts from %s.\n', ...
        patients(i).SampleID, height(q), quantFile);
end

save(fullfile(thisDir, 'patients.mat'), 'patients', '-v7.3');
fprintf('Saved patients.mat with %d Salmon TPM samples.\n', numel(patients));
end

function ids = strip_versions_local(ids)
ids = regexprep(string(ids), '\..*$', '');
end
