function [dati_micro, posizione_gene] = build_imat_dati_micro(patient, geni, geni_transcripts)
%BUILD_IMAT_DATI_MICRO Build transcript expression vector for one Salmon sample.

probe_genes_patient = regexprep(patient.TranscriptID, '\..*', '');
exprByTranscript = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(probe_genes_patient)
    key = probe_genes_patient{i};
    if isempty(key)
        continue;
    end
    exprValue = log2(max(patient.Exp(i), 0) + 1);
    if isKey(exprByTranscript, key)
        exprByTranscript(key) = mean([exprByTranscript(key), exprValue]);
    else
        exprByTranscript(key) = exprValue;
    end
end

posizione_gene = cell(length(geni),1);
dati_micro = ones(length(geni),1);
for i = 1:length(geni)
    current_gene = geni_transcripts{i};
    if iscell(current_gene)
        query = regexprep(current_gene, '\..*', '');
    else
        query = {regexprep(current_gene, '\..*', '')};
    end
    query(cellfun(@isempty, query)) = [];
    for j = 1:numel(query)
        if isKey(exprByTranscript, query{j})
            dati_micro(i) = exprByTranscript(query{j});
            break;
        end
    end
end
return;

if isfield(patient,'posizione_gene') && ~isempty(patient.posizione_gene) ...
        && numel(patient.posizione_gene) == numel(geni)
    posizione_gene = patient.posizione_gene;
else
    probe_genes_patient = regexprep(patient.TranscriptID, '\..*', '');
    patient_index = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(probe_genes_patient)
        key = probe_genes_patient{i};
        if isempty(key)
            continue;
        end
        if isKey(patient_index, key)
            patient_index(key) = [patient_index(key); i];
        else
            patient_index(key) = i;
        end
    end

    posizione_gene = cell(length(geni),1);
    for i = 1:length(geni)
        current_gene = geni_transcripts{i};
        if iscell(current_gene)
            query = regexprep(current_gene, '\..*', '');
        else
            query = {regexprep(current_gene, '\..*', '')};
        end
        query(cellfun(@isempty, query)) = [];
        matches = [];
        for j = 1:numel(query)
            if isKey(patient_index, query{j})
                matches = patient_index(query{j});
                break;
            end
        end
        posizione_gene{i} = matches;
    end
end

microarray_data = log2(max(patient.Exp, 0) + 1);
dati_micro = zeros(length(geni),1);
for i = 1:length(geni)
    if isempty(posizione_gene{i})
        dati_micro(i) = 1;
    else
        dati_micro(i) = mean(microarray_data(posizione_gene{i}));
    end
end
end
