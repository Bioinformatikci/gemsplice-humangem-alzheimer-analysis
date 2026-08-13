function [fbabreast, geni, reaction_expression, pos_genes_in_react_expr, geni_transcripts] = load_imat_gemsplice_model(groupName)
%LOAD_IMAT_GEMSPLICE_MODEL Prepare the common iMAT-reduced Human-GEM model for GEMsplice.

if nargin < 1 || isempty(groupName)
    groupName = 'Common';
end

groupKey = 'common';

cacheFile = sprintf('imat_gemsplice_model_%s.mat', groupKey);
cacheVersion = 'imat_humangem_common_tx2gene_AD_objectives_v1';
if exist(cacheFile, 'file')
    cache = load(cacheFile);
    if isfield(cache, 'cacheVersion') && strcmp(cache.cacheVersion, cacheVersion)
        fbabreast = cache.fbabreast;
        geni = cache.geni;
        reaction_expression = cache.reaction_expression;
        pos_genes_in_react_expr = cache.pos_genes_in_react_expr;
        geni_transcripts = cache.geni_transcripts;
        return;
    end
end

thisDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(thisDir);
imatDir = fullfile(workspaceDir, 'imat');
modelFile = fullfile(imatDir, 'outputs', 'imat_common_model_deadend_clean.mat');
if ~exist(modelFile, 'file')
    modelFile = fullfile(imatDir, 'outputs', 'imat_common_model.mat');
end
if ~exist(modelFile, 'file')
    error('iMAT model not found: %s', modelFile);
end

loaded = load(modelFile);
fbabreast = loaded.tissueModel;

nRxns = numel(fbabreast.rxns);
fbabreast.f = zeros(nRxns, 1);
fbabreast.g = zeros(nRxns, 1);
fbabreast.h = zeros(nRxns, 1);

ixMitoComplexI = find_imat_reaction(fbabreast, 'MAR06921', 'NADH:ubiquinone oxidoreductase');
ixLactateDehydrogenase = find_imat_reaction(fbabreast, 'MAR04388', '(S)-Lactate:NAD+ oxidoreductase');
ixGlutamateDehydrogenase = find_imat_reaction(fbabreast, 'MAR03802', 'L-glutamate:NAD+ oxidoreductase (deaminating)');
fbabreast.f(ixMitoComplexI) = 1;
fbabreast.g(ixLactateDehydrogenase) = 1;
fbabreast.h(ixGlutamateDehydrogenase) = 1;

tx2gene = readtable(fullfile(imatDir, 'prep_inputs', 'tx2gene.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
tx2gene.gene_id = strip_versions(string(tx2gene.gene_id));
tx2gene.transcript_id = strip_versions(string(tx2gene.transcript_id));

modelGenes = strip_versions(string(fbabreast.genes(:)));
gene_to_isoform_expr = containers.Map('KeyType','char','ValueType','char');
gene_to_isoform_ixs = cell(numel(modelGenes), 1);
geni = {};
geni_transcripts = {};

for i = 1:numel(modelGenes)
    tx = unique(tx2gene.transcript_id(tx2gene.gene_id == modelGenes(i)));
    tx(tx == "") = [];
    if isempty(tx)
        tx = modelGenes(i);
    end

    gene_to_isoform_ixs{i} = zeros(numel(tx), 1);
    txCell = cellstr(tx(:));
    for j = 1:numel(txCell)
        geni{end+1, 1} = txCell{j}; %#ok<AGROW>
        gene_to_isoform_ixs{i}(j) = numel(geni);
        geni_transcripts{end+1, 1} = txCell(j); %#ok<AGROW>
    end
    gene_to_isoform_expr(char(modelGenes(i))) = join_with_function(txCell, 'max');
end

reaction_expression = cell(nRxns, 1);
for i = 1:nRxns
    if isfield(fbabreast, 'grRules') && ~isempty(fbabreast.grRules{i})
        reaction_expression{i} = parse_gpr_to_expression(fbabreast.grRules{i}, gene_to_isoform_expr);
    else
        reaction_expression{i} = '';
    end
end

pos_genes_in_react_expr = cell(numel(geni), 1);
if isfield(fbabreast, 'rxnGeneMat') && ~isempty(fbabreast.rxnGeneMat)
    for i = 1:numel(modelGenes)
        rxnIxs = find(fbabreast.rxnGeneMat(:, i))';
        isoformIxs = gene_to_isoform_ixs{i};
        for j = 1:numel(isoformIxs)
            pos_genes_in_react_expr{isoformIxs(j)} = rxnIxs;
        end
    end
else
    for i = 1:numel(geni)
        pos_genes_in_react_expr{i} = find(contains(reaction_expression, ['/' geni{i} '/']))';
    end
end

save(cacheFile, 'cacheVersion', 'fbabreast', 'geni', 'reaction_expression', ...
    'pos_genes_in_react_expr', 'geni_transcripts', '-v7.3');
end

function ids = strip_versions(ids)
ids = regexprep(ids, '\..*$', '');
end

function ix = find_imat_reaction(model, rxnID, rxnName)
ix = find(strcmp(model.rxns, rxnID), 1);
if isempty(ix) && isfield(model, 'rxnNames')
    ix = find(strcmp(model.rxnNames, rxnName), 1);
end
if isempty(ix)
    error('Required objective reaction not found in iMAT model: %s / %s', rxnID, rxnName);
end
end

function expression = parse_gpr_to_expression(rule, gene_to_isoform_expr)
tokens = regexp(rule, '\(|\)|\band\b|\bor\b|[A-Za-z0-9_.:-]+', 'match');
pos = 1;
[expression, pos] = parse_or(tokens, pos, gene_to_isoform_expr);
if pos <= numel(tokens)
    error('Could not parse grRule: %s', rule);
end
end

function [node, pos] = parse_or(tokens, pos, gene_to_isoform_expr)
[node, pos] = parse_and(tokens, pos, gene_to_isoform_expr);
while pos <= numel(tokens) && strcmpi(tokens{pos}, 'or')
    pos = pos + 1;
    [rhs, pos] = parse_and(tokens, pos, gene_to_isoform_expr);
    node = ['max(' node ',' rhs ')'];
end
end

function [node, pos] = parse_and(tokens, pos, gene_to_isoform_expr)
[node, pos] = parse_factor(tokens, pos, gene_to_isoform_expr);
while pos <= numel(tokens) && strcmpi(tokens{pos}, 'and')
    pos = pos + 1;
    [rhs, pos] = parse_factor(tokens, pos, gene_to_isoform_expr);
    node = ['min(' node ',' rhs ')'];
end
end

function [node, pos] = parse_factor(tokens, pos, gene_to_isoform_expr)
if strcmp(tokens{pos}, '(')
    pos = pos + 1;
    [node, pos] = parse_or(tokens, pos, gene_to_isoform_expr);
    if pos > numel(tokens) || ~strcmp(tokens{pos}, ')')
        error('Unbalanced parentheses in grRule');
    end
    pos = pos + 1;
else
    token = regexprep(tokens{pos}, '\..*$', '');
    if isKey(gene_to_isoform_expr, token)
        node = gene_to_isoform_expr(token);
    else
        node = ['/' token '/'];
    end
    pos = pos + 1;
end
end

function expression = join_with_function(ids, fnName)
if numel(ids) == 1
    expression = ['/' ids{1} '/'];
else
    expression = ['/' ids{1} '/'];
    for i = 2:numel(ids)
        expression = [fnName '(' expression ',/' ids{i} '/'];
        expression = [expression ')'];
    end
end
end
