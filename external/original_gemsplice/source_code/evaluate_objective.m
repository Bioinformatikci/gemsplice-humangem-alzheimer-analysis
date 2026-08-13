function [v1_table f] = evaluate_objective(x, M, V, fbarecon, geni, reaction_expression, pos_genes_in_react_expr)

yt=x';      % x' is the transpose of x, that is the gene expression array associated with the child (two codon usages)
yt=yt(1:V);     % needed because sometimes, and especially with child_3 in genetic_operator, there is a bug and all the child is passed to this function here, including the final rank, crowding distance and objective functions (instead of passing only the V decision variables)


eval_reaction_expression = reaction_expression;

for i=1:numel(yt) %loop over the array of the non-1 gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values. All the gene set expressions with only 1s as gene values , will be left empty and at the end of this loop everything empty will be substituted with 1 anyway. This avoids looping over all the genes yt, which is very expensive
    posizioni_gene = pos_genes_in_react_expr{i}; 
    for j=1:length(posizioni_gene)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
        eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, ['/' geni{i} '/'], num2str(yt(i),'%.15f'));  %Matlab strangely truncates decimal digits when using num2str. Addimg %.12f at least ensures that 12 decimal digits are included in the number converted into string
    end
end
eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1.0'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally



% for i=1:numel(eval_reaction_expression)
%     for j=1:numel(geni) %in the eval_reaction_expression in which already all the non-1 genes have been substituted with their values, we now have to substitute all the remaining genes with the number 1
%         eval_reaction_expression{i} = strrep(eval_reaction_expression{i}, ['/' geni{j} '/'], '1.0');  %we do this job only for the non-1 reaction expressions, so we save a lot of computational time. Note that we need '1.0' because the regecpr below is looking for NUM.NUM strings, and therefore putting only '1' will give an infinite loop
%     end
% end


%eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions or reactions whose genes have all gene expression 1) with 1, i.e. gene expressed nomally

num_reaction_expression = zeros(1,length(eval_reaction_expression));
gamma = zeros(1,length(reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    
    while (numel(strfind(str,')')) > 32) %if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM)
        substrings_to_replace = regexp(str, to_replace, 'match');
        for j = 1:numel(substrings_to_replace)
            ss_rep = substrings_to_replace{j};
            str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
        end
    end
    
    num_reaction_expression(i) = eval(str);   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
    
    %     match_geneset = find(yt == num_reaction_expression(i));     %this instruction finds the location of the resultant geneset expression in the gene expression array, so as to understand the gene responsible for this geneset expression
    %
    %     if length(match_geneset) > 1
    %         gamma(i) = 1; %in this case, there are multiple matches, which meansthe expression was exactly 1, which means almost surely that the gene in the model was not present in the probe genes, and therefore its expression and its variance must be set to 1
    %     else
    %         gamma(i) = gene_importance(match_geneset);
    %     end
end

gamma = 100000.*ones(1,length(reaction_expression));


for i=1:length(num_reaction_expression)   %loop over the array of the geneset expressions
    
    %         fbarecon.lb(i) = fbarecon.lb(i)*num_reaction_expression(i)^gamma(i);
    %         fbarecon.ub(i) = fbarecon.ub(i)*num_reaction_expression(i)^gamma(i);
    %
    if num_reaction_expression(i)>=1
        fbarecon.lb(i) = fbarecon.lb(i)*(1+gamma(i)*log(num_reaction_expression(i)));
        fbarecon.ub(i) = fbarecon.ub(i)*(1+gamma(i)*log(num_reaction_expression(i)));
    else
        fbarecon.lb(i) = fbarecon.lb(i)/(1+gamma(i)*abs(log(num_reaction_expression(i))));
        fbarecon.ub(i) = fbarecon.ub(i)/(1+gamma(i)*abs(log(num_reaction_expression(i))));
    end
end
%END NEW CODE

if isfield(fbarecon, 'h')    %then trilevel FBA
    [v, fmax, fmin, fmax_max, fmax_min, fmin_max, fmin_min, v1max_max, v1max_min, v1min_max, v1min_min] = flux_balance_trilevel(fbarecon,true);
else    %bilevel FBA
    [v, v1max, v1min, fmax, fmin] = flux_balance(fbarecon,true);
end

% objective functions number M is 2
f(1) = fbarecon.f' * v; % Biomass
f(2) = fmax; % max of the 1st synthetic obj
f(3) = fmin; % min of the 1st synthetic obj

if isfield(fbarecon, 'h')   %then trilevel FBA
    f(4) = fmax_max;
    f(5) = fmax_min;
    f(6) = fmin_max;
    f(7) = fmin_min;
end

%format longG;
format compact;
format shortG;

%f
f(1) = f(1); %put - if we maximise, + if we minimise (because nsga is by default minimising the objectives, so the minus sign means it will minimise the negative objective, which means it will maximise the objective
f(2) = f(2);
f(3) = f(3); 

if isfield(fbarecon, 'h')
    f(4) = f(4);
    f(5) = f(5);
    f(6) = f(6);
    f(7) = f(7);
end

if isfield(fbarecon, 'h')
    v1_table = [fbarecon.rxns num2cell(v), num2cell(v1min_min), num2cell(v1min_max), num2cell(v1max_min), num2cell(v1max_max)];
else
    v1_table = [fbarecon.rxns num2cell(v) num2cell(v1min) num2cell(v1max)];
end

%[concentrationMatrix,excRxnNames,timeVec,biomassVec] = my_dynamicFBA_antibiotic(fbamodel,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);