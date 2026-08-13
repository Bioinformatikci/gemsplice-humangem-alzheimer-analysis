format long
if ~exist('patients','var'), load('patients.mat'); end
if ~exist('geni','var'), load('geni.mat'); end
if ~exist('reaction_expression','var'), load('reaction_expression.mat'); end
if ~exist('geni_transcripts','var'), load('geni_transcripts.mat'); end
if ~exist('pos_genes_in_react_expr','var'), load('pos_genes_in_react_expr.mat'); end
if ~exist('fbabreast_trilevel','var'), load('fbabreast_trilevel.mat'); end
addAttachedFiles(gcp,{'glpk.m','glpkcc'});  %adds these files to the independent workers of the parallel pool

ixs_ENS_transcripts = find(strncmpi(geni_transcripts,'ENS',3));     %these are the 130 genes where we have annotation in terms of transcript information in Recon

%% PERFORM METABOLIC CONTROL ANALYSIS (PUT 1 HERE to RE-RUN the MCA, LEAVE 0 if THE RESULTS HAVE BEEN ALREADY OBTAINED AND LOADED)
%%
if 1
    delta = 0.001;
    
    results_mca = cell(numel(patients),1);
    for j = 1:numel(patients)
        j
        results_plus = cell(numel(ixs_ENS_transcripts),1);
        results_minus = cell(numel(ixs_ENS_transcripts),1);
        mca_plus = cell(numel(ixs_ENS_transcripts),1);
        mca_minus = cell(numel(ixs_ENS_transcripts),1);
        
        
        patient = patients(j);
        posizione_gene = patient.posizione_gene;
        
        
        %% without perturbation (these can be imported from "patients_results.mat")
        microarray_data = patient.Exp;
        dati_micro = zeros(length(geni),1);
        parfor i=1:length(geni)
            if isempty(posizione_gene{i})
                dati_micro(i) = 1; %if data is available, put here the average gene expression of that gene (across condition)
            else
                dati_micro(i) = mean(microarray_data(posizione_gene{i})); %if the gene is mapped to more probe locations (as it happens with the isoforms for which we don't have annotation and we had to use the gene name), we average the gene expression level in these probe locations
            end
        end
        
        [out1, out2] = evaluate_objective(dati_micro,3,numel(geni),fbabreast,geni,reaction_expression, pos_genes_in_react_expr);
        out2
        
        
        
        %% with positive and negative perturbation
        parfor ii = 1:numel(ixs_ENS_transcripts)
            ii
            ix_current_transcript = ixs_ENS_transcripts(ii);
            
            dati_plus = dati_micro;
            dati_plus(ix_current_transcript) = dati_micro(ix_current_transcript) + delta*dati_micro(ix_current_transcript);
            [out1_plus, out2_plus] = evaluate_objective(dati_plus,3,numel(geni),fbabreast,geni,reaction_expression, pos_genes_in_react_expr);
            results_plus{ii} = out2_plus;
            
            dati_minus = dati_micro;
            dati_minus(ix_current_transcript) = dati_micro(ix_current_transcript) - delta*dati_micro(ix_current_transcript);
            [out1_minus, out2_minus] = evaluate_objective(dati_minus,3,numel(geni),fbabreast,geni,reaction_expression, pos_genes_in_react_expr);
            results_minus{ii} = out2_minus;
            
            aux_mca_plus = ((out2_plus - out2)./out2) ./ delta;    %the formula of metabolic flux analysis is {[f(x+dx) - f(x)] / f(x)} / {(x+dx - x) / x}
            aux_mca_plus(find(out2_plus - out2 == 0)) = 0;   %those whse numerator is zero, it means the mca is zero, no matter what the denominator out2 is equal to. This "overwrites" the fact that, if out2 is zro, the full ratio becomes NaN, while we can certainly say that it should be 0 because the numerator is zero and we are measuring the perturbation out2_plus - out2
            mca_plus{ii} = aux_mca_plus;
            
            aux_mca_minus = ((out2_minus - out2)./out2) ./ delta;
            aux_mca_minus(find(out2_minus - out2 == 0)) = 0;
            mca_minus{ii} = aux_mca_minus;
            
            [mca_plus{ii}; mca_minus{ii}]
        end
        results_mca{j}.results_minus = results_minus;
        results_mca{j}.results_plus = results_plus;
        results_mca{j}.mca_plus = mca_plus;
        results_mca{j}.mca_minus = mca_minus;
        
        clc;
        
    end
end



%% PLOT THE RESULTS OF THE METABOLIC CONTROL ANALYSIS
%%
if ~exist('results_mca_eps0.001.mat','var'), load('results_mca_eps0.001.mat'); end

for j = 1:numel(results_mca) %patients/cancers
    positive_mca = results_mca{j}.mca_plus;
    negative_mca = results_mca{j}.mca_minus;
    for i = 1:numel(positive_mca) %transcripts
        tmp1 = positive_mca{i};
        tmp2 = negative_mca{i};
        mca_biomass(i,j) = max( abs(tmp1(1)),abs(tmp2(1)) );  %we are interested in 1st, 2nd and 5th outputs in this case (see notes_breast_cancer.txt)
        mca_pyr_kinase(i,j) = max( abs(tmp1(2)),abs(tmp2(2)) );
        mca_lact_dehydrogenase(i,j) = max( abs(tmp1(5)),abs(tmp2(5)) );
    end
end

names_transcripts = geni_transcripts(ixs_ENS_transcripts);

% we now remove NaN columns
mca_biomass(:,any(isnan(mca_biomass)))=[];
mca_pyr_kinase(:,any(isnan(mca_pyr_kinase)))=[];
mca_lact_dehydrogenase(:,any(isnan(mca_lact_dehydrogenase)))=[];

% % Remove zero rows
% ixs_nonzero_biomass = find(~all(~mca_biomass,2));
% mca_biomass = mca_biomass(ixs_nonzero_biomass,:);
% names_transcripts_nonzero_biomass = names_transcripts(ixs_nonzero_biomass);
% 
% ixs_nonzero_pyr_kinase = find(~all(~mca_pyr_kinase,2));
% mca_pyr_kinase = mca_pyr_kinase(ixs_nonzero_pyr_kinase,:);
% names_transcripts_nonzero_pyr_kinase = names_transcripts(ixs_nonzero_pyr_kinase);
% 
% ixs_nonzero_lact_dehydrogenase = find(~all(~mca_lact_dehydrogenase,2));
% mca_lact_dehydrogenase = mca_lact_dehydrogenase(ixs_nonzero_lact_dehydrogenase,:);
% names_transcripts_nonzero_lact_dehydrogenase = names_transcripts(ixs_nonzero_lact_dehydrogenase);

% Compute average and standard deviation across 31 cancer cells, for each transcript
avg_mca_transcript_biomass = mean(mca_biomass');
std_mca_transcript_biomass = std(mca_biomass');
avg_mca_transcript_pyr_kinase = mean(mca_pyr_kinase');
std_mca_transcript_pyr_kinase = std(mca_pyr_kinase');
avg_mca_transcript_lact_dehydrogenase = mean(mca_lact_dehydrogenase');
std_mca_transcript_lact_dehydrogenase = std(mca_lact_dehydrogenase');

[sorted_avg_mca_transcript_biomass, ixs_sorted_biomass] = sort(avg_mca_transcript_biomass,'descend');
[sorted_avg_mca_transcript_pyr_kinase, ixs_sorted_pyr_kinase] = sort(avg_mca_transcript_pyr_kinase,'descend');
[sorted_avg_mca_transcript_lact_dehydrogenase, ixs_sorted_lact_dehydrogenase] = sort(avg_mca_transcript_lact_dehydrogenase,'descend');

ixs_union_top10 = unique([ixs_sorted_biomass(1:10), ixs_sorted_pyr_kinase(1:10), ixs_sorted_lact_dehydrogenase(1:10)]); % concatenation+unique computes the union of the top 10 taken from each set

%plot multiple bar plot
my_barvalues = [avg_mca_transcript_biomass(ixs_union_top10)' avg_mca_transcript_pyr_kinase(ixs_union_top10)' avg_mca_transcript_lact_dehydrogenase(ixs_union_top10)'];
my_errors = [std_mca_transcript_biomass(ixs_union_top10)'./sqrt(size(mca_biomass,2)) std_mca_transcript_pyr_kinase(ixs_union_top10)'./sqrt(size(mca_pyr_kinase,2)) std_mca_transcript_lact_dehydrogenase(ixs_union_top10)'./sqrt(size(mca_lact_dehydrogenase,2))]; %we give the standard error of the mean, namely std/sqrt(n)
figure;
handles = barweb(my_barvalues, my_errors, 1, names_transcripts(ixs_union_top10), [], [], [], cmap('teal'), [], {'Biomass control','Pyruvate kinase control','Lactate dehydrogenase control'}, 2, 'plot');
handles.ax.XTickLabelRotation=-35;
handles.legend.Location='northeast';
handles.legend.FontSize=handles.legend.FontSize - 2;
set(gca,'ygrid','on');
ylim([0 450]);

%exportfig(gcf, 'inset.pdf', 'color', 'cmyk', 'Width', '15', 'Height', '12', 'FontMode', 'scaled', 'FontSize', '0.8' );

%save('export_to_R.mat','avg_mca_transcript_biomass','std_mca_transcript_biomass','avg_mca_transcript_pyr_kinase','std_mca_transcript_pyr_kinase','avg_mca_transcript_lact_dehydrogenase','std_mca_transcript_lact_dehydrogenase','names_transcripts','ixs_union_top10');



% figure;
% h = barwitherr(std_mca_transcript_biomass, avg_mca_transcript_biomass);
% set(gca,'yscale','log');
% limsy = get(gca,'YLim');
% set(gca,'Ylim',[1e-10 limsy(2)]);
% set(gca,'Xlim',[0 numel(names_transcripts_nonzero_biomass)+1]);
% 
% figure;
% h = barwitherr(std_mca_transcript_pyr_kinase, avg_mca_transcript_pyr_kinase);
% set(gca,'yscale','log');
% limsy = get(gca,'YLim');
% set(gca,'Ylim',[1e-10 limsy(2)]);
% set(gca,'Xlim',[0 numel(names_transcripts_nonzero_pyr_kinase)+1]);
% 
% figure;
% h = barwitherr(std_mca_transcript_lact_dehydrogenase, avg_mca_transcript_lact_dehydrogenase);
% set(gca,'yscale','log');
% limsy = get(gca,'YLim');
% set(gca,'Ylim',[1e-10 limsy(2)]);
% set(gca,'Xlim',[0 numel(names_transcripts_nonzero_lact_dehydrogenase)+1]);
