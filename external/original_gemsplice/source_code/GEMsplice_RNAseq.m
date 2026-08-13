if ~exist('Recon_ID_trascript_conversion','var'), load('Recon_ID_trascript_conversion.mat'); end
if ~exist('geni_names','var'), load('geni_names.mat'); end
if ~exist('geni','var'), load('geni.mat'); end
if ~exist('reaction_expression','var'), load('reaction_expression.mat'); end
if ~exist('ID_conversion','var'), load('ID_conversion.mat'); end
if ~exist('ID_conversion_names','var'), load('ID_conversion_names.mat'); end
if ~exist('pos_genes_in_react_expr','var'), load('pos_genes_in_react_expr.mat'); end
if ~exist('fbabreast_trilevel','var'), load('fbabreast_trilevel.mat'); end

%We map the Ensembl and UCSC identifiers onto the Entrez ID notation of Recon
%we used BioMart tool in Ensembl \url{http://www.ensembl.org/biomart/martview/2a3c1aa45a4126aa9947f83d577eee2b}.


if ~exist('geni_transcripts','var')
    geni_transcripts = cell(numel(geni),1);
    
    parfor j=1:numel(geni)
        ix = find(strcmp(Recon_ID_trascript_conversion.EntrezIsoform,geni{j}));
        if ~isempty(ix)
            refseq = Recon_ID_trascript_conversion.Transcript{ix};
            geni_transcripts{j} = refseq2ensembl(refseq,ID_conversion);
        end
        if isempty(geni_transcripts{j}) || strcmp(geni_transcripts{j},'') %either 'ix' was empty, or a conversion refseq2ensembl was not possible
            geni_transcripts{j} = regexprep(geni{j},'\..*',''); %we will use the first part of the Entrez ID of the genes (after the dot character, indicated by \. in the regexp) if it's not an isoform, or when it is an isoform but its transcript ID is not available
        end
    end
end


if ~exist('patients','var'), load('patients.mat'); end

for j = 1:numel(patients)
    patient = patients(j);
    
    if ~isfield(patient,'posizione_gene') || isempty(patient.posizione_gene)
        %Thanks to this IF, the following instructions can indeed be avoided
        %if, for each patient, we know already the position of each the geni.mat
        %arrray into the probe_microarray of the patient. These positions are
        %stored in posizione_gene. If not, we can simply execute the following
        %instructions, ant the posizione_gene array will be computed for each
        %patient.
        
        probe_genes_patient = patient.TranscriptID;
        posizione_gene = cell(length(geni),1);
        
        for i = 1:length(geni)
            
            %this and the following instruction find the locations of the gene 'bXXXX'
            %(genes in the fbamodel) in the array probe_genes (sequence of genes appearing in the
            %microarray data available)
            
            current_gene = geni_transcripts{i};
            
            if numel(current_gene) >= 3 && strcmp(current_gene(1:3), 'ENS')   %then this Recon's gene was annotated with its transcript(s) ID
                matches = strfind(probe_genes_patient,current_gene);
                posizione_gene{i} = find(~cellfun('isempty', matches));
                
                if isempty(matches) %if we were unsuccessful, we try to search the same transcript but using the UCSC identifier instead of the Ensembl ID, because some patients' probes are UCSC
                    current_gene_ucsc = ensembl2ucsc(current_gene, ID_conversion);
                    matches = strfind(probe_genes_patient,current_gene_ucsc);
                    posizione_gene{i} = find(~cellfun('isempty', matches));
                end
                
            else    %then this Recon's gene was NOT annotated with its transcript, we need to use its gene Entrez ID
                gene_ensembl = entrez2ensembl(current_gene, ID_conversion_names);
                matches = many2many_strfind(probe_genes_patient,gene_ensembl);
                posizione_gene{i} = matches;
                
                if isempty(matches) %if we were unsuccessful, we try to search the same transcript but using the UCSC identifier instead of the Ensembl ID, because some patients' probes are UCSC
                    current_gene_ucsc = entrez2ucsc(current_gene, ID_conversion_names);
                    current_gene_ucsc = regexprep(current_gene_ucsc,'\..*',''); % like we did with the Entrez IDs, we will discard everything after the dot. Indeed, in something like uc011mwn.1, the .1 represents the revision number of the transcript (and DOES NOt indicate the transcript isoform). When a new version of UCSC Genes is released, a transcript like uc011mwn.1 could possibly remain the same, it could become uc011mwn.2, it could receive a new transcript ID entirely or it could disappear altogether from the new version of UCSC Genes.
                    matches = many2many_strfind(probe_genes_patient,current_gene_ucsc);
                    posizione_gene{i} = matches;
                end
                
            end
            
        end
        
        patients(j).posizione_gene = posizione_gene;
    end
end

for j = 1:numel(patients)
    j
    patient = patients(j);
    posizione_gene = patient.posizione_gene;
    
    microarray_data = patient.Exp;
    dati_micro = zeros(length(geni),1);
    parfor i=1:length(geni)
        if isempty(posizione_gene{i})
            dati_micro(i) = 1; %if data is available, put here the average gene expression of that gene (across condition)
        else
            dati_micro(i) = mean(microarray_data(posizione_gene{i})); %if the gene is mapped to more probe locations (as it happens with the isoforms for which we don't have annotation and we had to use the gene name), we average the gene expression level in these probe locations
        end
        
    end
    
    [out1 out2] = evaluate_objective(dati_micro,3,numel(geni),fbabreast,geni,reaction_expression, pos_genes_in_react_expr);
    patients_results(j).v1_table = out1;
    patients_results(j).f1 = out2;
    out2
    
    
end




