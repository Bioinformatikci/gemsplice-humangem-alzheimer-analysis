%% PLOT INVASIVE_VS_UNLABELED FIGURE FOR PUBLICATION
close all

if ~exist('fbabreast_trilevel','var'), load('fbabreast_trilevel.mat'); end
if ~exist('patients_results','var'), load ('patients_results.mat'); end
if ~exist('patients','var'), load('patients.mat'); end
if ~exist('geni','var'), load('geni.mat'); end
if ~exist('reaction_expression','var'), load('reaction_expression.mat'); end
if ~exist('pos_genes_in_react_expr','var'), load('pos_genes_in_react_expr.mat'); end




%% PLOT THE 3D OUTPUT OF TRILEVEL FBA; PLOT 3D SCATTER PLOT WITH SHADOWS ON XY, XZ, YZ PLANES
%%
ix_atp = find(strcmp(fbabreast.rxnNames,'ATP exchange'));
ix_glut = find(strcmp(fbabreast.rxnNames,'L-Glutamine exchange'));
ix_glc = find(strcmp(fbabreast.rxnNames,'D-Glucose exchange'));

for i=1:31
    biomass(i) = patients_results(i).f1(1);
    pyr_kin(i) = patients_results(i).f1(2);
    lac_dehy(i) = patients_results(i).f1(5);
    names{i} = patients(i).CancerType;
    %fourth(i) = patients_results(i).v1_table(ix_glut,5);
end

final_table = [names' num2cell(biomass') num2cell(pyr_kin') num2cell(lac_dehy')]


figure
X = biomass; Y = pyr_kin; Z = lac_dehy;
scatter3(X,Y,Z);
limsx = get(gca,'XLim');
limsy = get(gca,'YLim');
limsz = get(gca,'ZLim');
hold on;
plot3(X,Y,limsz(1)*ones(size(X)),'r.');
plot3(limsx(1)*ones(size(X)),Y,Z,'g.');
plot3(X,limsy(1)*ones(size(X)),Z,'k.');
set(gca,'Color',[0.8 0.8 0.8]);




%% PLOT SEQUENCE OF PLOTS AND PROJECTIONS
%%
index_invasive = find(not(cellfun('isempty', strfind(names,'invasive'))));
index_unlabeled = setdiff(1:numel(patients),index_invasive);

figure
a = cmap('teal'); colormap(a)

make_it_tight = false;
subplot = @(m,n,p) subtightplot (m, n, p, [0.11 0.11], 0.2, [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

points = [X' Y' Z'];

subplot(4,1,1);
x = points(index_invasive,1);
y = points(index_invasive,2);
z = points(index_invasive,3);
pl3 = plot3(x,y,z,'v','Color',rgb('Teal'));
set(pl3,'MarkerSize',1); 
%plot(x,y,'x','Color',[143/255 0/255 255/255])
hold on
xa = points(index_unlabeled,1);
ya = points(index_unlabeled,2);
za = points(index_unlabeled,3);
pl3 = plot3(xa,ya,za,'*','Color',rgb('deeppink'));
set(pl3,'MarkerSize',1); 

xlabel('Biomass [h^{-1}]')
ylabel('Pyruvate kinase [mmolh^{-1}gDW^{-1}]')
zlabel('Lactate dehydrogenase [mmolh^{-1}gDW^{-1}]')

%lg = legend('invasive','unlabeled','Location','NorthWest'); legend('boxoff');
%lg.Position=lg.Position + [-0.05 0.08 0 0];

[lg, hobj, hout, mout] = legend('invasive','unlabeled','Location','NorthWest'); legend('boxoff');
M = findobj(hobj,'type','Text');  %this reduces the huge length of the line representing the marker in the legend
M(1).FontSize = M(1).FontSize-5;
M(2).FontSize = M(2).FontSize-5;
M(1).Position(1) = M(1).Position(1)-0.15;
M(2).Position(1) = M(2).Position(1)-0.15;

lg.Location = 'NorthWest';
lg.Position(1)=lg.Position(1)-0.03;
lg.Position(2)=lg.Position(2)+0.074;

grid on

limsx = get(gca,'XLim');
limsy = get(gca,'YLim');
limsz = get(gca,'ZLim');

%for acetate:
set(gca,'box','on');
%set(gca,'dataaspectratio',[0.2 0.3 1.6],'projection','perspective','box','on') %cambiare i parametri/numeri per fissare il ratio, assi con numeri piu' piccoli devono avere ratio piccoli
%for succinate:
%set(gca,'dataaspectratio',[0.02 0.3 0.8],'projection','perspective','box','on') %cambiare i parametri/numeri per fissare il ratio, assi con numeri piu' piccoli devono avere ratio piccoli

h = rotate3d;
set(h,'ActionPostCallback',@align_axislabels) %to see the effect of the rotation of the labels, please click and interact with the output figure in Matlab (just a slight rotation is enough to update the figure and process the callbacks)

%%
subplot(4,1,2);
C = z;
p1 = scatter(x,y,20,C,'v');
hold on
non_dom_points = non_domination_sort_mod([-x +y],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2) * [-1 0; 0 1];     %keeps only the coordinates of the nondominated points and restore original signs
non_dom_points = sortrows(non_dom_points,1);
ln1 = plot(non_dom_points(:,1),non_dom_points(:,2),'-', 'Color',rgb('Teal'),'DisplayName','invasive');
set(ln1,'MarkerSize',1); 

% instructions for colorbars
hc=colorbar;
ylabel(hc,'Lactate dehydrogenase [mmolh^{-1}gDW^{-1}]')
hc.FontSize=hc.FontSize-5;
hc.Color=[153/255 0/255 0/255];
ax = gca;
axpos = ax.Position;
hc.Position(3) = hc.Position(3)*0.5;
ax.Position = axpos;
ax.Position(3)=1.05*ax.Position(3);

hold on
C = za;
p11 = scatter(xa,ya,10,C,'*');
non_dom_points = non_domination_sort_mod([-xa +ya],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2) * [-1 0; 0 1];     %keeps only the coordinates of the nondominated points and restore original signs
non_dom_points = sortrows(non_dom_points,1);
ln2 = plot(non_dom_points(:,1),non_dom_points(:,2),'-','Color',rgb('deeppink'),'DisplayName','unlabeled');
set(ln2,'MarkerSize',1); 

xlabel('Biomass [h^{-1}]')
ylabel('Pyruvate kinase [mmolh^{-1}gDW^{-1}]')

% aux = xlim;
% xlim([0 1.2]);
% aux = ylim;
% ylim([0 aux(2)]);

[lg, hobj, hout, mout] = legend([ln1 ln2]); legend('boxoff');
M = findobj(hobj,'type','Line');  %this reduces the huge length of the line representing the marker in the legend
M(1).XData(1) = 0.3;
M(3).XData(1) = 0.3;
M = findobj(hobj,'type','Text');  %this reduces the huge length of the line representing the marker in the legend
M(1).FontSize = M(1).FontSize-5;
M(2).FontSize = M(2).FontSize-5;
lg.Location = 'NorthEast';
lg.Position(1)=lg.Position(1)+0.025;
lg.Position(2)=lg.Position(2)+0.045;

%%
subplot(4,1,3);
C = x;
p2 = scatter(y,z,20,C,'v');
hold on
non_dom_points = non_domination_sort_mod([+y +z],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2);     %keeps only the coordinates of the nondominated points
non_dom_points = sortrows(non_dom_points,1);
ln1 = plot(non_dom_points(:,1),non_dom_points(:,2),'-', 'Color',rgb('Teal'),'DisplayName','invasive');

% instructions for colorbars
hc=colorbar;
ylabel(hc,'Biomass [h^{-1}]')
hc.FontSize=hc.FontSize-5;
hc.Color=[153/255 0/255 0/255];
ax = gca;
axpos = ax.Position;
hc.Position(3) = hc.Position(3)*0.5;
ax.Position = axpos;
ax.Position(3)=1.05*ax.Position(3);

hold on
C = xa;
p22 = scatter(ya,za,10,C,'*');
non_dom_points = non_domination_sort_mod([+ya +za],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2);     %keeps only the coordinates of the nondominated points
non_dom_points = sortrows(non_dom_points,1);
ln2 = plot(non_dom_points(:,1),non_dom_points(:,2),'-','Color',rgb('deeppink'),'DisplayName','unlabeled');
set(ln2,'MarkerSize',1); 

xlabel('Pyruvate kinase [mmolh^{-1}gDW^{-1}]')
ylabel('Lactate dehydrogenase [mmolh^{-1}gDW^{-1}]')
% aux = xlim;
% xlim([0 aux(2)]);
% aux = ylim;
% ylim([0 aux(2)]);

[lg, hobj, hout, mout] = legend([ln1 ln2]); legend('boxoff');
M = findobj(hobj,'type','Line');  %this reduces the huge length of the line representing the marker in the legend
M(1).XData(1) = 0.3;
M(3).XData(1) = 0.3;
M = findobj(hobj,'type','Text');  %this reduces the huge length of the line representing the marker in the legend
M(1).FontSize = M(1).FontSize-5;
M(2).FontSize = M(2).FontSize-5;
lg.Location = 'NorthEast';
lg.Position(1)=lg.Position(1)+0.025;
lg.Position(2)=lg.Position(2)+0.045;

%%
subplot(4,1,4);
C = y;
p3 = scatter(x,z,20,C,'v');
hold on
non_dom_points = non_domination_sort_mod([-x +z],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2) * [-1 0; 0 1];     %keeps only the coordinates of the nondominated points and restore original signs
non_dom_points = sortrows(non_dom_points,1);
ln1 = plot(non_dom_points(:,1),non_dom_points(:,2),'-', 'Color',rgb('Teal'),'DisplayName','invasive');

% instructions for colorbars
hc=colorbar;
ylabel(hc,'Pyruvate kinase [mmolh^{-1}gDW^{-1}]')
hc.FontSize=hc.FontSize-5;
hc.Color=[153/255 0/255 0/255];
ax = gca;
axpos = ax.Position;
hc.Position(3) = hc.Position(3)*0.5;
ax.Position = axpos;
ax.Position(3)=1.05*ax.Position(3);


hold on
C = ya;
p33 = scatter(xa,za,10,C,'*');
non_dom_points = non_domination_sort_mod([-xa +za],2,0);  %ranks points according to non_domination (remember that non_domination_sort_mod is set to look for min-min by default
non_dom_points = non_dom_points(find(non_dom_points(:,3)==1),1:2) * [-1 0; 0 1];     %keeps only the coordinates of the nondominated points and restore original signs
non_dom_points = sortrows(non_dom_points,1);
ln2 = plot(non_dom_points(:,1),non_dom_points(:,2),'-','Color',rgb('deeppink'),'DisplayName','unlabeled');
set(ln2,'MarkerSize',1); 

xlabel('Biomass [h^{-1}]');
ylabel('Lactate dehydrogenase [mmolh^{-1}gDW^{-1}]');

% aux = xlim;
% xlim([0 1.2]);
% aux = ylim;
% ylim([0 aux(2)]);

[lg, hobj, hout, mout] = legend([ln1 ln2]); legend('boxoff');
M = findobj(hobj,'type','Line');  %this reduces the huge length of the line representing the marker in the legend
M(1).XData(1) = 0.3;
M(3).XData(1) = 0.3;
M = findobj(hobj,'type','Text');  %this reduces the huge length of the line representing the marker in the legend
M(1).FontSize = M(1).FontSize-5;
M(2).FontSize = M(2).FontSize-5;
lg.Location = 'SouthEast';
lg.Position(1)=lg.Position(1)+0.025;
lg.Position(2)=lg.Position(2)-0.055;

if exist('exportfig','file')
    exportfig(gcf, 'inset.pdf', 'color', 'cmyk', 'Width', '10', 'Height', '30', 'FontMode', 'scaled', 'FontSize', '0.8' );
end








%% PLOT THE AVERAGE FLUX IN THE SUBSYSTEMS
%% In mean_flux_pathway, we compute the average absolute flux rate through every pathway and for every patient

dati_micro = ones(1, numel(geni));
[out1 out2] = evaluate_objective(dati_micro,3,numel(geni),fbabreast,geni,reaction_expression, pos_genes_in_react_expr);
[out2(1) out2(2) out2(5)]

[subs_unique, ia, rxn_to_subs_unique] = unique(fbabreast.subSystems);
mean_flux_pathway = zeros(numel(subs_unique),numel(patients));
initial_mean_flux_pathway = zeros(numel(subs_unique),1);
pathway_occurrences = histc(rxn_to_subs_unique, 1:numel(subs_unique));

for i = 1:numel(subs_unique)
    ixs_reactions_in_this_pathway = find(rxn_to_subs_unique == i);
    %initial_mean_flux_pathway(i) = mean(abs(cell2mat(out1(ixs_reactions_in_this_pathway,5))));
    initial_mean_flux_pathway(i) = mean(cell2mat(out1(ixs_reactions_in_this_pathway,5)));
    
    for j = 1:numel(patients)
        %mean_flux_pathway(i,j) = mean(abs(cell2mat(patients_results(j).v1_table(ixs_reactions_in_this_pathway,5)))); %the fifth column of v1_table is the one we are interested in, as it is v1min_max (we maximize the third objective and minimize the second
        mean_flux_pathway(i,j) = mean(cell2mat(patients_results(j).v1_table(ixs_reactions_in_this_pathway,5))); %the fifth column of v1_table is the one we are interested in, as it is v1min_max (we maximize the third objective and minimize the second
    end
end


%% for each subsystem, we compute mean and standard deviation *across* patients/cancers
avg_by_pathway = mean(mean_flux_pathway');
std_by_pathway = std(mean_flux_pathway');
ixs_nonzero = find(avg_by_pathway ~= 0);
names_nonzero = subs_unique(ixs_nonzero);
names_nonzero(find(isempty(names_nonzero)))='';

%we normalize the pathway fluxes with respect the "initial" pathway flux in a general cancer model.
normalized_mean_flux_pathway = bsxfun (@rdivide, mean_flux_pathway, initial_mean_flux_pathway); %divide (element by element) all columns of mean_flux_pathway by the array initial_mean_flux_pathway

normalized_avg_by_pathway = zeros(1,numel(avg_by_pathway));
normalized_std_by_pathway = zeros(1,numel(std_by_pathway));

normalized_avg_by_pathway(ixs_nonzero) = avg_by_pathway(ixs_nonzero)./initial_mean_flux_pathway(ixs_nonzero)';
normalized_std_by_pathway(ixs_nonzero) = std(normalized_mean_flux_pathway(ixs_nonzero,:)');

derivative_avg_by_pathway = nan(1,numel(subs_unique));
derivative_avg_by_pathway(ixs_nonzero) = (avg_by_pathway(ixs_nonzero)-initial_mean_flux_pathway(ixs_nonzero)')./avg_by_pathway(ixs_nonzero);

avg_by_pathway_nonzero = avg_by_pathway(ixs_nonzero);
std_by_pathway_nonzero = std_by_pathway(ixs_nonzero);
normalized_avg_by_pathway_nonzero = normalized_avg_by_pathway(ixs_nonzero);
normalized_std_by_pathway_nonzero = normalized_std_by_pathway(ixs_nonzero);

%save('export_to_R.mat','avg_by_pathway_nonzero','std_by_pathway_nonzero','normalized_avg_by_pathway_nonzero','normalized_std_by_pathway_nonzero','names_nonzero');


%% we repeat everything in the last section but now applied to the two
%%subcases, i.e. we split the 31 samples in invasive and unlabeled and
%%repeat the analysis
unlabeled_avg_by_pathway = mean(mean_flux_pathway(:,index_unlabeled),2)';
unlabeled_std_by_pathway = std(mean_flux_pathway(:,index_unlabeled),0,2)';
invasive_avg_by_pathway = mean(mean_flux_pathway(:,index_invasive),2)';
invasive_std_by_pathway = std(mean_flux_pathway(:,index_invasive),0,2)';

unlabeled_derivative_avg_by_pathway = nan(1,numel(subs_unique));
invasive_derivative_avg_by_pathway = nan(1,numel(subs_unique));
unlabeled_derivative_avg_by_pathway(ixs_nonzero) = (unlabeled_avg_by_pathway(ixs_nonzero)-initial_mean_flux_pathway(ixs_nonzero)')./unlabeled_avg_by_pathway(ixs_nonzero);
invasive_derivative_avg_by_pathway(ixs_nonzero) = (invasive_avg_by_pathway(ixs_nonzero)-initial_mean_flux_pathway(ixs_nonzero)')./invasive_avg_by_pathway(ixs_nonzero);

pathway_table = [subs_unique(:) num2cell(unlabeled_derivative_avg_by_pathway(:)) num2cell(invasive_derivative_avg_by_pathway(:))]

ixs_nonzero_reduced = setdiff(ixs_nonzero, [find(isnan(unlabeled_derivative_avg_by_pathway)) find(isnan(invasive_derivative_avg_by_pathway)) find(strcmp(subs_unique,'exchange'))  find(strcmp(subs_unique,'Transport, Nuclear'))] ); %for visualizaion purposes, we also remove NaN, transport and exchange pathways 
names_nonzero_reduced = subs_unique(ixs_nonzero_reduced);
x = unlabeled_derivative_avg_by_pathway(ixs_nonzero_reduced);
y = invasive_derivative_avg_by_pathway(ixs_nonzero_reduced);
r = pathway_occurrences(ixs_nonzero_reduced);
c = log(std_by_pathway(ixs_nonzero_reduced)'./sqrt(pathway_occurrences(ixs_nonzero_reduced))); %log (for visualization purposes) of the normalized_error = st_deviation/square_root_of_size_of_the_pathway (this resembles the standard error of the mean, which however has the sample size at the denominator, namely the number of samples used to compute the mean)
fig = figure;
scatter(x, y, r, c, 'filled', 'MarkerEdgeColor', 'k');
xlabel('Unlabeled - normalized relative pathway flux change');
ylabel('Invasive - normalized relative pathway flux change');
colormap(parula(256))

% instructions for horizontal colorbars
ax=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
pos=get(gca,'pos');
hc = colorbar('location','northoutside','position',[pos(1) pos(2)+pos(4)+0.03 pos(3) 0.02]);
set(hc,'xaxisloc','top');
ylabel(hc,'log(NE)');
hc.FontSize=hc.FontSize-2;
hc.Color=[153/255 0/255 0/255];
%set(hc,'YTick', 10:10:max(C));

[spearman_corr, pval_spearman] = corr(x',y','type','Spearman')
[pearson_corr, pval_pearson] = corr(x',y','type','Pearson')
%lsline
refline(1,0)    % plot the y=x line (slope 1, intercept 0)

%% SAVE ORIGINAL FIGURES, THEN POLISHED FIGURES AND POINT INFORMATION
outdir = fullfile('plot_outputs','invasive_vs_unlabeled');
save_original_gemsplice_plot_outputs(outdir, 'Invasive vs Unlabeled', ...
    patients, patients_results, final_table, pathway_table, names, biomass, pyr_kin, lac_dehy);


% figure;
% h = barwitherr(std_by_pathway(ixs_nonzero), avg_by_pathway(ixs_nonzero));
% set(gca,'yscale','log');
% limsy = get(gca,'YLim');
% set(gca,'Ylim',[1e-10 limsy(2)]);
% set(gca,'Xlim',[0 numel(ixs_nonzero)+1]);
% 
% figure;
% h = barh(avg_by_pathway(ixs_nonzero));
% set(gca,'xscale','log');
% limsy = get(gca,'XLim');
% set(gca,'Xlim',[-10 limsy(2)]);
% set(gca,'Ylim',[0 numel(ixs_nonzero)+1]);
% 
% exp1 = avg_by_pathway(ixs_nonzero);
% exp2 = normalized_avg_by_pathway(ixs_nonzero);
% exp3 = std_by_pathway(ixs_nonzero);
% exp4 = normalized_std_by_pathway(ixs_nonzero); 
% 
% 
% 
% figure;
% h = barwitherr(log(normalized_std_by_pathway), log(normalized_avg_by_pathway));
% %set(gca,'yscale','log');
% limsy = get(gca,'YLim');
% set(gca,'Ylim',[-10 limsy(2)]);
% set(gca,'Xlim',[0 101]);
% 
% figure;
% h = barh(log(normalized_avg_by_pathway));
% set(gca,'xscale','log');
% limsy = get(gca,'XLim');
% set(gca,'Xlim',[-10 limsy(2)]);
% set(gca,'Ylim',[0 102]);




