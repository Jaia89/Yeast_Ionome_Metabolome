%% Import yeast ionome-metabolome multiplex

% cosM layers 
edge_list_MX = readtable('yeast_ionome_metabolome_multiplex_cosM.txt');
n = max(max(edge_list_MX.node_i), max(edge_list_MX.node_j));

% unweighted adjacency matrices
l1_MX = sparse(edge_list_MX.node_i, edge_list_MX.node_j, (edge_list_MX.ionome_ko > 0), n, n);
l2_MX = sparse(edge_list_MX.node_i, edge_list_MX.node_j, (edge_list_MX.ionome_oe > 0), n, n);
l3_MX = sparse(edge_list_MX.node_i, edge_list_MX.node_j, (edge_list_MX.metabolome_aa > 0), n, n);

%% multilink analysis 

Multilinks{1}=triu(ones(size(l1_MX)),1).*l1_MX.*(1-l2_MX).*(1-l3_MX);
Multilinks{2}=triu(ones(size(l1_MX)),1).*(1-l1_MX).*l2_MX.*(1-l3_MX);
Multilinks{3}=triu(ones(size(l1_MX)),1).*(1-l1_MX).*(1-l2_MX).*l3_MX;
Multilinks{4}=triu(ones(size(l1_MX)),1).*l1_MX.*l2_MX.*(1-l3_MX);
Multilinks{5}=triu(ones(size(l1_MX)),1).*l1_MX.*(1-l2_MX).*l3_MX;
Multilinks{6}=triu(ones(size(l1_MX)),1).*(1-l1_MX).*l2_MX.*l3_MX;
Multilinks{7}=triu(ones(size(l1_MX)),1).*l1_MX.*l2_MX.*l3_MX;

for l=1:7
Multilink_Stat(l)=sum(sum(Multilinks{l}));
end
   
%% figure: overall multilink stat  

figure;
bar(Multilink_Stat(4:end));
xticks(1:4);
xticklabels(["ion.ko-ion.oe","ion.ko-met.aa","ion.oe-met.aa","ion.ko-ion.oe-met.aa"]);
xtickangle(45);
set(gca,'box','on','FontSize',16,'Fontname','Arial');
ylabel("# Multilinks");
set(gca,'YScale','log');
title("Ionome-Metabolome Multiplex (cosM)");


%%  Multilayer neighborhood of TPI1

all_nodes = readtable('nodes_list.txt');

AGG_COSM = Multilinks{4} + Multilinks{5} .* 3 + Multilinks{6} .* 5 + Multilinks{7} .* 7;
AGG_COSM = AGG_COSM + AGG_COSM';

G_COSM = graph(AGG_COSM,cellstr(all_nodes.geneName));

gene_id = find(all_nodes.geneName == "YDR050C"); 
neighbourhood = neighbors(G_COSM, gene_id);

H_COSM = subgraph(G_COSM,[gene_id;neighbourhood]);

figure;
h=plot(H_COSM);
h.EdgeLabel=H_COSM.Edges.Weight;
highlight(h,1, 'Nodecolor', 'r');
title('TPI1 (YDR050C) multilink neighborhood (cosM)')
annotation('textbox',...
    [0.2 0.2 0.2 0.2],...
    'String',{'1 ion.ko-ion.oe','3 ion.ko-met.aa',...
    '5 ion.oe-met.aa','7 ion.ko-ion.oe-met.aa'});




