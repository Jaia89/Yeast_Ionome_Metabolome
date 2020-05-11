%% Import yeast ionome-metabolome multiplex

% cosM network only

MX_COSM=cell(1,3);

MX_COSM{1} = importdata('PCC_multiplex/PCC_ionome_ko.txt');
MX_COSM{2} = importdata('PCC_multiplex/PCC_ionome_oe.txt');
MX_COSM{3} = importdata('PCC_multiplex/PCC_metabolome_aa.txt');

%% multilink analysis 

Multilinks_COSM=cell(1,7);
Multilink_Stat_COSM=zeros(1,7);

% 100K top-scored links
thr = 100000;     
  
L1_COSM = ( MX_COSM{1} <= thr) .* (MX_COSM{1} >0 );
L2_COSM = ( MX_COSM{2} <= thr) .* (MX_COSM{2} >0 );
L3_COSM = ( MX_COSM{3} <= thr) .* (MX_COSM{3} >0 );

Multilinks_COSM{1}=triu(ones(size(L1_COSM)),1).*L1_COSM.*(1-L2_COSM).*(1-L3_COSM);
Multilinks_COSM{2}=triu(ones(size(L1_COSM)),1).*(1-L1_COSM).*L2_COSM.*(1-L3_COSM);
Multilinks_COSM{3}=triu(ones(size(L1_COSM)),1).*(1-L1_COSM).*(1-L2_COSM).*L3_COSM;
Multilinks_COSM{4}=triu(ones(size(L1_COSM)),1).*L1_COSM.*L2_COSM.*(1-L3_COSM);
Multilinks_COSM{5}=triu(ones(size(L1_COSM)),1).*L1_COSM.*(1-L2_COSM).*L3_COSM;
Multilinks_COSM{6}=triu(ones(size(L1_COSM)),1).*(1-L1_COSM).*L2_COSM.*L3_COSM;
Multilinks_COSM{7}=triu(ones(size(L1_COSM)),1).*L1_COSM.*L2_COSM.*L3_COSM;

for l=1:7
Multilink_Stat_COSM(l)=sum(sum(Multilinks_COSM{l}));
end
   
%% figure: overall multilink stat  

figure;
bar(Multilink_Stat_COSM(4:end));
xticks(1:4);
xticklabels(["ion.ko-ion.oe","ion.ko-met.aa","ion.oe-met.aa","ion.ko-ion.oe-met.aa"]);
xtickangle(45);
set(gca,'box','on','FontSize',16,'Fontname','Arial');
ylabel("# Multilinks");
set(gca,'YScale','log');
title("Ionome-Metabolome Multiplex (cosM)");


%%  Multilayer neighborhood of TPI1

all_nodes = readtable('nodes_list.txt');

AGG_COSM = Multilinks_COSM{4} + Multilinks_COSM{5} .* 3 + Multilinks_COSM{6} .* 5 + Multilinks_COSM{7} .* 7;
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




