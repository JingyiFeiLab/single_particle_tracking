file = '/Users/reyer/Data/STORM/81520/';
close all
X = {'SgrS t=30, -ptsG','SgrS t=60, +ptsG','RyhB t=30, -sodB','RyhB t=60, +sodB'};
%X = reordercats(X,{'WT rne, WT SgrS','WT rne, WT RyhB','WT rne, mut SgrS','rne701, WT SgrS'});
%y = [3.308433723 0.9733033419 1.104049937 0.1572210462; 4.186459433 0.966771806 2.634121304 1.033102583; 0 0 1.370988622 0.7119615517; 5.445199532 1.031210025 0.8030783902 0.5380187228];
y = [mean([0.9508104946 0.9688955166]) mean([0.9863893134 0.9811016599]); mean([1.037200795 0.9698600346]) mean([0.9482083104 1.02853082]); mean([0.9730083306 1.003876099]) mean([1.125412 0.9313498963]); mean([1.015505189 1.080646737]) mean([1.03780685 0.9195394716])];
err  = [std([0.9508104946 0.9688955166]) std([0.9863893134 0.9811016599]); std([1.037200795 0.9698600346]) std([0.9482083104 1.02853082]); std([0.9730083306 1.003876099]) std([1.125412 0.9313498963]); std([1.015505189 1.080646737]) std([1.03780685 0.9195394716])];
b = bar(y, 'FaceColor','flat');
set(b(1),'FaceColor','b')
set(b(2),'FaceColor','r')
hold on
errorbar([.85 1.15; 1.85 2.15; 2.85 3.15; 3.85 4.15],y,err,'LineStyle','none','Color','k','LineWidth',1)
set(gca,'xtick',[1:4],'xticklabel',X,'FontSize',14)
lgd = legend('Cytoplasm','DAPI');
ylabel('Enrichment','FontSize',14)
title('sRNA Localization','FontSize',16)
set(gca, 'FontSize', 10)
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
lgd.Location = 'northeastoutside';
set(gcf,'position',[318,235,556,247])
%set(gca, 'FontWeight', 'bold')
file1 = strcat([file,'STORM_summary']);
set(gcf,'PaperPositionMode','auto')
print(file1,'-dpng','-r0')
print(file1,'-painters','-depsc','-r0')
file1_fig = strcat([file,'STORM_summary.fig']);
savefig(gcf,file1_fig)
