clear all
clc
%time = [0,5,10,15,20];
parentDir = '/Users/reyer/Data/FRET';

Date = {'October_12_2018'};
condition = {'cm','puro','spec','tbox','tet'};
sample ={[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]};

for i =1:length(condition)
    cy5_combined = [];
    cy3_combined = [];
    for j = 1:length(sample{i})
        
        dateDir = strcat([parentDir,'/',Date{1}]);
        if exist(dateDir,'dir')~=7
            mkdir(parentDir,Date{1})
        end
        
        conditionDir = strcat([dateDir,'/',condition{i}]);
        if exist(conditionDir,'dir')~=7
            mkdir(dateDir,condition{i})
        end
        
        filename = strcat([conditionDir,'/sample_00',num2str(sample{i}(j)),'.mat']);
        vbFRETfile = strcat([conditionDir,'/sample_00',num2str(sample{i}(j)),'.txt']);
        
        pathFile = strcat([conditionDir,'/sample',num2str(sample{i}(j)),'_PATH.dat']);
        
        d = getDwell(pathFile);
        D = purifyDwell(d,.1);
        close all
        graph_title = strcat(['Date{1} ', 'condition{i}']);
        figure(1)
        plotTDP(D,24)
        title(graph_title,'FontSize',32)
        colorbar
        set(gcf,'position',[835,883,868,667])
        file1 = strcat([conditionDir,'/sample_00',num2str(sample{i}(j))]);
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-painters','-depsc','-r0')
        set(gcf,'PaperPositionMode','auto')
        print(file1,'-dpng','-r0')
        file1_fig = strcat([conditionDir,'/sample_00',num2str(sample{i}(j)),'.fig']);
        savefig(gcf,file1_fig)
        
        save(filename)
        clearvars d D pathFile cy3c cy5c cy3 cy5 dateDir conditionDir filename vbFRETfile cy3file cy5file good_bleach good_non_bleach total_good
        
        end
    
    vbFRETCombinedfile = strcat([conditionDir,'/CombinedSamples.txt']);
    
    pathFile = strcat([conditionDir,'/CombinedSamples_PATH.dat']);
        
    d = getDwell(pathFile);
    D = purifyDwell(d,.1);
    close all
    graph_title = strcat(['Date{1} ', 'condition{i}']);
    figure(1)
    plotTDP(D,24)
    title(graph_title,'FontSize',32)
    colorbar
    set(gcf,'position',[835,883,868,667])
    file1 = strcat([conditionDir,'/CombinedSamples']);
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([conditionDir,'/CombinedSamples.fig']);
    savefig(gcf,file1_fig)
    
    
    CombinedFilename = strcat([conditionDir,'/CombinedSamples.mat']);
    save(CombinedFilename)
        
    clearvars d D pathFile cy3c cy5c cy3 cy5 dateDir conditionDir filename vbFRETfile cy3file cy5file good_bleach good_non_bleach total_good
        
    
end
