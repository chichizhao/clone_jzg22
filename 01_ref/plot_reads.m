% matlab2022b
% author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
% plot the distribution of reads length and quality score
% Input: read.csv
% the function is used to plot the distribution of reads length and quality score
% it has three subplots, the first one is the distribution of reads length
% the second one is the distribution of quality score
% the third one is the scatter plot of reads length and quality score

read_length = Reads.VarName1 ;
read_quality = reads.VarName2 ;
subplot(2,2,1)
% use log10 to make the distribution more clear
histogram(read_length,200)
set(gca,'Position',[0.1 0.7 0.6 0.15])
% note that the x axis is log10 of reads length
set(gca,'xtick',0:1:6)
set(gca,'yaxislocation','right')
%ylabel('Number of reads')
% xlabel('log10(Reads length)')
% set the y label to be on the right


% set x axis form 0 to 6
set(gca,'xlim',[0 200000])
set (gca,"XTickLabel",[])
subplot(2,2,4)
% plot the distribution of quality score
score = histogram(read_quality,300)
set(gca,'Position',[0.7 0.1 0.15 0.6])
% reverse the x axis
set(gca,'xdir','reverse')
set(gca,'xtick',0:10:80)
set(gca,'ytick',0:10000:20000)
set(gca,"XTickLabel",[])
set(gca,'yaxislocation','left')
ylabel('Number of reads')
camroll(270)
set(gca,'xlim', [5 35])


subplot(2,2,3)
scatter(read_length,read_quality,2,"filled","MarkerFaceColor","flat")
xlabel('Reads length')
ylabel('Quality score')
set(gca,'xlim', [0 200000])
    % set y axis form 0 to 6,start from 2
set(gca,'ytick',0:10:80)
set(gca,'ylim', [5 35])
set(gca,'xtick',0:20000:200000)
% set the subplot to be tight
set(gca,'Position',[0.1 0.1 0.599 0.599])