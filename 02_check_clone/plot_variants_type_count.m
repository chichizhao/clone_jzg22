% matlab2022b
% author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
% function: plot the variants type distribution of each sample

mutationsitedistributionsnp = fopen('mutationsitedistributionsnp.csv','r');
mutationsitedistributionindel = fopen('mutationsitedistributionindel.csv','r');
popsnp_01 = table2array(mutationsitedistributionsnp(:,5));
popsnp_11 = table2array(mutationsitedistributionsnp(:,6));
popsnp_00 = table2array(mutationsitedistributionsnp(:,7));
popsnp_dot = table2array(mutationsitedistributionsnp(:,8));
popindel_01 = table2array(mutationsitedistributionindel(:,5));
popindel_11 = table2array(mutationsitedistributionindel(:,6));
popindel_00 = table2array(mutationsitedistributionindel(:,7));
popindel_dot = table2array(mutationsitedistributionindel(:,8));

ax1 = axes('Position',[0.12 0.12 0.18 0.25],'Visible','off');
h =histogram(popsnp_01, 'BinWidth', 1);
h.FaceColor = [0,0.6,0.3];
title('Population Shared heterozygous SNP');
xlabel('number of shared samples');
ylabel('number of sites');
xlim([0 73]);
xticks(0:6:73);
annotation('arrow',[0.26,0.28],[0.16,0.14]);
ax2 = axes('Position',[0.17 0.18 0.10 0.16],'Visible','off');
h = histogram(popsnp_01, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; 
xlim([48 73]);
ylim([0 300]);

ax3 = axes('Position',[0.34 0.12 0.18 0.25],'Visible','off');
h=histogram(popsnp_11, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; % Set to gray color
title('Population Shared homozygous SNP');
xlabel('number of shared samples');
ylabel('number of sites');
xlim([0 73]);
xticks(0:6:72);

ax5 = axes('Position',[0.56 0.12 0.18 0.25],'Visible','off');
h = histogram(popindel_01, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; 
title('Population Shared heterozygous Indel');
xlabel('number of shared samples');
ylabel('number of sites');
xlim([0 73]);
xticks(0:6:73);

annotation('arrow',[0.70,0.72],[0.16,0.14]);
ax4 = axes('Position',[0.61 0.18 0.10 0.16],'Visible','off');
h = histogram(popindel_01, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; 
% set the x axis limit to 48 to 73
xlim([48 73]);
% set the y axis limit to 0 to 500
ylim([0 20]);
ax4 = axes('Position',[0.78,0.12 0.18 0.25],'Visible','off');
h =histogram(popindel_11, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; 
title('Population Shared homozygous Indel');
xlabel('number of shared samples');
ylabel('number of sites');
xlim([0 73]);
xticks(0:6:73);
annotation('arrow',[0.92,0.94],[0.16,0.14]);
ax2 = axes('Position',[0.83 0.18 0.10 0.16],'Visible','off');
h = histogram(popindel_11, 'BinWidth', 1);
h.FaceColor = [0.5 0.5 0.5]; 
% set the x axis limit to 48 to 73
xlim([48 73]);
% set the y axis limit to 0 to 800
ylim([0 20]);