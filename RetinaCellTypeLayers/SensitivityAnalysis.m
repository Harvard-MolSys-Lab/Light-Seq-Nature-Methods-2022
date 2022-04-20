% Import list of bipolar marker genes detected by smFISH in (West et al., 2022)
Bipolar_Markers=[
    {'Pcdh17' }
    {'Wls'    }
    {'Neto1'  }
    {'Chrm2'  }
    {'Erbb4'  }
    {'Grik1'  }
    {'Slitrk5'}
    {'Col11a1'}
    {'Meis2'  }
    {'Lrrtm1' }
    {'Vsx1'   }
    {'Grm6'   }
    {'Igfn1'  }
    {'Cpne9'  }
    {'Prkca'  }
    {'Vsx2'   }]
%% Import Drop-Seq counts per layer of all marker genes
D_counts=readtable('./DropSeq/DropSeqLayerFrequenciesMinusBC1B.csv');

% This is after removing the BC1B cluster, since this subtype was not captured in the Light-Seq
% data due to it's positioning within the amacrine cell layer rather than
% the bipolar cell layer. 
% See Figure S6 of (Shekhar et al., 2016).
D_counts(ismember(D_counts.Gene,Bipolar_Markers),8:13)

% Subset for only the marker genes listed above for smFISH comparison
for i=1:16
    tD(i,:)=D_counts(ismember(D_counts.Gene,Bipolar_Markers(i)),[1,8:13])
end

%% Import Light-Seq counts per layer of all marker genes
L_counts=readtable('./LightSeq/ReorderedLightSeq.csv');

% Subset for only the marker genes listed above for smFISH comparison
for i=1:16
t(i,:)=L_counts(ismember(L_counts.Gene,Bipolar_Markers(i)),[1,6:9])
end

%% Import smFISH data from (West et al. 2022)

% Mean TPC from FISH data, pooled across the entire bipolar cell population
% at their measured ratios pulled from West et al. 2022. Cells of subtype identity BC1B was
% removed from the data because their cell body location is below (basal to) the regions targeted by Light-Seq

MeanTPC_FISHarray=[1.235267575	0.066871221	3.043768613	0.709322263	0.844237885	20.57738471	1.252233553	0.68730259	1.970670517	1.084739644	7.948650844	35.50482808	4.722858948	0.774027615	19.43272268	20.0822128]

% Number of BPs within Light-Seq ROIs for each replicate (A2, B1, B3, and
% C2, respectively)
n_BPs=[233 238 209 188];

%no BC1B
n_BPs_D=[2902	3244	3466	3631	6535	5869];

% Define expected transcript number as [mean transcripts per cell by FISH]*[# Bipolar Cells captured for sequencing]
Transcripts_Expected=MeanTPC_FISHarray.*n_BPs';
Transcripts_Expected_D=MeanTPC_FISHarray.*n_BPs_D';


%% Import gene lengths (based on EnsembleIDs)
l=readtable('./Gene_Lengths.xlsx')
% Pulled lengths of bipolar markers from the featurecounts output file
length=l.Length

%% Compute Sensitivity

READS=table2array(removevars(t,'Gene'));
READS_D=table2array(removevars(tD,'Gene'));

% Define sensitivity as [Raw Reads] / [Expected Reads based on FISH] * 100
Sensitivity=(READS./Transcripts_Expected'*100)'
Sensitivity_D=(READS_D./Transcripts_Expected_D'*100)'


%% Create boxplot Light-Seq vs. Drop-Seq
% Each point represents the estimated sensitivity of detecting one bipolar
% subtype marker gene gene in one replicate. This includes 16 genes x 4
% replicates for Light-Seq and 16 genes x 6 replicates for Drop-Seq

figure;
subplot(1,2,1)
% Plot Light-Seq sensitivity
boxplot(Sensitivity(:),'Colors',[0 0 0],'whisker', 1.5);
hold on; scatter(ones(size(Sensitivity(:))),Sensitivity(:),'filled','MarkerFaceColor',[0.5 0.6 0.7],'MarkerFaceAlpha',0.4,'jitter',0.02)
title('Sensitivity Compared to smFISH')
set(gca,'FontSize',16)
ylim([-0.1 25])
ylabel('Sensitivity (%)')
set(gcf,'color','w')
xticklabels('Light-Seq')

% Plot Drop-Seq sensitivity
subplot(1,2,2)
boxplot(Sensitivity_D(:),'Colors',[0 0 0],'whisker', 1.5);
hold on; scatter(ones(size(Sensitivity_D(:))),Sensitivity_D(:),'filled','MarkerFaceColor',[0.5 0.6 0.7],'MarkerFaceAlpha',0.4,'jitter',0.02)
set(gca,'FontSize',16)
ylabel('Sensitivity (%)')
ylim([-0.1 26])
set(gcf,'color','w')
xticklabels('Drop-Seq')
title('Sensitivity Compared to smFISH')

%% Re-index for ordering genes by increasing length
index=[1:1:16]'
li=[length index]
sortedi= sortrows(li,1,'ascend');
sind=sortedi(:,2)
slen=3

%% Define viridis colormap based on length
inf=min(length)
sup=max(length)

% Uses viridis.m from: https://github.com/moffat/matlab/blob/master/colormaps/viridis.m
VC=viridis((sup-inf+1))
VC=flipud(VC)
colorind=[1:1:size(VC,1)]

% Normalize lengths to map to viridis colormap
lnorm2 = round(normalize(length, 'range', [1 256]))
figure;heatmap(lnorm2);colormap(VC)
 
% Construct colormap
 for i=1:16 
     VCD(i,:)=VC(lnorm2(i),:)
 end
 

%% Plot per gene for Light-Seq
figure;
for j=1:16
    i=sind(j)
hold on; scatter(j*ones(size(Sensitivity(:,i))),Sensitivity(:,i),70,'filled','CData',VCD(i,:))
%hold on; scatter(j,mean(Sensitivity(:,i)),150,'filled','CData',VCD(i,:))
%hold on; scatter(j,mean(Sensitivity(:,i)),length(i)/100,'filled','CData',VCD(i,:))
hold on; errorbar(j,mean(Sensitivity(:,i)),std(Sensitivity(:,i)),'LineStyle','None','Color',VCD(i,:),'LineWidth',1)
end
xticks([1:1:16])
xticklabels(t.Gene(sind))
xlim([0 17])
ylim([0 45])
set(gcf,'color','w')
ylabel('Sensitivity (%)')

%% Plot per gene for Drop-Seq
figure
for j=1:16
    i=sind(j)
hold on; scatter(j*ones(size(Sensitivity_D(:,i))),Sensitivity_D(:,i),70,'filled','CData',VCD(i,:))
hold on; errorbar(j,mean(Sensitivity_D(:,i)),std(Sensitivity_D(:,i)),'LineStyle','None','Color',VCD(i,:),'LineWidth',1)
end

xticks([1:1:16])
xticklabels(t.Gene(sind))
xlim([0 17])
ylim([0 45])
set(gcf,'color','w')
ylabel('Sensitivity (%)')


%% Plot colormap
figure;heatmap(1:1:size(VC),'GridVisible','off')
colormap(VC)

%% Plot sensitivity difference between Light-Seq and Drop-Seq

ML=mean(Sensitivity)
MD=mean(Sensitivity_D)
 
figure;
plot([0,17], [0,0],'Color',[0.5 0.5 0.5])
err=sqrt((std(Sensitivity).^2/size(Sensitivity,1))+(std(Sensitivity_D).^2/size(Sensitivity_D,1)))

for j=1:16
    i=sind(j)
    hold on
scatter([j],(ML(i)-MD(i)),70,VCD(i,:),'filled')
errorbar([j],(ML(i)-MD(i)),err(i),'LineStyle','None','Color',VCD(i,:),'LineWidth',1)
end

set(gcf,'color','w')
ylabel('Light-Seq Sensitivity - Drop-Seq Sensitivity (%)')
ylim([-15 15])
xlim([0 17])
xticks([1:1:16])
xticklabels(t.Gene(sind))
title('Sensitivity Comparison Between Light-Seq and Drop-Seq vs. smFISH')


%% Create boxplot Light-Seq vs. Drop-Seq, colored by replicate
% Each point represents the estimated sensitivity of detecting one bipolar
% subtype marker gene gene in one replicate. This includes 16 genes x 4
% replicates for Light-Seq and 16 genes x 6 replicates for Drop-Seq

% Colored by replicates

figure;
subplot(1,2,1)
% Plot Light-Seq sensitivity
boxplot(Sensitivity(:),'Colors',[0 0 0],'whisker', 1.5);
for i = 1:4
hold on; scatter(ones(size(Sensitivity(i,:))),Sensitivity(i,:),100,'filled','MarkerFaceColor',[0 0.2 0.2]*i,'MarkerFaceAlpha',0.4,'jitter',0.1)
end
title('Sensitivity Compared to smFISH')
set(gca,'FontSize',16)
ylim([-0.1 26])
ylabel('Sensitivity (%)')
set(gcf,'color','w')
xticklabels('Light-Seq')

% Plot Drop-Seq sensitivity
subplot(1,2,2)
boxplot(Sensitivity_D(:),'Colors',[0 0 0],'whisker', 1.5);
for i = 1:6
hold on; scatter(ones(size(Sensitivity_D(i,:))),Sensitivity_D(i,:),100,'filled','MarkerFaceColor',[0 0.15 0.15]*i,'MarkerFaceAlpha',0.4,'jitter',0.1)
end

set(gca,'FontSize',16)
ylabel('Sensitivity (%)')
ylim([-0.1 26])
set(gcf,'color','w')
xticklabels('Drop-Seq')
title('Sensitivity Compared to smFISH')
%% Plot Colormap for boxplot replicates
figure;
hold on
for i=1:6
scatter(i,1,100,[0 0.15 0.15]*i,'filled','MarkerFaceAlpha',0.4)
end

%% Plot gene lengths on line
figure;
plot([1,2],[0,0])
yticks([length(sind)])
yticklabels(t.Gene(sind))
ylim([min(length),max(length)])