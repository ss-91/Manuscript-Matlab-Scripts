%% Rapid Gemcitabine Breakdown Assay Analysis, Serkan Sayin, 2020-03-06 Massachusetts
%% Read in the Data, Labels and Definitions
%Definitions
DataFileName1 = '2021-09-11-Rapid-Gem-15min-Rep1-10ul-manual.xlsx';
DataFileName2 = '2021-09-11-Rapid-Gem-15min-Rep2-10ul-manual.xlsx';
DataFileName3 = '2021-09-11-Rapid-Gem-15min-Rep3-10ul-manual.xlsx';
DataRange = 'D29:CU149';
LabelFileName = 'Rapid_Gemcitabine_Degradation_Strains.xlsx';
LabelFileRange = 'B2:M9';
warning('off','all')

% Read in OD equalization data for all replicates
eq_od={}    
eq_od1=xlsread('2021-09-08-equal-ODs-Rep1.xlsx','C23:N30');
eq_od1(7:8,12)=NaN;
eq_od{1,1}=eq_od1;
eq_od2=xlsread('2021-09-09-equal-ODs-Rep2-90ul.xlsx','C23:N30');
eq_od2(7:8,12)=NaN;
eq_od{1,2}=eq_od2;
eq_od3=xlsread('2021-09-10-equal-ODs-Rep3.xlsx','C23:N30');
eq_od3(7:8,12)=NaN;
eq_od{1,3}=eq_od3;

% Read in OD Data for all replicates
ODdata1 = xlsread(DataFileName1, DataRange);
ODdata2 = xlsread(DataFileName2, DataRange);
ODdata3 = xlsread(DataFileName3, DataRange);
ODdata_combined = horzcat(ODdata1, ODdata2, ODdata3);
[x, y,Labels] = xlsread(LabelFileName,'Sheet2', LabelFileRange);
clear x y
dataFormat = 1;           % 1 if goes by row by column or 0 if goes by column by row
kineticCycle = 20.2;        % x hours
timeInterval = 10 ;       % x minutes

% Parameters to Set
select_tp=43; % This is for analysis up to 10 hours
pval = 0.05; %p-val cutoff

%%%%Plate Format Info%%%%
rowNumber = 8;
colNumber = 12;
rawOD8x12 = cell(rowNumber,colNumber);
normOD8x12 = cell(rowNumber,colNumber);
%%%%%%%%%%%%%%%%%%%%%%%%
time = [0:10/60:((kineticCycle*60)-timeInterval)/60];
%Background Subtraction
OD0 = median(ODdata_combined([1:3],:),1); % treat the first 3 time points as blanks
normOD = ODdata_combined-OD0;
[r,c]=find(normOD<0)
normOD(r,c)=0;
n = size(ODdata_combined,1);
%% Quality Control of OD Equalization
labels={'Rep1','Rep2','Rep3'}
figure; 
for i=1:3;
subplot(1,3,i)
histogram(eq_od{1,i},20)
xlim([0 0.2])
ylim([0 50]);
ylabel('Absorbance');
grid on;
title(labels{1,i})
end

figure;
for i=1:3
subplot(1,3,i)
imagesc(eq_od{1,i},[0 0.18])
xticks([1:12]);
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12'});
yticklabels({'A','B','C','D','E','F','G','H'})
title(labels{1,i})
end
%% Organize the Data Format into 96/384-well Cell Array and Perform Background Subtraction
if dataFormat == 0; 
    i = 1;
    for c=1:colNumber;
        for r=1:rowNumber;
            curOD = ODdata1(:,i);
            curNormOD = normOD(:,i);
            normOD8x12{r,c} = curNormOD;
            rawOD8x12{r,c} = curOD;
            i = i+1;
        end
    end
else if dataFormat == 1;
        i=1;
        for r=1:rowNumber;
            for c=1:colNumber;
                curOD = horzcat(ODdata_combined(:,i), ODdata_combined(:,i+96), ODdata_combined(:,i+192));
                curNormOD = horzcat(normOD(:,i), normOD(:,i+96), normOD(:,i+192));
                rawOD8x12{r,c} = curOD;
                normOD8x12{r,c} = curNormOD;
                i = i+1;
            end
        end
        
    end
end

%% Gather all wt data in one array
norm_wt=[];
j=1;
for i=1:6
     temp=normOD8x12{i,12}
     temp2=rawOD8x12{i,12}
     norm_wt(:,j:j+2)=temp;
     raw_wt(:,j:j+2)=temp2;
     j=j+3
end

mean_normwt=mean(norm_wt,2)
mean_rawwt=mean(raw_wt,2)
%% plot the raw and normalized OD graphs
H1=figure('color','w','Position', get(0, 'Screensize')); 
i = 1;
maxOD = max(max(rawOD8x12{8,12}(1:select_tp,:)));
minOD = 0;
colors={'g','r','b'};

for r=1:rowNumber;
    for c=1:colNumber;
        H1=subplot(rowNumber,colNumber,i);
        for j=1:3;
        hold on; grid minor;
        plot(time(1:select_tp),rawOD8x12{r,c}(1:select_tp,j),colors{1,j},'Linewidth',2);
        set(gca,'ylim',[0 maxOD],'xlim',[min(time) time(select_tp)]); box on;
        xticks([0 time(select_tp)])
        yticks([0 maxOD])
        end
        plot(time(1:select_tp),mean_rawwt(1:select_tp,1),'--k','Linewidth',2)
        i = i+1;
        title(Labels(r,c))
    end
end
suptitle({'del-cdd raw growth curves (supernatant 15min)','green=rep1, red=rep2, blue=rep3, dashed=wt'});


H2=figure('color','w','Position', get(0, 'Screensize')); 
i = 1;
maxOD = max(max(normOD8x12{8,12}(1:select_tp,:)));
minOD = 0;
colors={'k','r','b'};

for r=1:rowNumber;
    for c=1:colNumber;
        H2=subplot(rowNumber,colNumber,i);
        for j=1:3;
        hold on; grid minor;
        plot(time(1:select_tp),normOD8x12{r,c}(1:select_tp,j),colors{1,j},'Linewidth',2);
        set(gca,'ylim',[0 maxOD],'xlim',[min(time) time(select_tp)]); box on;
        xticks([0 time(select_tp)])
        yticks([0 maxOD])
        end
        plot(time(1:select_tp),mean_normwt(1:select_tp,1),'--k','Linewidth',2)
        i = i+1;
        title(Labels(r,c))
    end
end
suptitle({'del-cdd normalized growth curves (supernatant 15min)','green=rep1, red=rep2, blue=rep3, dashed=wt'});


%% Remove data from bad replicates due to unequal OD while performing the experiment (from NormOD)
%normOD8x12{3,1}(:,1)=NaN %remove Rep1 C1
%normOD8x12{1,7}(:,1)=NaN %remove Rep1 A7
%normOD8x12{1,4}(:,2)=NaN %remove Rep2 A4
%normOD8x12{7,11}(:,3)=NaN %remove Rep3 G11
%% Remove well because one replicate behaves differently
%normOD8x12{2,4}=1;
%normOD8x12{3,7}=1;

%% Calculate the Mean and Std for the Replicates
for r=1:rowNumber;
    for c=1:colNumber;
        if normOD8x12{r,c}~=1
        normOD8x12{r,c}(:,4) = nanmean(normOD8x12{r,c}(:,1:3),2);
        normOD8x12{r,c}(:,5) = nanstd(normOD8x12{r,c}(:,1:3),[],2);
        end
    end
end

%% Calculate the Area Under the curve in 7 hours of growth for each replicate seperately
AUC_Array = nan(96,1); %AUC_8x12 = nan(rowNumber,colNumber);
i = 1;
x = time(1,1:select_tp)'; %only up to 7 hrs
for r=1:rowNumber;
    for c=1:colNumber
        if normOD8x12{r,c}~=1  
        y1 = normOD8x12{r,c}(:,1);
        y2 = normOD8x12{r,c}(:,2);
        y3 = normOD8x12{r,c}(:,3);
            AUC1 = trapz(x, y1(1:select_tp));
            AUC2 = trapz(x, y2(1:select_tp));
            AUC3 = trapz(x, y3(1:select_tp));
            AUC=[AUC1 AUC2 AUC3];
            AUC_8x12{r,c} = AUC;
            i = i+1;
        else
            i=i+1;
        end
    end
end

% Group wt AUC data together (total of 18 values)
wt_AUC_array=[];
for i=1:6
    wt_AUC_array=[wt_AUC_array,AUC_8x12{i,12}]
end
wt_AUC_array=wt_AUC_array'

%% Perform one tailed Student's t-test
%Fast Degraders
for r=1:8
    for c=1:12
      [h1(r,c) stat_fast_degraders(r,c)]=ttest2(AUC_8x12{r,c},wt_AUC_array,'tail','right');
    end
end

p_cutoff=0.1
stat_vector = stat_fast_degraders(:)
Labels_vector = Labels (:)
stat_vector_adjusted = mafdr(stat_vector)
hits= Labels_vector(stat_vector_adjusted<p_cutoff)
hits_padj=stat_vector(stat_vector_adjusted<p_cutoff)


%% Color the Slow and Fast Gemcitabine Degraders according to statistics
% Generate Index Matrix
inx_matrix =zeros(8,12);
for i=1:length(hits)
        tmp_matrix=strcmp(Labels, hits(i,1));
        inx_matrix = inx_matrix+tmp_matrix;
end

%% Mark the fast degraders
H4= figure('Position', get(0, 'Screensize'));
i = 1;
maxOD = max(max(normOD8x12{8,12}(1:select_tp,:)));
minOD = min(min(ODdata1));

for r=1:rowNumber;
    for c=1:colNumber;
      if normOD8x12{r,c}~=1
        H4= subplot(rowNumber,colNumber,i); hold on; box on;
        if inx_matrix(r,c)==1
        plot(time,normOD8x12{r,c}(:,4),'-g','LineWidth',2);
        plot(time,mean_normwt,'--k','Linewidth',2);
 
        %xticks([0 :2:time(select_tp)]);
        xticks([0 2 4 6 8])
        xticklabels([])
        yticks([0 :0.1:maxOD]);
        grid on;
        i = i+1;
        set(gca,'ylim',[0 maxOD],'xlim',[min(time) time(select_tp)]);
        else
        plot(time,normOD8x12{r,c}(:,4),'-r','LineWidth',2);
        plot(time,mean_normwt,'--k','Linewidth',2);
        xticks([0 2 4 6 8])
        xticklabels([])
        %num2str(round(tao8x12(r,c)))
        %xticks([0 :2:time(select_tp)]);
        yticks([0 :0.1:maxOD]);
        grid on;
        i = i+1;
        set(gca,'ylim',[0 maxOD],'xlim',[min(time) time(select_tp)]);       
        end
      else
      H4= subplot(rowNumber,colNumber,i); hold on;
      plot(0,0)
      i=i+1;
      end
      title(Labels{r,c});
    end
end
suptitle({'Significant Fast Degraders(FDR p-adj<' pval,'AUC-' time(select_tp) 'hrs'});

%%











