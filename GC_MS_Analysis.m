%% Analyze GC-MS Data
time1 =[0 20 75];
time2 =[0 20 45 75];
evo_gem = xlsread('GC-MS.xlsx','GEM','E2:G46');
evo_dfdU= xlsread('GC-MS.xlsx','dFdU','E2:G46');
hits_gem = xlsread('GC-MS.xlsx','GEM','E48:H95');
hits_dfdu = xlsread('GC-MS.xlsx','dFdU','E48:H95');
controls_gem = xlsread('GC-MS.xlsx','GEM','E97:H108');
controls_dfdu = xlsread('GC-MS.xlsx','dFdU','E97:H108');

[x,y,BW_labels] = xlsread('GC-MS.xlsx','GEM','B2:B16');
labels={'BW Evo1','BW Evo2','BW Evo3','BW Control','BW Ancestor';
        'EcN Evo1','EcN Evo2','EcN Evo3','EcN Control','EcN Ancestor';
        'F18 Evo1','F18 Evo2','F18 Evo3','F18 Control','F18 Ancestor'};
labels2 = {'ompR','glnG','envZ','yfjG','yohK','cytR','ihfB','phoR','ubiF','bioA','rfaG','cdd','nupC','BW25113','rfaG(geno)','rfaH(geno)'};
labels3={'Gemcitabine','Only BW25113','Only EcN','Only F18'};
stat_gem=[];
stat_dfdu=[];

%% Convert Peak Areas to Concentrations
evo_gem = getconc_gem(evo_gem);
evo_dfdU = getconc_dfdu(evo_dfdU);
hits_gem = getconc_gem(hits_gem);
hits_dfdu = getconc_dfdu(hits_dfdu);
controls_gem = getconc_gem(controls_gem);
controls_dfdu=getconc_dfdu(controls_dfdu);

%% Plot the Gemcitabine Breakdown for Evolved Strains
k=1;
l=13;
for j=1:3
    figure;
for i=1:5
    subplot(1,5,i)
    hold on; grid on; box on;
    tmp1(:,1:3) = evo_gem(k:(k+2),:);
    tmp2(:,1:3) = evo_dfdU(k:(k+2),:);
    errorbar(time1,mean(tmp1),std(tmp1),'-o','Linewidth',1,'MarkerSize',5); %gem plot
    %errorbar(time1,mean(tmp2),std(tmp2),'-o','Linewidth',1,'MarkerSize',5); %dfdu plot
    tmp3_wt(:,1:3)=evo_gem(l:(l+2),:); %gem_wt
    errorbar(time1,mean(tmp3_wt),std(tmp3_wt),'--ok','Linewidth',1,'MarkerSize',5);%
    %dfdu wt
    %tmp3(:,1:3)=evo_dfdU(l:(l+2),:);
    %errorbar(time1,mean(tmp3),std(tmp3),'--ok','Linewidth',1,'MarkerSize',5)
    k=k+3;
    ylim([-10 120]);
    xlim([0 80]);
    title(labels{j,i});
    xticks([0 20 75]);
    xlabel('Time (mins)');
    ylabel('Gemcitabine Concentration (uM)');
    %ylabel('dFdU Concentration (uM)');
    set(gca,'FontSize',15);
    set(gcf,'Position',[31 522 1440 120]);
end
legend('GEM','wt');
%legend('dFdU','wt');
l=l+15;
%legend('GEM');
end


%% Plot Gemcitabine Breakdown for ASKA Hits Strains
k=1;
figure;
for i=1:16
    subplot(4,4,i); hold on; grid on; box on;
    tmp1(:,1:4) = hits_gem(k:(k+2),:);
    tmp2(:,1:4) = hits_dfdu(k:(k+2),:);
    errorbar(time2,mean(tmp1),std(tmp1),'-o','Linewidth',2,'MarkerSize',5); %gem
    errorbar(time2,mean(tmp2),std(tmp2),'-o','Linewidth',2,'MarkerSize',5);
    %dfdu
    tmp3_wt(:,1:4) = hits_gem(40:42,:);
    errorbar(time2,mean(tmp3_wt),std(tmp3_wt),'--ok','Linewidth',1,'MarkerSize',5);
    tmp4_wt(:,1:4) = hits_dfdu(40:42,:);
    errorbar(time2,mean(tmp4_wt),std(tmp4_wt),'--ok','Linewidth',1,'MarkerSize',5);  
    k=k+3;
    ylim([-10 100]);
    xlim([0 80]);
    title(labels2{1,i})
    xticks([0 20 45 75]);
    xlabel('Time (mins)');
    ylabel('gemcitabine (uM)');
    set(gca,'FontSize',10);
    set(gcf, 'Position', [238 356 1171 600]);
         
        j=2
        for j=2:4
        [h1(i,j-1) stat_gem_fast(i,j-1)]=ttest2(tmp1(:,j),tmp3_wt(:,j),'tail','left');
        [h2(i,j-1) stat_dfdu_fast(i,j-1)]=ttest2(tmp2(:,j), tmp4_wt(:,j),'tail','right');
         
        [h3(i,j-1) stat_gem_slow(i,j-1)]=ttest2(tmp1(:,j),tmp3_wt(:,j),'tail','right');
        [h4(i,j-1) stat_dfdu_slow(i,j-1)]=ttest2(tmp2(:,j), tmp4_wt(:,j),'tail','left');
        end
end
legend('GEM','dFdU','WT Gem');
labels2=labels2'

%% FDR Correct the p-values
for i=1:3
stat_gem_fast_cor(:,i)=mafdr(stat_gem_fast(:,i),'BHFDR',1)
stat_dfdu_fast_cor(:,i)=mafdr(stat_dfdu_fast(:,i),'BHFDR',1)
end

for i=1:3
stat_gem_slow_cor(:,i)=mafdr(stat_gem_slow(:,i),'BHFDR',1)
stat_dfdu_slow_cor(:,i)=mafdr(stat_dfdu_slow(:,i),'BHFDR',1)
end


%% Plot Gemcitabine Breakdown for Controls
k=1;
figure;
for i=1:4
    subplot(1,4,i); hold on; grid on;
    tmp1(:,1:4) = controls_gem(k:(k+2),:);
    tmp2(:,1:4) = controls_dfdu(k:(k+2),:);
    errorbar(time2,median(tmp1),std(tmp1),'-o','Linewidth',2,'MarkerSize',5);
    errorbar(time2,median(tmp2),std(tmp2),'-o','Linewidth',2,'MarkerSize',5);
    k=k+3;
    ylim([0 120]);
    xlim([0 80]);
    title(labels3{1,i});
    xticks([0 20 45 75]);
    xlabel('Time (mins)');
    ylabel('Peak Area');
    set(gca,'FontSize',10);
end
legend('GEM','dFdU');
