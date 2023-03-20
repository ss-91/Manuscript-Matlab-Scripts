%% Definitions and Loading Whole Data
doses = [0.001 2 3.2 5.12 8.20 13.11 20.98 33.57 53.71 85.94 137.50 220];
yconditions={'No Bacteria','2M','4M','8M','16M','32M','64M','128M'};
n_strains = 13;
files = dir;
MetaFile = readtable('MetaData.xlsx');
Data ={};
IC50_response=[];

%% Organize the Data
for i=1:13
    Data{i,1} = MetaFile.Strains{i};
    Data{i,2} = readmatrix(files(i+69).name,'Range',MetaFile.T0_Range{i});
    Data{i,3} = readmatrix(files(i+69).name,'Range',MetaFile.Tend_Range{i});
end

%% Perform the Analysis
figure; hold on;
for i=1:n_strains
    if i==3 || i==2
       [ IC50_response(:,i) area(:,i)]= aSpheroid_forcdd(Data{i,2}, Data{i,3}, doses, Data{i,1}, yconditions)
    else
       [ IC50_response(:,i) area(:,i)]= aSpheroid_forrest(Data{i,2}, Data{i,3}, doses, Data{i,1}, yconditions)
    end
end

%% Quick IC50 Overview
figure; hold on;
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,1),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,2),'LineWidth',3,'Color','#f70202')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,3),'--xg','LineWidth',3)
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,4),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,5),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,6),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,7),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,8),'LineWidth',3,'Color','#0076a8')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,9),'--xk','LineWidth',3)
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,10),'LineWidth',3,'Color','#f70202')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,11),'LineWidth',3,'Color','#f70202')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,12),'LineWidth',3,'Color','#f70202')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,13),'LineWidth',3,'Color','#f70202')
plot([1 2 3 4 5 6 7 8 ], IC50_response(:,14),'LineWidth',3,'Color','#f70202')
legend(Data(:,1)')
grid on;
xticklabels({'No Bacteria','2M','4M','8M','16M','32M','64M','128M'});
ylabel('IC50 (uM)');
%%
figure;
IC50_response_Norm = log2(IC50_response ./ IC50_response(:,9))
imagesc(IC50_response_Norm,[-2 2])
xticks([1:1:19])
xticklabels(Data(:,1)')
xtickangle(45);
yticklabels({'No Bacteria','2M','4M','8M','16M','32M','64M','128M'});
c=redblue(100)
colormap(c)
hold off;
figure;
heatmap(IC50_response_Norm)

%% Plot the IC 50 Responses fro No Bacteria Condition
figure;
bar(IC50_response(1,:))
xticks([1:18])
xticklabels(Data(:,1)')
ylabel('Gemcitabine(uM)')
ylim([0 20]);
xtickangle(45)

%% Plot the Area under the IC50 Curve
strains = categorical({'cytR','yohK','rfaG','envZ','phoR','glnG','wild type','ihfB','yfjG','nupC','ubiF','ompR','cdd'});
strains = reordercats(strains,cellstr(strains)');
sorted_area= [159473	219820	240353	250751	254021	254785	262455	280632	289688	292240	296343	298429	414727];
figure;
bar(strains,sorted_area(1,:))

%% Plot the Area under the IC50 Curve as Fold Change
wt_area = 262455;
sorted_area_fc= sorted_area ./ wt_area;
sorted_area_fc= log2(sorted_area_fc);
figure;
imagesc(sorted_area_fc,[-1 1])
colormap(redgreencmap)

figure;
bar(sorted_area_fc)
grid on;
