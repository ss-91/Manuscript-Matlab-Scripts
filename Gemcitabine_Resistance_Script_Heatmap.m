%% Script to analyze results from IC50 portion of Gemcitabine Resistance

filenames = ["2022-05-02_Gemcitabine_Day0.xlsx", "2022-05-03_Gemcitabine_Day1.xlsx", "2022-05-04_Gemcitabine_Day2.xlsx", "2022-05-05_Gemcitabine_Day3.xlsx", "2022-05-06_Gemcitabine_Day4.xlsx", "2022-05-07_Gemcitabine_Day5.xlsx", "2022-05-08_Gemcitabine_Day6.xlsx", "2022-05-09_Gemcitabine_Day7.xlsx"]; 
day =[0, 1, 2, 3, 4, 5, 6, 7];
dayvar = ["Day 0", "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7"]; %days as strings to use for heatmap 
dilutionFN = ["Dilutions Day0_1.xlsx", "Dilutions Day0_1.xlsx", "Dilutions Day2_on.xlsx", "Dilutions Day2_on.xlsx", "Dilutions Day2_on.xlsx", "Dilutions Day5.xlsx", "Dilutions Day2_on.xlsx","Dilutions Day2_on.xlsx"];
% Concentration = 0 has been replaced with Concentration = 1E-10 in order
% for it to work in log 
%Put in index of data location in the excel 
% row start, row end, column start, column end 
rinds = [71 180; 184 293; 297 406; 410 519]; 
cinds = [2 99; 4 99; 4 99; 4 99 ]; %Skipping time and temperature columns of last 3 tables 
timeinterval = 10/60; %Time interval in hr
tp = [ 8, 10, 12, 14, 16, 18]; %Timepoints to make heatmap  
rows = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"];
columns = 1:24; 
% Wells to remove due to errors when making the plate 
wells_remove = {["G22", "A15", "C15", "E15", "G15", "I15", "K15", "M15", "O15"];...
    ["A13", "C13", "E13", "G13", "I13", "K13", "M13", "O13", "A19", "C19", "E19", "G19", "I19", "K19", "M19", "O19"];...
    ["B6", "D6", "F6", "H6", "J6", "L6", "N6", "P6", "I1","I2","I3","I4","I5","I6","I7","I8", "I9", "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", "I20", "I21", "K21", "M21"];...
    ["I21", "K21", "M21"];...
    ["I7", "B12", "D12", "F12", "H12", "J12", "L12", "N12", "P12"];...
    ["B12", "D12", "F12", "H12", "J12", "L12", "N12", "P12", "M18"];...
    ["H24"]; ["C20", "A23","A12", "C12", "E12", "G12", "I12", "K12", "M12", "O12"]};

numcolonies = 18; %Number of colonies
%% Importing in data as tables and storing it into a cell

all_data = {};
for i = 1: length(filenames)
    [num,txt,raw] = xlsread(filenames(i));
    data = [];

    for j = 1:length(rinds)
        dat = raw(rinds(j,1):rinds(j,2), cinds(j, 1): cinds(j,2));
        data = [data, dat];
    end

    %Removing temperature column
    data(:,2) = [];

    %Making table
    data_table = cell2table(data(2:end,:));
    data_table.Properties.VariableNames = data(1,:);
    r = height(data_table);
    h = 0:timeinterval:((r-1)*timeinterval);
    data_table.Properties.VariableNames{'Time'} = 'Time';
    data_table{:,"Time"} =h(1:end)'; % Replacing time values with 0: last time point in hours
    
    % Importing Dilution Scheme & Colony location
    [~,~,rawD] = xlsread(dilutionFN(i));
    y=rawD(2:18, 1:25)';
    dil_table = cell2table(y(2:end, :));%Gemcitabine dilution table
    dil_table.Properties.VariableNames = y(1,:);
    colony_table = cell2table(rawD(22:38, 2:25)'); %Colony Table
    colony_table.Properties.VariableNames = rawD(22:38, 1)';
    fig = figure('color','w','Position', get(0, 'Screensize')); 
    ip=1;
    wellnames = [];
    dilutions = [];
    colonies = [];
    for j = 1:length(rows)
        for k = 1:length(columns)
            wellname = sprintf("%s%i",rows(j), columns(k));
            subplot(16, 24,ip);
            plot(data_table.Time, data_table{:,wellname}, 'LineWidth',2);
            ip = ip+1;
            ylim([0 1.1])
            hold on
            x=sprintf("Day %i : Growth Curves of All Wells (X= Time(HR), Y=Absorbance at OD600)", day(i));
            sgtitle(x)
            ftitle = sprintf("Day_%i_GrowthCurve.jpg", day(i));
            % Making table of dilution and colony per well
            wellnames = [wellnames, wellname];
            dilutions = [dilutions, dil_table{k, rows(j)}];
            colonies = [colonies, colony_table{k, rows(j)}];
        
        end
    end
    hold off 
    saveas(fig,ftitle)
    reference_table = table(wellnames', colonies', dilutions');
    reference_table.Properties.VariableNames =["Wellname",  "ColonyType","Dilution_uM"];
    % Removing wells 
    rwells = wells_remove{i};
    for i1 = 1:length(rwells)
        data_table{:,rwells(i1)} = nan;
    end 
    % Storing data 
    all_data{i,1} = day(i);
    all_data{i,2} = data_table;
end 


%% Normalizing the data by subtracting OD at time = 0 
% Taking OD600 values at timepoints 
% Storing in all_data
timepoints = []; 

for q = 1:length(filenames)
    plate_data = all_data{q,2}; 
    timepoints = [];
    strain = [];
    dilution =[];
    normdata = {};
    max_growth = [];
    rawdata={};
    for i=1:length(columns)
        for j=1:length(rows)
            wellname = sprintf("%s%i",rows(j), columns(i));
            growth=plate_data{:,wellname};
            if q == 6
                ngrowth = growth-growth(9);%Normalizing the data by subtracting OD at timepoint 9 for file 6 because plate moved 
            else 
                ngrowth = growth-growth(1); %Normalizing the data by subtracting OD at time=0
            end 
            time=plate_data{:,"Time"};
           
            tps= []; 
            for h = 1:length(tp)
                timepoint= ngrowth(time==tp(h));
                tps = [tps, timepoint]; 
            end 

             % Storing data
            rind = find(reference_table.Wellname==wellname);
            strain = [strain, reference_table.ColonyType(rind) ];
            dilution =[dilution, reference_table.Dilution_uM(rind)];
            normdata = [normdata, [time, ngrowth]] ;
            rawdata = [rawdata, growth]; 
            timepoints = [timepoints; tps];
        end
    end
   

Growth = table(strain','VariableNames', {'ColonyType'});
Growth.Dilution = dilution'; 
Growth.Normalized_Data = normdata';
Growth.RawData = rawdata'; 
Growth.TimePoints = num2cell(timepoints);
all_data{q,3} = Growth; 
end 


%% Making a Heatmap of the OD at Highest Concentration 
for q1 = 1:length(filenames)
Growth = all_data{q1,3};
dilution = Growth.Dilution;
unidilution = unique(dilution); 
unidilution = unidilution(~isnan(unidilution));
colonytypes = unique(Growth.ColonyType); 
colonytypes(strcmp(colonytypes,'Nothing'))=[];

% Collecting the OD at Highest Concentration 
% Storing the data in all_data 
OD_data = {}; 
 for tp1 = 1:length(tp)
     highMG = []; 
     for t = 1:length(colonytypes)
         max_growth = cell2mat(Growth.TimePoints(:,tp1));
         datmg=[];
         ind = find(strcmp(strain,colonytypes{t}));
         dil = dilution(ind);
         mg = max_growth(ind);
         mymat=[dil,mg];
         for t1 = 1:length(unidilution)
             datmg(:,t1)=mg(dil==unidilution(t1))';
         end

         MeanGrowth = mean(datmg, 'omitnan'); %getting mean across column
         stdErrResponse = std(datmg,0, 'omitnan'); % getting standard error

         indL = find(unidilution== min(unidilution)); %finding index of lowest drug conc
         normMG = MeanGrowth./MeanGrowth(indL); %Normalizing by lowest drug conc

         indH = find(unidilution== max(unidilution)); %finding index of highest drug conc
         highMG = [highMG, normMG(indH)];
     end
    OD_data{tp1} = highMG;
end
all_data{q1, 4}=OD_data; 
end

%% Making Heatmap
fig1 = figure('color','w','Position', get(0, 'Screensize'));
inx=[1 2 3 6 7 8 9 12 13 14 15 18];
mycolormap = customcolormap(linspace(0,1,8), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d'});
colorbar('southoutside');
colormap(mycolormap);
axis off;
for t=1:length(tp)
    oneday = [];
    for d = 1:length(day)
        daydata = all_data{d,4};
        oneday = [oneday; daydata{t}];
    end
    subplot(2,length(tp)/2, t);
    oneday=oneday'
    h=heatmap(day, colonytypes(inx,1), oneday(inx,:), 'CellLabelColor','none', 'ColorLimits',[0 1]);
    title(sprintf("Time(hr)= %d", tp(t)));
    h.YLabel = 'Colony Type';
    h.XLabel = 'Day'; 
    h.Colormap = mycolormap;
end 
sgtitle('Normalized OD at Highest Concentration')

saveas(fig1,'HighestOD_Heatmap.svg')


















