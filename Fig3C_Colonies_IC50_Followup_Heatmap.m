%% Script to Create Heatmap of IC50 results from colonies from each strain 
filenames = ["2022-05-11_Plate1_followup.xlsx", "2022-05-12_Plate2_followup.xlsx","2022-06-01_Gemcitabine_Followup_Ic50_2.xlsx" ]; 
plate =[1, 2, 3];
dayvar = ["Plate 1", "Plate 2", "Plate 3"]; %days as strings to use for heatmap 
dilutionFN = ["Dilutions FollowUp_IC50._plate1.xlsx", "Dilutions FollowUp_IC50_plate2.xlsx", "Dilutions FollowUp_IC50_tak2.xlsx"];
% Concentration = 0 has been replaced with Concentration = 1E-10 in order
% for it to work in log 
%Put in index of data location in the excel 
% row start, row end, column start, column end 
rinds = [71 180; 184 293; 297 406; 410 519]; 
cinds = [2 99; 4 99; 4 99; 4 99 ]; %Skipping time and temperature columns of last 3 tables 
timeinterval = 10/60; %Time interval in hr
tp = [8]; %Timepoints to make heatmap  
rows = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"];
columns = 1:24; 
% Wells to remove due to errors when making the plate 
wells_remove = {["H19", "G19", "P12"], ["E15"], ["A8","C8","E8","G8"]};

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
    %fig = figure('color','w','Position', get(0, 'Screensize')); 
    ip=1;
    wellnames = [];
    dilutions = [];
    colonies = [];
    lines = [];
    for j = 1:length(rows)
        for k = 1:length(columns)
            wellname = sprintf("%s%i",rows(j), columns(k));
            %Graphing of 384 well growth curve supressed 
%             subplot(16, 24,ip);
%             plot(data_table.Time, data_table{:,wellname}, 'LineWidth',2);
%             ip = ip+1;
%             ylim([0 1.1])
%             hold on
%             x=sprintf("Plate %i : Growth Curves of All Wells (X= Time(HR), Y=Absorbance at OD600)", plate(i));
%             sgtitle(x)
%             ftitle = sprintf("Plate_%i_GrowthCurve.jpg", plate(i));
            % Making table of dilution and colony per well
            wellnames = [wellnames, wellname];
            dilutions = [dilutions, dil_table{k, rows(j)}];
            colonies = [colonies, colony_table{k, rows(j)}];
            mycolony = cell2mat(colony_table{k, rows(j)});
            k1 = strfind(mycolony, " Colony");
            lines = [ lines, {mycolony(1:k1)}];
        end
    end
    %hold off 
    %saveas(fig,ftitle)
    reference_table = table(wellnames', lines', colonies', dilutions');
    reference_table.Properties.VariableNames =["Wellname",  "Line", "ColonyType", "Dilution_uM"];
    % Removing wells
    rwells = wells_remove{i};
    for i1 = 1:length(rwells)
        data_table{:,rwells(i1)} = nan;
    end 
    % Storing data 
    all_data{i,1} = plate(i);
    all_data{i,2} = data_table;
    all_data{i,5} = reference_table; 
end 


%% Normalizing the data by subtracting OD at time = 0 
% Taking OD600 values at timepoints 
% Storing in all_data
timepoints = []; 

for q = 1:length(filenames)
    plate_data = all_data{q,2}; 
    reference_table = all_data{q,5}; 
    timepoints = [];
    strain = [];
    dilution =[];
    normdata = {};
    max_growth = [];
    line=[];
    rawdata={};
    for i=1:length(columns)
        for j=1:length(rows)
            wellname = sprintf("%s%i",rows(j), columns(i));
            growth=plate_data{:,wellname};
            ngrowth = growth-growth(1); %Normalizing the data by subtracting OD at time=0
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
            line = [line, reference_table.Line(rind)];
            normdata = [normdata, [time, ngrowth]] ;
            rawdata = [rawdata, growth]; 
            timepoints = [timepoints; tps];
        end
    end
   

Growth = table(line','VariableNames', {'Line'});
Growth.ColonyType = strain';
Growth.Dilution = dilution'; 
Growth.Normalized_Data = normdata';
Growth.RawData = rawdata'; 
Growth.TimePoints = num2cell(timepoints);
all_data{q,3} = Growth; 
end 


%% Fitting a curve to the growth at timepoint to determine IC50
ic50_cell = {}; 
strainnames = [];
strainlines = [];
for q1 = 1:length(filenames)
    Growth= all_data{q1,3};
    dilution = Growth.Dilution;
    strain = Growth.ColonyType;
    unidilution = unique(dilution);
    unidilution = unidilution(~isnan(unidilution));
    colonytypes = unique(Growth.ColonyType);
    colonytypes(strcmp(colonytypes,'PBS'))=[];
    lines = Growth.Line;
    %%
    for t = 1:length(colonytypes)
        myline=[];
        ic50_vector = zeros(1,length(tp));
        for tp1 = 1:length(tp)
            max_growth = cell2mat(Growth.TimePoints(:,tp1));
            ind = find(strcmp(strain,colonytypes{t}));
            dil = dilution(ind);
            mg = max_growth(ind);
            mymat=[dil,mg];
            myline = unique(lines(ind));
            % remove Nan
            indN = find(isnan(mg));
            mg(indN) = [];
            dil(indN) = [];
            indL = find(dil== min(dil)); %finding index of lowest drug conc

            %compute mean control response
            controlResponse=mean(mg(indL),'omitnan');
            %remove controls from dose/response curve
            response=mg/controlResponse;
            dose=dil;
            indhalf = find(response<0.5); %Finding if response is greater than half of max no drug growth
            % Response at all dosages is greater then 1/2 of max no drug
            % growth then it is resistant at levels greater than highest
            % drug concentration(6000) so using a best fit line instead of
            % sigmoid to determine IC50
            if isempty(indhalf)
                p = polyfit(dose,response, 1);
                ic50 = (0.5-p(2))/p(1);
                xpred = logspace(log10(0.1), log10(max(dose)),5000);
                ypred = polyval(p, xpred);
%                                  figure
%                                  hold on
%                                  semilogx(xpred, ypred)
%                                  plot(dose, response, "b*")
%                                  plot(ic50, 0.5,"rO")
%                                  ylim([0,1])
%                                  title (colonytypes{t})
%                                  hold off 
                if ic50<0 %Correcting if IC50 value is a negative number
                    ic50 = 6000; % Assuming IC50 is 10^5 because it is resistent at all drug concentrations (growth>0.5) and so IC50 must be greater than largest drug concentration
                    %>Serkan, corrected this value. Bec. the max drug dose
                    %in the experiment is 6mM.
                end
            else
                % Using a Sigmoid to determine IC50
                % mymat= [dose, response]
                % 50 = IC50
                [hillCoeff, ic50, mymat, coeffs]=doseResponse_sigmoid(dose,response,[50]);
            end
            ic50_vector(tp1) = ic50;
        end
        %     ftitle = sprintf("%s_IC50.jpg", colonytypes{t});
        %     saveas(fig,ftitle)
        % t = colony type; q1= Day
        ic50_cell{t,q1} = ic50_vector;
        strainlines = [strainlines; myline];
    end
    strainnames = [strainnames; colonytypes];

end 


%% Making IC50 Heatmap 
x= find(cellfun(@isempty, ic50_cell)); %Finding index of where there is no IC50 data because plate 3 has less data
ic50_cell(x) = {[nan, nan]};
daydata = [];
%>Serkan added custom colormap
mycolormap = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
colormap(mycolormap);
ic50day={};
for t=1:length(tp)
    oneday = cellfun(@(v)v(t),ic50_cell);
    daydata = [daydata, [oneday(:,1);oneday(:,2); oneday(:,3)]];
end 

    sortedStrainLines = sortrows(unique(strainlines));
    for i = 1:length(sortedStrainLines)
        lind = find(strcmp( sortedStrainLines(i),strainlines));%Getting index for Line
        %fig1 = figure('color','w');
        if sortedStrainLines{i} == "BW Line 3 " 
            % Using the data from 6/1 for BW Line 3 as we choose the colony for WGS for that
            subcolonies = 1:8;
            ic50day(:,i) = num2cell(daydata(lind(lind>96),:)); %Finding OD values corresponding to the line from Plate 3 
        elseif sortedStrainLines{i} == "F-18 Line 3 "
            % Using the data from first two plates for F-18 Line 3 as we choose the colony for WGS for that
            subcolonies = 1:8;
            ic50day(:,i) = num2cell(daydata(lind(lind<96),:)); %Finding OD values corresponding to the line from Plate 1 & 2
        elseif sortedStrainLines{i} == "F-18 Line 1 "
            % Using the data from first two plates for F-18 Line 3 as we choose the colony for WGS for that
            subcolonies = 1:8;
            ic50day(:,i) = num2cell(daydata(lind(lind<96),:)); %Finding OD values corresponding to the line from Plate 1 & 2
        else 
            subcolonies = 1:length(lind);
            ic50day(:,i) = num2cell(daydata(lind,:));
        end 
            
        %h=heatmap(tp(1:2), subcolonies, ic50day ./ max(ic50day) ,'CellLabelColor','none', 'ColorLimits', [0,1]);
        %title(sortedStrainLines(i));
        %h.Colormap = flipud(redblue);
        %h.Colormap = mycolormap;
        %h.YLabel = 'Colony Number';
        %h.XLabel = 'Time (HR)';
        %title(sprintf("%s log10(IC50) Colonies", sortedStrainLines{i}))
        %saveas(fig1,sprintf("%s_IC50_Heatmap_Colonies.jpg", sortedStrainLines{i}))
        %savefig(fig1,sprintf("%s_IC50_Heatmap_Colonies", sortedStrainLines{i}))
    end 

    %Normalization of the IC50 Data(to the max value of each strain's most
    %resistant)
ic50day=cell2mat(ic50day)
%ic50day(:,1:4) = ic50day(:,1:4) ./ max(max(ic50day(:,1:4)))
%ic50day(:,5:8) = ic50day(:,5:8) ./ max(max(ic50day(:,5:8)))
%ic50day(:,9:12) = ic50day(:,9:12) ./ max(max(ic50day(:,9:12)))

ic50day=log10(ic50day)
%ic50day(:,1:4) = log2(ic50day(:,1:4) ./ median(ic50day(:,1)))
%ic50day(:,5:8) = log2(ic50day(:,5:8) ./ median(ic50day(:,5)))
%ic50day(:,9:12) = log2(ic50day(:,9:12) ./ median(ic50day(:,9)))


%% Making IC50 Heatmap as subplots
%% Making pie charts as subplots
x= find(cellfun(@isempty, ic50_cell)); %Finding index of where there is no IC50 data because plate 3 has less data
ic50_cell(x) = {[nan, nan]};
daydata = [];
%>Serkan added custom colormap
%mycolormap = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
mycolormap = customcolormap(linspace(0,1,10), {'#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});

colormap(mycolormap);

explode=[0,0,0,0,1,0,0,0;
         0,0,0,0,0,1,0,0;
         0,0,0,0,0,1,0,0;
         1,0,0,0,0,0,0,0;
         0,0,0,0,0,0,1,0;
         0,0,0,0,0,1,0,0;
         0,0,0,0,0,1,0,0;
         0,0,0,0,0,1,0,0;
         0,0,0,0,0,0,1,0;
         0,1,0,0,0,0,0,0;
         0,1,0,0,0,0,0,0;
         0,0,0,1,0,0,0,0]


     
for t=1:length(tp)
    oneday = log10(cellfun(@(v)v(t),ic50_cell));
    daydata = [oneday(:,1);oneday(:,2); oneday(:,3)];
fig1 = figure('color','w','Position', get(0, 'Screensize'));
    sortedStrainLines = sortrows(unique(strainlines));
    for i = 1:length(sortedStrainLines)
        subplot(3,4,i)
        lind = find(strcmp( sortedStrainLines(i),strainlines));%Getting index for Line
        h=heatmap(tp(t), subcolonies, ic50day(:,i),'ColorLimits', [0,3.778]);
        %h=heatmap(tp(t), subcolonies, ic50day(:,i),'CellLabelColor','none', 'ColorLimits', [0,15]);
        %h=pie([12.5 12.5 12.5 12.5 12.5 12.5 12.5 12.5], explode(i,:))
        title(sortedStrainLines(i));
        %h.Colormap = flipud(gray);
        h.Colormap = mycolormap;
        h.YLabel = 'Colony Number';
        h.XLabel = 'Time (HR)';
        title(sortedStrainLines{i})
    end 
    sgtitle(sprintf("log10(IC50) for Individual Colonies (Normalized to Each Strain's Max) T=%i hrs", tp(t)))
    %saveas(fig1,sprintf('IndividualColonies_IC50_Heatmap_%i_hrs.jpg', tp(t)))
    %savefig(fig1, sprintf('IndividualColonies_IC50_Heatmap_%i_hrs', tp(t)))
end 
























