%% Script to graph IC50 Growth Curves of the LastDay
filenames = ["2022-06-28_BW_Nissle_Ic50_LastDay.xlsx","2022-06-23_F18_Lastday_anc.xlsx"]; 
day =[ 7, 7];
strain_name = [ "BW_Nissle", "F18"]; %days as strings to use for heatmap 
dilutionFN = [  "Dilutions FollowUp_IC50_BW_Nissle.xlsx","Dilutions FollowUp_IC50_F18.xlsx"];
% Concentration = 0 has been replaced with Concentration = 1E-10 in order
% for it to work in log 
%Put in index of data location in the excel 
% row start, row end, column start, column end 
rinds = [71 180; 184 293; 297 406; 410 519]; 
cinds = [2 99; 4 99; 4 99; 4 99 ]; %Skipping time and temperature columns of last 3 tables 
timeinterval = 10/60; %Time interval in hr
tp = [8 10]; %Timepoints to make heatmap  
rows = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"];
columns = 1:24; 
% Wells to remove due to errors when making the plate 
%wells_remove = {["G13", "G16", "E15", "F15", "H18"], ["A11", "M20", "N20", "O20"]};
% Serkan> Added more wells to remove due to weird high OD compared to other
% wells. Maybe a bubble in the well? Reading error?
wells_remove = {["G13", "G16", "E15", "F15", "H18","H9","K10","G13","H13","E15","F15","G16","H18"], ["A11", "M20", "N20", "O20","F6","H9"]};
%% Importing in data as tables and storing it into a cell
% Commented out: Graphing growth curves of 96 well plate 
all_data = {};
for i = 1:length(filenames)
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
    
    % Importing Dilution Scheme & Colony location from Dilution excel 
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
            ylim([0 1.2])
            hold on
            x=sprintf("%s: Growth Curves of All Wells (X= Time(HR), Y=Absorbance at OD600)", strain_name(i));
            sgtitle(x)
            ftitle = sprintf("%s: Day_%i_GrowthCurve.jpg", strain_name(i));
            % Making table of dilution and colony per well
            wellnames = [wellnames, wellname];
            dilutions = [dilutions, dil_table{k, rows(j)}];
            colonies = [colonies, colony_table{k, rows(j)}];
        
        end
    end
    hold off 
    %saveas(fig,ftitle)
    reference_table = table(wellnames', colonies', dilutions');
    reference_table.Properties.VariableNames =["Wellname",  "ColonyType","Dilution_uM"];
    % Removing wells 
    rwells = wells_remove{i};
    for i1 = 1:length(rwells)
        data_table{:,rwells(i1)} = nan;
    end 
    % Storing data 
    all_data{i,1} = day(i); %Day
    all_data{i,2} = data_table; % Data organized into table with Wellname, Time, OD 
    all_data{i,4} = reference_table; % Wellname, Colony Type, Dilution 
end 

%% Normalizing the data by subtracting OD at time = 0 
% Organizing the Data into Growth table of ColonyType, Dilution, Normalized
% Data, RawData, OD at Timepoints 
% Taking OD600 values at timepoints 
% Storing in all_data
timepoints = []; 

for q = 1:length(filenames)
    plate_data = all_data{q,2}; 
    reference_table = all_data{q,4}; 
    timepoints = [];
    strain = [];
    strain_type = [];
    dilution =[];
    normdata = {};
    max_growth = [];
    rawdata={};
    for ic=1:length(columns)
        for jr=1:length(rows)
            wellname = sprintf("%s%i",rows(jr), columns(ic));
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
            cs = reference_table.ColonyType(rind);
            strain = [strain,  cs ];
                if contains(cs, 'BW') == true
                    st = "BW";
                elseif contains(cs, 'Nissle') == true
                    st = "Nissle";
                elseif contains(cs, 'F-18') == true
                    st = "F-18";
                end
            strain_type = [strain_type, st];
            dilution =[dilution, reference_table.Dilution_uM(rind)];
            normdata = [normdata, [time, ngrowth]] ;
            rawdata = [rawdata, growth]; 
            timepoints = [timepoints; tps];
        end
    end
   

Growth = table(strain','VariableNames', {'ColonyType'});
Growth.Strain = strain_type'; 
Growth.Dilution = dilution'; 
Growth.Normalized_Data = normdata';
Growth.RawData = rawdata'; 
Growth.TimePoints = num2cell(timepoints);
all_data{q,3} = Growth; 
end 

% %% Plotting growth curves on top of each other
% mydat = all_data{1,3};
% ct = mydat.ColonyType;
% rd = mydat.RawData;
% dil = mydat.Dilution;
% udil = unique(dil);
% udil(isnan(udil))=[];
% uct = unique(ct);
% for iu = 1:length(uct)
%     uind = find(strcmp(uct(iu), ct));
%     urd =rd(uind);
%     cdil = dil(uind);
%     figure
%     for p = 1:7
%         subplot(1,7,p)
%         hold on 
%         indd = find(cdil==udil(p));
%         plot(data_table.Time, urd{indd(1)}, 'LineWidth',2)
%         plot(data_table.Time, urd{indd(2)}, 'LineWidth',2)
%         plot(data_table.Time, urd{indd(3)}, 'LineWidth',2)
%         legend(["Rep1", "Rep2", "Rep3"], "Location", "best")
%         title(udil(p))
%     end 
%     sgtitle(uct(iu))
% end 

%% Fitting a curve to the growth at timepoint 
% Fits a sigmoid curve and finds the IC50 using the function: doseResponse_sigmoid
% If all OD are greater than 0.5 (normalized OD) then a line is fitted to
% find the IC50
% If the fitted line has a positive slope (higher drug conc= more growth then no drug), instead of negative slope as it
% should be then IC50 is assumed to be greater then highest drug conc
% tested, so assumed 10^5 
for q1 = 1:length(filenames)
    ic50_cell = {};
    % looping over each day's data 
    Growth= all_data{q1,3};
    dilution = Growth.Dilution;
    strain = Growth.ColonyType;
    unidilution = unique(dilution);
    unidilution = unidilution(~isnan(unidilution));
    colonytypes = unique(Growth.ColonyType);
    colonytypes(strcmp(colonytypes,'PBS'))=[];
    strain_type = Growth.Strain; 
    for t =1:length(colonytypes)
        ic50_vector = zeros(1,length(tp));
        graph_curves = {}; 
        for tp1 = 1:length(tp)
            max_growth = cell2mat(Growth.TimePoints(:,tp1));
            ind = find(strcmp(strain,colonytypes{t}));
            dil = dilution(ind);
            s_name = unique(strain_type(ind));
            mg = max_growth(ind);
            mymat=[dil,mg];
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
            mat_response = [dose, response];
            indhalf = find(response<0.5); %Finding if response is greater than half of max no drug growth
            
            % If response at all dosages is greater then 1/2 of max no drug
            % growth then it is resistant at levels greater than highest
            % drug concentration(6000)
            if isempty(indhalf)
                p = polyfit(dose,response, 1);
                ic50 = (0.5-p(2))/p(1);
                xpoints= logspace(log10(10^-1), log10(max(dose)),5000);
                ypoints = polyval(p, xpoints);
                doses = unique(dose);
                for i=1:length(doses)
                    responses=response(dose==doses(i));
                    meanResponse(i)=mean(responses);
                    %stdErrResponse(i)=std(responses)/sqrt(length(responses));
                    stdErrResponse(i)=std(responses);
                end
                %Graphing is supressed 
                %                 figure
                %                 hold on
                %                 semilogx(xpred, ypred)
                %                 plot(dose, response, "b*")
                %                 plot(ic50, 0.5,"rO")
                %                 ylim([0,1])
                if ic50<0 %Correcting if IC50 value is a negative number
                    ic50 = 10^5;
                end
            else
                %Graphing of IC50 curve is suppresed in the function 
                [hillCoeff, ic50, mymat, coeffs, meanResponse, stdErrResponse,doses, xpoints,ypoints]=doseResponse_sigmoid_LastDay(dose,response,[50]);
            end
            ic50_vector(tp1) = ic50;
            graph_curves{tp1} = {doses,meanResponse, stdErrResponse, xpoints,ypoints};
        end
        %     ftitle = sprintf("%s_IC50.jpg", colonytypes{t});
        %     saveas(fig,ftitle)
        % t = colony type; q1= Day
        ic50_cell{t,1} = ic50_vector; % Storing the Ic50 values 
        ic50_cell{t,2} = graph_curves(1); % Storing growth curve for timepoint 1
        ic50_cell{t,3} = graph_curves(2); % Storing growth curve for timepoint 2
        ic50_cell{t,4} = colonytypes(t); % Storing colony type
        ic50_cell{t,5} = s_name; %Storing strain name 
    end
    all_data{q1,4} = ic50_cell; 
end 

%% Graphing the IC50 Growth Curve of Each Colony
% graph_curves = {doses,meanResponse, stdErrResponse, xpoints,ypoints};
% Combining IC50 data 
ic50_total = [all_data{1,4}; all_data{2,4}];
legtitle = ic50_total(:,4);
strains= cellstr(ic50_total(:,5));
%Serkan>Modified strains structure to exclude the unnecessary Ancestors and
%control lines
strains(2:3,1)={'NaN'}
strains(5:6,1)={'NaN'}
strains(11:12,1)={'NaN'}
strains(14:15,1)={'NaN'}
strains(20:21,1)={'NaN'}
strains(23:24,1)={'NaN'}
colony = ic50_total(:,4);
strain_type = ["BW", "Nissle", "F-18"];
%Serkan> Added shapes and custom colors
shapes={'o','v','^','>','d'}
cols=[[35 31 32];[35 31 32];[200 120 178];[200 120 178];[200 120 178]]
cols=cols./256
for h1 = 1:length(strain_type)
    sind = find(strcmp(strain_type(h1), strains)); %Getting index for day 7 data for strain
    lt = [];
    for il = 1:length(sind)
        st = erase(legtitle{sind(il)},strain_type(h1));
        lt =[lt,st, st];
    end
    D= char(lt);
    for tp2 = 1:length(tp)
        f=figure('color','w');
        c_map = jet(9);
        leg = [];
        for f1 = 1:length(sind)
            cind = sind(f1);
            growth_curve = ic50_total{cind,3}{:,:};
            doses = growth_curve{1,1};
            doses(doses==min(doses)) = 0.1; %Changing 0 to 0.1 for logscale
            meanResponse = growth_curve{1,2};
            stdErrResponse = growth_curve{1,3};
            xpoints= growth_curve{1,4};
            ypoints = growth_curve{1,5};
            semilogx(xpoints, ypoints, 'Color', cols(f1,:), 'LineWidth',2)
            hold on
            %errorbar(doses,meanResponse,stdErrResponse,shapes{1,f1},'Color',c_map(f1,:),'LineWidth',2,'MarkerSize',15,'MarkerFaceColor',c_map(f1,:),'MarkerEdgeColor','k')
            errorbar(doses,meanResponse,stdErrResponse,shapes{1,f1},'Color',cols(f1,:),'LineWidth',4,'MarkerSize',50,'MarkerFaceColor',cols(f1,:),'MarkerEdgeColor','k')
            set(gcf,'position',[[-196,1249,1146,555]]) %hard coded figure size
            xlim([0.01 10000])
            ylim([0 1.2])
        end
        legend(char(D), "Location", 'eastoutside')
        title(sprintf("%s Evolution Gemcitabine IC50 Curves, T-%i hrs", strain_type(h1), tp(tp2)))
        savefig(f,sprintf("%s_IC50_Curves_%i_hrs.fig", strain_type(h1), tp(tp2)));
        saveas(f,sprintf("%s_IC50_Curves_%i_hrs", strain_type(h1), tp(tp2)),'jpg');
        hold off
        clear f
    end
end
%% Heatmap of IC50 
daydata = ic50_total(:,1); %Getting ic50 data
strainlines = cellstr(ic50_total(:,5)); %getting strain
lines= ic50_total(:,4); %getting line 
sortedStrainLines = cellstr(sortrows(unique(cellstr(strainlines))));
    for i = 1:length(sortedStrainLines)
        fig1 = figure('color','w');
        lind = find(strcmp(sortedStrainLines(i),strainlines));%Getting index for strain
        ic50day = cell2mat(daydata(lind,:)); %Finding IC50 values corresponding to the strain
        slines = lines(lind); %Getting Lines corresponding to the strain
        h=heatmap(tp(1:2), slines, log10(ic50day),'CellLabelColor','none', 'ColorLimits', [0.5,5]);
        title(['IC50 Last Day & Ancestor Line:' sortedStrainLines(i)]);
        h.YLabel = 'Line';
        h.XLabel = 'Time (HR)';
        h.Colormap = flipud(redblue);
        stitle = sprintf('LastDay_Ancestor_IC50_Heatmap_%s.jpg', sortedStrainLines{i});
        saveas(fig1,stitle)
    end 






