function [IC50_Response,area_pixels] = aSpheroid(Response_T0,Response_Tend, doses, strain, ycondition)
%This functions analyzes the spheroid IC50 from Spheroid-Bacteria-Drug
%Coincubation experiments, It requires estIC50 function in the directory
%(Serkan Sayin, 2020-12-28)
%   The following are the plots made by this function
%   1- Plot the heatmap of the response ( subplots, ResponseT0,
%   ResponseTend, Response_Normalized)
%   2-Calculate IC50 for every ycondition, plot and pass them to
%   IC50_Response output
%   3-Fit a surface to the data
%   4-Plot surface, real datapoints and IC50 line

%% Normalization of the Spheroid Size to T0
%Response_Norm = Response_Tend ./ Response_T0; %correcting for initial spheroid size differences
%Response_Norm2 = Response_Norm ./ max(max(Response_Norm));
Response_Norm=[];
Response_Norm2=[];
    Response_Norm = Response_Tend ./ Response_T0 %Normalization to Timepoint Zero
    sorted = Response_Norm(:,1);
    sorted = sort(sorted,'descend');
    sorted = sorted(~isnan(sorted));
    sorted = sorted(2,1);
    Response_Norm2 = ((Response_Norm - min(min(Response_Norm)))  ./ (sorted-min(min(Response_Norm))));
    for r=1:8
        for c=1:12
            if Response_Norm2(r,c) > 1
               Response_Norm2(r,c) = 1;
            end
        end
    end
%% Plotting Heat Maps
% Heatmap 1
figure;
H=subplot(3,1,1)
imagesc(Response_T0,[0 400000])
title(strain)
yticks([1 2 3 4 5 6 7 8]);
yticklabels(ycondition)
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
xticklabels(doses);
xtickangle(90);
ylabel('Bacteria Concentration')
xlabel('Drug Concentration (uM)');
title([strain ' Response at T0 (Area,micrometer sq)'])

%Heatmap 2
H=subplot(3,1,2)
imagesc(Response_Tend,[300000 600000])
title(strain)
yticks([1 2 3 4 5 6 7 8]);
yticklabels(ycondition)
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
xticklabels(doses);
xtickangle(90);
ylabel('Bacteria Concentration')
xlabel('Drug Concentration (uM)');
title([strain ' Response at T endpoint (Area,micrometer sq)'])

%Heatmap 3
H=subplot(3,1,3)
imagesc(Response_Norm2,[0 1])
title(strain)
yticks([1 2 3 4 5 6 7 8]);
yticklabels(ycondition)
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
xticklabels(doses);
xtickangle(90);
ylabel('Bacteria Concentration')
xlabel('Drug Concentration (uM)');
title([strain  ' Response Normalized (Tn / T0)'])
saveas(H, ['01_Heatmaps_' strain])

%% Fitting IC50 Curves
Fits_Response=cell(1,8)
figure;hold on;
for i=1:8
subplot(2,4,i)
idx = isnan(Response_Norm2(i,:)');
Fits_Response{i}= fit(doses(~idx)', Response_Norm2(i,~idx)','smoothingspline','SmoothingParam',0.01);
plot(Fits_Response{i},doses',Response_Norm2(i,:))
xlabel('Doses')
ylabel('Area D6 / Area D0')
ylim([0 1]);
title([strain ' ' ycondition{1,i}])
set(gca, 'XScale', 'log')
end

mat_Response = [];
for i = 1:8
      mat_Response(i,:) = Fits_Response{i}(doses)
      IC50_Response(i,1) = estIC50(Fits_Response{i},doses)     
end

figure;
s=mesh(mat_Response)
s.FaceColor = 'flat';
xlabel('Drug Doses');
ylabel('Bacteria Titration');
zlabel('Fitness');
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
xticklabels({'0','2','3.2','5.12','8.20','13.11','20.98','33.57','53.71','85.94','137.50','220'});
yticks([1 2 3 4 5 6 7 8]);
yticklabels(ycondition)
title([strain ' 3D IC50 Fit'] );
zlim([0 1])
saveas(s, ['02_IC50_3DPlot_' strain]);


%% Fit a surface to the data 
[X,Y]=ndgrid(1:12,1:8);
[X_hi,Y_hi]=ndgrid(1:0.01:12,1:0.01:8);
   
   % Plot the Mesh Surface
   z=mat_Response'
   sf=fit([X(:),Y(:)],z(:),'poly44');
   Z_hi = sf(X_hi,Y_hi);
   %Z_hi_norm = Z_hi ./ max(max(Z_hi))
   figure; hold on;
   q=surf(X_hi,Y_hi,Z_hi);
   caxis([0 1])
   %colormap(hot)
   colormap(redblue)
   shading interp
  
   % Plot the IC50 Line
   %mid_res =  max(max(Z_hi_norm))/2
   mid_res = 0.5;
   for i=1:size(Z_hi,1)
       min_val= abs(Z_hi(i,:)-mid_res);
       tmp=find(min_val==min(min_val));
       if min(min_val)>0.01
          inx_mat(i,1)=i;
          inx_mat(i,2)=NaN;  
       else
       inx_mat(i,1)=i;
       inx_mat(i,2)=tmp;   
       end
   end
   
   s=isnan(inx_mat)
   s=sum(s)
   inx_mat=inx_mat((s(1,2)+1):end,:);
   
   if any(isnan(inx_mat(:,2)))==1
    for i=1:size(Z_hi,1)
       min_val= abs(Z_hi(i,:)-mid_res);
       tmp=find(min_val==min(min_val));
       inx_mat(i,1)=i;
       inx_mat(i,2)=tmp;   
    end 
   end
   
   
   area_pixels = sum(inx_mat(find(inx_mat(:,2)>100),2))
   
   norm_z= sf(X_hi(inx_mat(:,1))',Y_hi(1,inx_mat(:,2)));
   plot3(X_hi(inx_mat(:,1)),Y_hi(1,inx_mat(:,2)),norm_z,'-k','LineWidth',2);
   %xlabel('Drug Doses(uM)');
   %ylabel('Bacteria Titration');
   %zlabel('Response');
   xticks([1:12]);
   yticks([1:8]);
   xticklabels(doses)
   yticklabels(ycondition)
   %title([strain ' 3D Poly Surface Fit'] );
   
   
   %Plot the Actual Data
   %c=colormap(hot(96));
   c=redblue(96)
   Response_Norm2 = Response_Norm2'
   %[B I] = sort(tmp5,'ascend')
   %c=c(I,:);
   zz=zeros(12,8)
   zz= zz + 1.5
   i=1;
   hold on;
   for x=1:12
       for y=1:8
           scaleFactor = 1;
           %if Response_Norm2(x,y) < 0.2
            %  Response_Norm2(x,y) = 0.2
            %  tmp=1
           %else
           tmp=round((Response_Norm2(x,y)-(1-scaleFactor))/scaleFactor*96)
           if tmp==0
               tmp=1
           end
           if  isnan(tmp)==1
           plot3(X(x,y),Y(x,y),zz(x,y),'x','MarkerEdgeColor','#D3D3D3','MarkerSize',10,'Linewidth',4)
           elseif Response_Norm2(x,y)<0.5
           curColor = c(tmp,:);
           plot3(X(x,y),Y(x,y),zz(x,y),'o','MarkerFaceColor',curColor,'LineWidth',2,'MarkerEdgeColor','#ffffff','MarkerSize',Response_Norm2(x,y)*20+5)
           else
           curColor = c(tmp,:);
           plot3(X(x,y),Y(x,y),zz(x,y),'o','MarkerFaceColor',curColor,'LineWidth',2,'MarkerEdgeColor','#ffffff','MarkerSize',Response_Norm2(x,y)*20+5)  
           i=i+1;
           end
       end
   end
   axis off
   
   saveas(q, ['03_PolySurf_3DPlot_' strain]);  
end
