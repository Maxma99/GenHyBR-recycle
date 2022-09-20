function overlay2figures(Zprime,Z)
% %create data
% [X,Y]=meshgrid(1:21,1:41);
% Z =zeros(41,21);
% Z(10,10) = 10;
% Zprime = fluo_g_disp(:,:,9);
[X,Y]=meshgrid(1:size(Zprime,2),1:size(Zprime,1));
contourmin = min(Z(:));
contourmax = max(Z(:));
pcolormin = min(Zprime(:));
pcolormax = max(Zprime(:));

%create figure and store handle
hF = figure;
colormap(hF,bone);
%create axes for pcolor and store handle
hAxesP = axes;
%set colormap for pcolor axes
colormap(hAxesP,bone);

%plot pcolor for gradient
pcolorPlot = pcolor(X,Y,Zprime);
% set(pcolorPlot,'ydir','normal');
set(pcolorPlot,'FaceColor','interp','EdgeColor','interp');
%set(h2,'color','none','visible','off');
%create color bar and set range for color
cbP = colorbar(hAxesP,'Location','west');
caxis(hAxesP,[pcolormin pcolormax]);

%create axes for the countourm axes
hAxesCM = axes;
%set visibility for axes to 'off' so it appears transparent
axis(hAxesCM,'off')
%set colormap for contourm axes
colormap(hAxesCM,bone);

%plot contourm
contourmPlot = contourm(Y,X,Z);
% pcolorPlot2 = pcolor(X,Y,Z);
% set(pcolorPlot2,'color','none','visible','off')
%create color bar and set range for color
% colormap(pcolorPlot2,'jet');
%     set(pcolorPlot2,'ydir','normal');
cbCM = colorbar(hAxesCM,'Location','east');
% set(h2,'ydir','normal');
caxis(hAxesCM,[contourmin contourmax]);

%link the two overlaying axes so they match at all times to remain accurate
linkaxes([hAxesP,hAxesCM]);
end
