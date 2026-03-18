function CreateRadarEchogram(CS,FIELD)

%Creates an radar echogram

defval('FIELD','data') %'data'/'coherence'/'phase_difference'
defval('FntSz',14)     %Set fontsize to be used in figures

if isequal(CS.GEO.OPERATION_MODE,'SIR_L1B_SARIN')
    WD = reshape(CS.SIN.(FIELD),size(CS.SIN.data,1),size(CS.SIN.data,2)*size(CS.SIN.data,3));
elseif isequal(CS.GEO.OPERATION_MODE,'SIR_L1B_SAR')
    WD = reshape(CS.SAR.data,size(CS.SAR.data,1),size(CS.SAR.data,2)*size(CS.SAR.data,3));
elseif isequal(CS.GEO.OPERATION_MODE,'SIR_L1B_LRM')
    WD = reshape(CS.LRM.data,size(CS.LRM.data,1),size(CS.LRM.data,2)*size(CS.LRM.data,3));
else
    error('Mode %s not reckognized!',CS.GEO.OPERATION_MODE)
end

if strcmp(FIELD,'data')
    WD = WD./repmat(max(WD),size(WD,1),1);
end

ScreenSize = get (0,'Screensize');
figure ('Position',ScreenSize);
imagesc(1:size(WD,2),1:size(WD,1),WD)
set(gca,'YDir','normal')
if strcmp(FIELD,'data')
    caxis([0 1])
else
    caxis(minmax(WD))
end
xlabel('Record','fontsize',FntSz)
ylabel('Bin','fontsize',FntSz)
set(gca,'FontSize',14)
colorbar
set(gcf, 'Color',[1 1 1])
if strcmp(FIELD,'data')
    set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String', 'Normalized power','FontSize',14);
else
    set(get(findobj(gcf,'tag','Colorbar'),'ylabel'),'String', FIELD,'FontSize',14);
end
% axis image

end