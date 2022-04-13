function [] = plotGS(GS_cell)
    
% Color = [0 0 0]/1;

for k = 1:length(GS_cell)
    GS = GS_cell{k};
    plot(GS.lon,GS.lat,'o','LineWidth', 2,'MarkerSize', 5, 'MarkerEdgeColor', [1 1 1]) 
    text(GS.lon, GS.lat+7*sign(GS.lat),GS.name,'FontSize', 14, 'Color', [1 1 1], 'LineWidth',2,'FontWeight','bold', 'HorizontalAlignment', 'right', 'Interpreter', 'tex')
end

end

