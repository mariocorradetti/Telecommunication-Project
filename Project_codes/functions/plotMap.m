function [] = plotMap(lon, lat, DATA, alphaData)

myCMap = jet;

if nargin < 4
    imagesc(lon, lat, DATA)
elseif nargin == 4
    imagesc(lon, lat, DATA, 'AlphaData', alphaData) 
end
set(gca,'XTick',[-180:30:180],'XTickMode','manual');
set(gca,'YTick',[-90:30:90],'YTickMode','manual');
set(gca,'YDir','normal')
colormap(gca, myCMap); shading interp;
axis equal
axis tight
grid on
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
hold on

load coastlines 
geoshow(coastlat, coastlon, 'Color', 'k', 'LineWidth', 1.0)

end

