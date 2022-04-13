% function used in previous cours to obtain coordinate in the correct range
function [lonshift, datashift] = MapData_centerGreenwich(lon, data)
lon = lon(:);
shiftl = floor(length(lon)/2);
datashift = circshift(data, [0, shiftl]);
lonshift = circshift(wrapTo180(lon), shiftl);
lonshift(shiftl) = [];
lonshift = [-lonshift(end); lonshift];
datashift(:, shiftl) = [];
datashift = [datashift(:,end), datashift];
end

