
function [NU_t] = count_NU(t, map_sheet, NU_0)

% map sheet = [tkl, RTk, wck, Tck]

NU_t = NU_0;
U = 0;
size_map_sheet = size(map_sheet);
for k = 1 : size_map_sheet(1)
    if t < map_sheet(k,1)
        U = 0;
    elseif t > map_sheet(k,1) && map_sheet(k,3) == 0
        U = 1;
        continue;
    elseif t > map_sheet(k,1) && map_sheet(k,3) ~= 0
        U = 2;
    end
    NU_t = NU_t + (ceil((t-map_sheet(k,1))/(2*pi/map_sheet(k,3)))*U*map_sheet(k,2));
end
end