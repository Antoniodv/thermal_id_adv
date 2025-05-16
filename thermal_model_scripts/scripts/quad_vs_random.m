% functions to move in the grid
function neighbors = get_neighbor_ids(core_id, n_rows, n_cols)
 
    [x_core, y_core] = id_to_coordinates(core_id, n_rows, n_cols);

    [x_north, y_north] = find_north_coordinate(x_core, y_core, n_rows, n_cols);
    id_north = coordinates_to_id(x_north, y_north, n_rows, n_cols);

    [x_south, y_south] = find_south_coordinate(x_core, y_core, n_rows, n_cols);
    id_south = coordinates_to_id(x_south, y_south, n_rows, n_cols);

    [x_east, y_east] = find_east_coordinate(x_core, y_core, n_rows, n_cols);
    id_east = coordinates_to_id(x_east, y_east, n_rows, n_cols);

    [x_west, y_west] = find_west_coordinate(x_core, y_core, n_rows, n_cols);
    id_west = coordinates_to_id(x_west, y_west, n_rows, n_cols);

    neighbors = [ id_north, id_south, id_east, id_west ]; 
    neighbors = neighbors(neighbors ~= -1); 
end

function [x, y] = id_to_coordinates(core_id, n_rows, n_cols)
    % this function returns x, y coordinates on a grid given core id
    % x is the horizontal coordinate
    % y is the vertical coordinate
    % numeration of slots (core_id) has the lower number at the beginning of
    % each row on the left and increases from top to bottom

    x = mod(core_id - 1, n_cols) + 1;
    y = floor((core_id - 1) / n_cols) + 1;
end

function core_id = coordinates_to_id(x_core, y_core, n_rows, n_cols)
    % this function returns core_id on a grid given its coordinates
    % x is the horizontal coordinate
    % y is the vertical coordinate

    if (x_core ~= -1) && (y_core~= -1)
        y_contrib = (y_core - 1) * n_cols;
        x_contrib = x_core;
    
        core_id = x_contrib + y_contrib;    
    else
        core_id = -1;
    end
end

function [x, y] = find_north_coordinate(x_core, y_core, ~, ~)
    if y_core > 1
        x = x_core;
        y = y_core - 1;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_south_coordinate(x_core, y_core, n_rows, ~)
    if y_core < n_rows
        x = x_core;
        y = y_core + 1;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_east_coordinate(x_core, y_core, ~, n_cols)
    if x_core < n_cols
        x = x_core + 1;
        y = y_core;
    else
        x = -1;
        y = -1;
    end
end

function [x, y] = find_west_coordinate(x_core, y_core, ~, ~)
    if x_core > 1
        x = x_core - 1;
        y = y_core;
    else
        x = -1;
        y = -1;
    end
end


