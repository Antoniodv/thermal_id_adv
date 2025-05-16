function row_values = get_row_from_file(row_index, file_path)

    % open file
    fid = fopen(file_path, 'r');
    if fid == -1
        error('File opening error');
    end
    
    % read row by row
    numeric_rows = {};  % Cell array for numeric rows
    while ~feof(fid)
        line = fgetl(fid);
        
        % Check if row is only numeric
        if ~isempty(line) && all(ismember(line, '0123456789 '))
            % convert to numbers and save
            numeric_rows{end+1} = sscanf(line, '%d')';
        end
    end
    
    % close file
    fclose(fid);
    
    % Checks on input index
    if row_index < 1 || row_index > length(numeric_rows)
        error('Index out of range. The file contains %d numeric rows', length(numeric_rows));
    end
    
    % Return the selected row
    row_values = numeric_rows{row_index};
end
