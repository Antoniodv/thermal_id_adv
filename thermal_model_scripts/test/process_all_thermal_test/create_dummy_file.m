function create_dummy_file(full_path)
    % Extract directory from path
    folder_path = fileparts(full_path);

    % Create dir if needed
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end

    % Create file and write dummy content
    fid = fopen(full_path, 'w');
    if fid == -1
        error(['Could not create file: ', full_path]);
    end
    fwrite(fid, 'dummy');
    fclose(fid);
end
