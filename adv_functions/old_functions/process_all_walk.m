function find_and_process_files(root_dir)
    % Avvia il motore MATLAB
    fprintf("Avvio motore MATLAB...\n");

    % Ottieni tutte le sottocartelle di root_dir
    allFolders = strsplit(genpath(root_dir), ';'); % Lista di tutte le sottocartelle

    for i = 1:length(allFolders)
        currentFolder = allFolders{i};

        if contains(currentFolder, filesep + "power_model") % Cerca solo cartelle power_model
            % Costruisci il percorso del file gb_core_uncore_tot_temp.pkl
            pkl_file = fullfile(currentFolder, "gb_core_uncore_tot_temp.pkl");

            if exist(pkl_file, 'file')
                % Determina la cartella di output
                parent_dir = fileparts(currentFolder); % Cartella che contiene power_model
                output_folder = fullfile(parent_dir, "thermal_model");

                % Assicurati che la cartella di output esista
                if ~exist(output_folder, 'dir')
                    mkdir(output_folder);
                end

                fprintf("Processing: %s\n", pkl_file);
                fprintf("Output Folder: %s\n", output_folder);

                % Chiamata alla funzione MATLAB thermal_sa
                thermal_sa(pkl_file, output_folder);
            end
        end
    end

    fprintf("Processo completato.\n");
end
