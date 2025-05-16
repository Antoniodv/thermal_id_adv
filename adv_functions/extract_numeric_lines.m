function numericMatrix = extract_numeric_lines(filename)
    % Apri il file per la lettura
    fid = fopen(filename, 'r');
    if fid == -1
        error('Impossibile aprire il file.');
    end
    
    numericMatrix = [];
    
    % Leggi il file riga per riga
    while ~feof(fid)
        line = fgetl(fid);
        
        % Controlla se la riga è composta solo da numeri (interi o decimali)
        if is_valid_numeric_line(line)
            % Converti la stringa in un array numerico
            numericArray = str2num(line); %#ok<ST2NM>
            if ~isempty(numericArray)
                numericMatrix = [numericMatrix; numericArray]; %#ok<AGROW>
            end
        end
    end
    
    % Chiudi il file
    fclose(fid);
end

function isValid = is_valid_numeric_line(line)
    % Controlla se la linea è vuota o contiene caratteri non numerici
    isValid = ~isempty(line) && all(ismember(line, '0123456789+-. eE')); 
    % ismember controlla che siano solo numeri, segni, punti o spazi
end
