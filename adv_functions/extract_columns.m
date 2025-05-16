function subMatrix = extract_columns(matrix, indices)
    % Estrarre tutte le righe e solo le colonne specificate da 'indices'
    subMatrix = matrix(:, indices);
end