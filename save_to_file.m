function r = save_to_file(filename, data, mode)
    if nargin < 3
        mode = 'a';  
    end
    fileID = fopen(filename, mode);
    fprintf(fileID, [data '\n']);
    fclose(fileID);
end