function data = parseKatzSolutionFiles(folderPath,type)
% Initialize a structure to hold the extracted data
data = struct();

% Get a list of all files in the specified folder
files = dir(fullfile(folderPath, type));

% Loop through each file
for i = 1:length(files)
    if ~files(i).isdir
        % Open the file for reading
        fileID = fopen(fullfile(folderPath, files(i).name), 'r');
        if fileID == -1
            warning('Could not open file: %s', files(i).name);
            continue;
        end

        % Initialize variables for the table
        matrixValue = [];
        timeToSolve = [];
        alphaValue = [];

        % Read the file line by line
        while ~feof(fileID)
            line = fgetl(fileID);

            % Extract relevant data using regular expressions
            if contains(line, 'Using PSBLAS version:')
                data(i).psblasVersion = extractAfter(line, 'Using PSBLAS version: ');
            elseif contains(line, 'Using AMG4PSBLAS version:')
                data(i).amg4psblasVersion = extractAfter(line, 'Using AMG4PSBLAS version: ');
            elseif contains(line, 'Bandwidth and profile:')
                data(i).bandwidthProfile = sscanf(line, 'Bandwidth and profile: %f %f %f');
            elseif contains(line, 'Preconditioner:')
                data(i).preconditioner = extractAfter(line, 'Preconditioner: ');
            elseif contains(line, 'Iterations to convergence:')
                data(i).iterationsToConvergence = sscanf(line, 'Iterations to convergence: %d');
            elseif contains(line, 'Time to solve system               :')
                timeToSolve = sscanf(line, 'Time to solve system               : %e');
                data(i).timetosolve = timeToSolve;
            elseif contains(line, 'Residual 2-norm:')
                data(i).residual2Norm = sscanf(line, 'Residual 2-norm: %f');
            elseif contains(line, 'Matrix:')
                matrixValue = extractAfter(line, 'Matrix: ');
                % Parse the matrix field to extract the string between
                % '../testmatrices/mm/' and '.mtx'
                matrixValue = extractBetween(matrixValue, '../testmatrices/mm/', '.mtx');
                if ~isempty(matrixValue) 
                    matrixValue = matrixValue{1};
                    % Sanitize the string substituting any "_" with "\_"
                    matrixValue = strrep(matrixValue, '_', '\_');
                end
            elseif contains(line, 'Value of alpha                     :')
                alphaValue = sscanf(line, 'Value of alpha                     : %f');
            end
        end

        % Store the extracted values in the data structure
        data(i).matrix = matrixValue;
        data(i).alpha = alphaValue;

        % Close the file
        fclose(fileID);
    end
end