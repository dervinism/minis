function [fileName, filePath, fileList] = findFile(attribute, order)
% Find the newest to the oldest file or the first to the last file
% alphabetically.
%
% The function finds the short name and the full name (including a path)of
% a file of a specified saving date order or alphabetical order in a given
% directory path. It also lists all the files in a directory sorted
% according to this particular attribute starting with the newest one or
% the first in the alphabetical order.
%
% The input variables are string 'attribute' (either 'dateSaved' or 'name')
% and order (from the newest to the oldest/
% from the first to the last alphabetically: 1,2,...,'end').
%
% If the specified folder contains no files, the function returns a 1x3
% array of zeros.


D = dir; % Getting the folder structure
dirIndex = [D.isdir]; % Index array of directories
noFiles = numel(dirIndex)- sum(dirIndex); % Number of files in the directory
if ~noFiles
    fileName = 0; filePath = 0; fileList = 0;
    return
end
if strcmpi(order, 'end')
    order = noFiles;
end

%% Finding the short and the full names of a file:
if strcmpi(attribute, 'dateSaved')
    filesDate = zeros(noFiles, 2); % The list of file order and date saved
elseif strcmpi(attribute, 'name')
    filesName = cell(noFiles, 1); % The list of file names
end
iFile = 0; % The row order of variable 'Files'.
for iDir = 3:length(D) % The element of variable 'D'
    if ~D(iDir).isdir
        iFile = iFile + 1;
        if strcmpi(attribute, 'dateSaved')
            filesDate(iFile,1) = iDir; % Indexing files
            filesDate(iFile,2) = D(iDir).datenum; % Extracting the date when the file was saved
        elseif strcmpi(attribute, 'name')
            filesName{iFile} = D(iDir).name; % Extracting the file name
        end
    end
end
format long

if strcmpi(attribute, 'dateSaved')
    sortedFiles = flipud(sortrows(filesDate,2)); % Sorting files according to the date they were saved.
    fileName = (D(sortedFiles(order,1)).name); % Finding the file of a particular saving date order
    filePath = fullfile(pwd, fileName); % Finding the path of the file
elseif strcmpi(attribute, 'name')
    sortedFiles = sort(filesName); % Sorting files alphabetically.
    fileName = char(sortedFiles(order)); % Finding the file of a particular alphabetical order
    filePath = fullfile(pwd, fileName); % Finding the path of the file
end

%% Finding the list of the files sorted according to the date saved:
if strcmpi(attribute, 'dateSaved')
    fileList = cell(noFiles, 1); % The list of files sorted according to the date saved
    for jFile = 1:size(sortedFiles,1) % The row order of variable 'sortedFiles'
        fileList{jFile,1} = D(sortedFiles(jFile,1)).name; % Extracting the name of the file
    end
elseif strcmpi(attribute, 'name')
    fileList = sortedFiles;
end