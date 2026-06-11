% compileMinisStandalone.m
% Compile the minis GUIDE app as a standalone Windows application and
% package it as a web installer (minisInstaller_web.exe).
%
% Requirements: MATLAB Compiler toolbox must be installed.
%
% Output: minisStandalone\for_redistribution\minisInstaller_web.exe

%% 1. Define paths
sourceFolder  = fileparts(mfilename('fullpath'));
if isempty(sourceFolder)
    sourceFolder = pwd;
end
repoFolder    = fileparts(sourceFolder);           % …\minis
outputFolder  = fullfile(repoFolder, 'minisStandalone');

appFile       = fullfile(sourceFolder, 'minis.m');
splashFile    = fullfile(sourceFolder, 'splashScreen.jpg');

fprintf('Source folder : %s\n', sourceFolder);
fprintf('Repo folder   : %s\n', repoFolder);
fprintf('Output folder : %s\n', outputFolder);

%% 2. Collect all additional files in the source_code folder
contents = dir(sourceFolder);
additionalFiles = {};
for iFile = 1:numel(contents)
    fname = contents(iFile).name;
    % Skip hidden / current-dir entries and the main app file itself
    if fname(1) == '.' || strcmpi(fname, 'minis.m')
        continue
    end
    additionalFiles{end+1, 1} = fullfile(sourceFolder, fname); %#ok<AGROW>
end

fprintf('Number of additional files: %d\n', numel(additionalFiles));

%% 3. Build standalone application
fprintf('\n--- Building standalone application ---\n');
opts = compiler.build.StandaloneApplicationOptions(appFile, ...
    'ExecutableName',        'minis', ...
    'AdditionalFiles',       additionalFiles, ...
    'OutputDir',             outputFolder, ...
    'ExecutableIcon',        splashFile, ...
    'ExecutableSplashScreen',splashFile, ...
    'ExecutableVersion',     '1.0', ...
    'Verbose',               'On');

buildResults = compiler.build.standaloneApplication(opts);
fprintf('Build complete.\n');

%% 4. Package as a web installer
fprintf('\n--- Packaging as web installer ---\n');
installerFile = fullfile(outputFolder, 'minisInstaller_web.exe');

compiler.package.installer(buildResults, ...
    'InstallerName',         'minisInstaller_web', ...
    'OutputDir',             outputFolder, ...
    'RuntimeDelivery',       'web', ...
    'Verbose',               'On');

%% 5. Verify output
if ispc
    candidate = fullfile(outputFolder, 'for_redistribution', 'minisInstaller_web.exe');
    if exist(candidate, 'file')
        fprintf('\nSuccess! Installer created at:\n  %s\n', candidate);
        % Copy to the top level of minisStandalone as requested
        copyfile(candidate, fullfile(outputFolder, 'minisInstaller_web.exe'), 'f');
        fprintf('Also copied to:\n  %s\n', fullfile(outputFolder, 'minisInstaller_web.exe'));
    else
        % Some MATLAB versions put it directly in OutputDir
        if exist(installerFile, 'file')
            fprintf('\nSuccess! Installer created at:\n  %s\n', installerFile);
        else
            warning('Could not find minisInstaller_web.exe – check the outputFolder for the generated files.');
            disp(dir(outputFolder));
        end
    end
end
