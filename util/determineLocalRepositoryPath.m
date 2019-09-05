function LOCALPATH = determineLocalRepositoryPath
% This function determines the system-local path to the StereoPIV project
% repository folder.

% Determine the path to the current directory
currentDirectory = pwd;

% Determine the location in the current-directory string of the word
% 'StereoPIV', which we know is common to this project. Add 8 because this
% function finds the location of the letter 'S' in 'StereoPIV', and we want to
% find the location of the end of the word 'StereoPIV'
repositoryLoc = regexpi(currentDirectory, 'syntheticImageGeneration') + 23;

% This is the system-local path to the StereoPIV folder.
LOCALPATH = currentDirectory(1:repositoryLoc(1));

end