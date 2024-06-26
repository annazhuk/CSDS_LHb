function [Sample,scope,camera] = importfile_ImagingDAQ(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [SAMPLE,SCOPE,CAMERA] = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   [SAMPLE,SCOPE,CAMERA] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads
%   data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [Sample,scope,camera] = importfile('MatlabImportTool_paste_2590417401237086283tmp.txt',9, 134842);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2022/03/03 15:17:59

%% Initialize variables.
delimiter = {'\t',','};
if nargin<=2
    startRow = 9;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column3: double (%f)
%   column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%*s%f%f%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Sample = dataArray{:, 1};
scope = dataArray{:, 2};
camera = dataArray{:, 3};


