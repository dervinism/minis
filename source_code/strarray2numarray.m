function numarray = strarray2numarray(strarray)
% strarray2numarray returns a row vector of numbers converted from a column
% vector of strings. The elements of the string vector have to be
% separated. For this purpose use Matlab "char" function instead of
% "horzcat" or "vertcat".
%

numarray = zeros(1,length(strarray));
for num = 1:length(numarray)
    numarray(num) = str2double(strarray(num,:));
end