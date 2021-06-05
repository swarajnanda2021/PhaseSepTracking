function fmerge( file1, file2 )
%merge file2 to file1

% get data file1
[data2]=fload(file2);

% save to file2
fsave(file1,data2,'a')

end