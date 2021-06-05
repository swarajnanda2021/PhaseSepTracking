function [ data ,siz ] = fload( file )
%get file data

i=find(file=='.',1,'last');
fileID=fopen([file(1:i-1),'.bin'],'r'); % write and discard existing content
siz=fread(fileID, [1 2],'double');
fclose(fileID); % write and discard existing content

fileID=fopen(file,'r'); % write and discard existing content
data=fread(fileID, siz,'double');
fclose(fileID); % write and discard existing content

end

