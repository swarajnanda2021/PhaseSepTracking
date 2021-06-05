function [ mem ] = fmemmap( file )
%create_memmap create memory map using a the bin and dat file
%   input, file name incl dir root. file is open for writing.

i=find(file=='.',1,'last');
fileID=fopen([file(1:i-1),'.bin'],'r'); % write and discard existing content
siz=fread(fileID, [1 2],'double');
fclose(fileID); % write and discard existing content

mem=memmapfile(file,...
    'Format',{'double',siz,'dat'},...
    'Writable',true);

end