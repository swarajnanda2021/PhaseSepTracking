function fsave( file, data, perm )
%get file data

% get size
siz=size(data);
if strcmp(perm,'a') && exist(file,'file')~=0
    [~,siz_old]=fload(file);
    I=siz~=siz_old;
    siz(I)=siz_old(I)+siz(I);
end

%save .dat .opt .temp .whatever
fileID=fopen(file,perm); % write and discard existing content
fwrite(fileID, data ,'double');
fclose(fileID); % write and discard existing content

% size save in .bin
i=find(file=='.',1,'last');
fileID=fopen([file(1:i-1),'.bin'],'w'); % write and discard existing content
fwrite(fileID,siz,'double');
fclose(fileID); % write and discard existing content

end