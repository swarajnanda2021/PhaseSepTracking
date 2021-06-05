function fbin( file , rows)
%fbin Create .bin that contains shape .dat

fileID=fopen([file,'.bin'],'w'); % write and discard existing content

Cpln=memmapfile([file,'.dat'],...
    'Format','double','Writable',false);

fwrite(fileID, [rows length(Cpln.Data)/rows],'double');
fclose(fileID); % write and discard existing content

end

