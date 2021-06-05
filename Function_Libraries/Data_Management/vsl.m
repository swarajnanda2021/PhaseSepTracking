function var_slash=vsl
%op_slash is variable slash by operating system
%
%   Junaid
%
%   other system can be added if exist..?

if isunix
    var_slash='/' ;
elseif ispc
    var_slash='\' ;
end

end
