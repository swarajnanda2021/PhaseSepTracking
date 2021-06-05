function k = kmat2cpar( K )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

K=K/K(end,end); % normalize, if was not, essential!

fx=K(1,1);
fy=K(2,2);

px=K(1,3);
py=K(2,3);

s=K(1,2);

k=[fx,s,px,fy,py];

end

