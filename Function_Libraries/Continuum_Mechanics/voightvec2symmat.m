function M = voightvec2symmat(voightvec)
%symmat2voightvec Symmetric 3x3 Matrix to Voight vector notation
%
%   Input
%       voightvec Voightvector notation for 3xN Matrix
%   
%   Output
%       M Some (2x2) 3x3xN Matrix

% correct input
if size(voightvec,1)==3
    voightvec=[voightvec(1:2,:)
        zeros(3,size(voightvec,2))
        voightvec(3,:)]; % zeropad 3
end

% reshape data for matrix
voightvec=reshape(voightvec,size(voightvec,1),1,[]);

% strain matrix
M=[ voightvec(1,:,:)   voightvec(6,:,:)/2 voightvec(5,:,:)/2
    voightvec(6,:,:)/2 voightvec(2,:,:)   voightvec(4,:,:)/2
    voightvec(5,:,:)/2 voightvec(4,:,:)/2 voightvec(3,:,:)  ];

end

