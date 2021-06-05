function [ Pmat ] = krtm2pmat( Kmat,Rmat,tvec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(Kmat) || ~iscell(Rmat) || ~iscell(tvec) 
    Pmat=Kmat*[Rmat tvec];
else
    Pmat=cell(size(Kmat));
    for c=1:length(Kmat)
        Pmat{c}=Kmat{c}*[Rmat{c} tvec{c}];
    end
end

end

