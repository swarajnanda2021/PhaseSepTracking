function [ Kmat,Rmat,tvec  ] = pmat2krtm( Pmat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(Pmat)
    [Kmat,Rmat]=rq(Pmat(:,1:3));
    s=diag(sign(diag(Kmat)));%sign(Kmat(end));
    Rmat=s\Rmat;%/s;
    Kmat=Kmat*s;%/s;
    tvec=Kmat\Pmat(:,4);
else
    Kmat=cell(size(Pmat));
    Rmat=cell(size(Pmat));
    tvec=cell(size(Pmat));
    for c=1:length(Pmat)
        [Kmat{c},Rmat{c}]=rq(Pmat{c}(:,1:3));
        s=diag(sign(diag(Kmat)));%sign(Kmat{c}(end));
        Rmat{c}=s\Rmat{c};%/s;
        Kmat{c}=Kmat{c}*s;%/s;
        tvec{c}=Kmat{c}\Pmat{c}(:,4);
    end
end

end

