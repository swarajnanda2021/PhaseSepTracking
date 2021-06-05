clear
close all

a=readimx('2D-PIV particle.imx');


size(a.Frames)


figure
Frame=a.Frames{1};
Components = Frame.Components;
for i = 1:size(Components,1),
    Planes = Components{i}.Planes;
    nz = size(Planes,1);
    if nz==0,
        disp('no Planes')
    elseif nz==1,
        %D = showPlane( Planes{1}, Frame.Scales );
        imagesc( Planes{1})
    else
        %D = showVolume( Planes, Frame.Scales );
    end
end


figure
Frame=a.Frames{2};
Components = Frame.Components;
for i = 1:size(Components,1),
    Planes = Components{i}.Planes;
    nz = size(Planes,1);
    if nz==0,
        disp('no Planes')
    elseif nz==1,
        %D = showPlane( Planes{1}, Frame.Scales );
        imagesc( Planes{1})
    else
        %D = showVolume( Planes, Frame.Scales );
    end
end


