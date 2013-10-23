function [ get_layer ] = get_layer(z,Thicknesses)
%The bounds of the bottom, (highest z) yes it's supposed to be z down in
%composites, are both counted as part of the last layer, everywhere else
%the bound counts as part of the layer below it, higher z.

zloc=zeros(1,numel(Thicknesses));
zloc(1)=-sum(Thicknesses)/2;
for i=1:numel(Thicknesses)
    zloc(i+1)=zloc(i)+Thicknesses(i);
end
for i=1:numel(Thicknesses)
    if z>zloc(i) && z<=zloc(i+1)
        get_layer=i;
        return
    end
end
if (z==zloc(1))
    get_layer=1;
    return
end
get_layer=0;
return
end

