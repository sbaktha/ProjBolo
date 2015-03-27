T=36;
S=10;
tot=T*S;

z=1;
for k=1:S
for y=1:tot
    if(y==1)
        state(z)=1;
        z=z+1;
    elseif(y>=1 && y<=11)
        state(z)=2;
        z=z+1;
    elseif(y>11 && y<=18)
        state(z)=3;
        z=z+1;
    elseif(y>18 && y<=28)
        state(z)=4;
        z=z+1;
    elseif(y>28 && y<=36)
        state(z)=5;
        z=z+1;
    end
end
end

    