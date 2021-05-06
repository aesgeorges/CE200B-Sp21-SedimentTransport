x = linspace(0, H);
y = linspace(0, H);

function x = xpos(x, up)
    t = 0.01;
    x = x + t*up;
end

function z = zpos(z, up)
    t = 0.01;
    z = z + t*up;
end