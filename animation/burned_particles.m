function Cp = burned_particles(Xp, Yp, Cp, Xc, Yc, flame)

% Is any change to Cp needed? 
if flame == 0; return;
else

    % Which particles are inside the cylinder?
    [in, ~] = inpolygon(Xp, Yp, Xc, Yc);

    % Loop over all particles & set color
    for i = 1 : length(Cp)
        if in(i); Cp(i) = max(Cp(i), flame); end

    end
end
