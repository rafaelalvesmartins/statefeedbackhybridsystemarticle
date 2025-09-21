function matInfSup = INTLAB2InfSup(mat)
    
    % Not infsup, then nominal value
    if(isnumeric(mat))
        matInfSup.inf = mat;
        matInfSup.sup = mat;
    else % infsup type
        matInfSup.inf = inf(mat);
        matInfSup.sup = sup(mat);
    end
    
end