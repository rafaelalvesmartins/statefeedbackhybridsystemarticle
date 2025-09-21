function matVec = gerarVecRepInt(rep,mat)

    for i=1:rep
        if(i==1)
            matVec = mat;
        else
            matVec = [matVec mat];
        end
    end
    
end