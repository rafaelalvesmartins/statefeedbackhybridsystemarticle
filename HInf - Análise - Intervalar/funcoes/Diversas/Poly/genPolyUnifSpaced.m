function [out] = genPolyUnifSpaced(numPointsUniSpaced, mat)
    nVertice = length(mat);

    for i=1:(numPointsUniSpaced)
        
        alphaVec = genPoytopicAlphaVec(nVertice);
       
        polytopicMatrix = genStructure(mat);
%         Convex combination
        for j=1:nVertice
            polytopicMatrix = sumStructure(polytopicMatrix, multStructure(mat{j}, alphaVec(j)));
        end
        
        alphaVecs{i} = alphaVec;
        polytopicMatrices{i} = polytopicMatrix;
        
    end   
    
    out.alphaVecs = alphaVecs;
    out.polytopicMatrices = polytopicMatrices;

end