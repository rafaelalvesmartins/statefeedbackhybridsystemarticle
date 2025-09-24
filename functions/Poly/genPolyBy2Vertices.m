function [out] = genPolyBy2Vertices(mat,combinationVertices,numPoints)
    
    increment = 1/(numPoints-1);
    
    alpha1 = 0;
    for i=1:numPoints
        alpha2 = 1 - alpha1;
       
        polytopicMatrix = sumStructure(multStructure(mat{combinationVertices(1)}, alpha1), multStructure(mat{combinationVertices(2)}, alpha2));
        
        alphaVec = zeros(1,length(mat));
        
        alphaVec(1, combinationVertices(1)) = alpha1;
        alphaVec(1, combinationVertices(2)) = alpha2;
        
        alphaVecs{i} = alphaVec;
        polytopicMatrices{i} = polytopicMatrix;
        
        alpha1 = alpha1 + increment;
    end
    
    out.alphaVecs = alphaVecs;
    out.polytopicMatrices = polytopicMatrices;
end