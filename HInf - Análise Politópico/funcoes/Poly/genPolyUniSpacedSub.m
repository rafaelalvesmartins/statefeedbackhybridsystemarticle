function [out] = genPolyUniSpacedSub(mat,vecIndex,numPoints)
   
   for i=1:1:numPoints     
        alphaVec = genPoytopicAlphaVec(3);    
        polytopicMatrix = genStructure(mat);
        alphaVecs2Save = zeros(1, length(mat));
        
        % Convex combination
        for j=1:3
            polytopicMatrix = sumStructure(polytopicMatrix, multStructure(mat{vecIndex(j)}, alphaVec(j)));
            
            alphaVecs2Save(1,vecIndex(j)) = alphaVec(j);
        end
        
        alphaVecs{i} = alphaVecs2Save;
        polytopicMatrices{i} = polytopicMatrix;
    end   
    
    out.alphaVecs = alphaVecs;
    out.polytopicMatrices = polytopicMatrices; 
end