function [out] = genCombPoly(mat, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub, param)

    if nargin ~= 5
       param.onlyVerticesUniSpaced = false; 
    end
    nVertice = length(mat);
    
%   Add vertices and Uni Spaced to output
    if(isfield(param,'onlyVerticesUniSpaced') && param.onlyVerticesUniSpaced)
       
        
        for i=1:nVertice
            outVert.alphaVecs{i} = zeros(1,nVertice);
            outVert.alphaVecs{i}(1,nVertice)  = 1;
            outVert.polytopicMatrices{i} = mat{i};
        end
        
         if(nVertice > 1 && numPointsUniSpaced> 0)
             out = genPolyUnifSpaced(numPointsUniSpaced, mat);
             
             out.alphaVecs = [out.alphaVecs outVert.alphaVecs];
             out.polytopicMatrices = [out.polytopicMatrices outVert.polytopicMatrices];
         else
             out = outVert;
         end
        
        
    else
        
        if(nVertice > 1)
            %     Generate Matrixes Distributed on the polytope
            out = genPolyUnifSpaced(numPointsUniSpaced, mat);

            % Generate Matrixes Uniformed Spaced between two vertices
            verticesIndexs = 1:1:nVertice;
            combinationVertices = nchoosek(verticesIndexs,2);
            [m, ~] = size(combinationVertices);
            for i=1:1:m
                outPoly2Vertice = genPolyBy2Vertices(mat,combinationVertices(i,:),numPointsBy2Points);

                out.alphaVecs = [out.alphaVecs outPoly2Vertice.alphaVecs];
                out.polytopicMatrices = [out.polytopicMatrices outPoly2Vertice.polytopicMatrices];
            end

            if nVertice > 3 && numbPointsUniSpacedSub ~= 0
            %   Generate Matrixes distributed 3 vertices on the polytope
                verticesIndexs = 1:1:nVertice;
                combinationVertices = nchoosek(verticesIndexs,3);
                [m, ~] = size(combinationVertices);
                for i=1:1:m
                    outPolyUniSpacedSub = genPolyUniSpacedSub(mat,combinationVertices(i,:),numbPointsUniSpacedSub);

                    out.alphaVecs = [out.alphaVecs outPolyUniSpacedSub.alphaVecs];
                    out.polytopicMatrices = [out.polytopicMatrices outPolyUniSpacedSub.polytopicMatrices];
                end
            end 
        else
            out.alphaVecs{1} = 1;
            out.polytopicMatrices{1} = mat{1};
        end
    end
end