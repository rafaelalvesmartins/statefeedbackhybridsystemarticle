function [out] = genCombPoly(mat, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub,param)

    nVertice = length(mat);
    
    if(param.onlyVertice)
         for i=1:nVertice
             out.alphaVecs{i} = i;
             out.polytopicMatrices{i} = mat{i};
         end
    else
    
        if(nVertice > 1)
            %     Generate Matrixes Distributed on the polytope
            out = genPolyUnifSpaced(numPointsUniSpaced, mat);

            % Generate Matrixes Uniformed Spaced between two vertices
            verticesIndexs = 1:1:nVertice;
            combinationVertices = nchoosek(verticesIndexs,2);
            [m, ~] = size(combinationVertices);
            numPointsBy2Points4Each2Verts = round(numPointsBy2Points/m);
            if(numPointsBy2Points4Each2Verts<1)
                numPointsBy2Points4Each2Verts = 3;
                m = round(numPointsBy2Points/3);
            end
            for i=1:1:m
                outPoly2Vertice = genPolyBy2Vertices(mat,combinationVertices(i,:),numPointsBy2Points4Each2Verts);

                out.alphaVecs = [out.alphaVecs outPoly2Vertice.alphaVecs];
                out.polytopicMatrices = [out.polytopicMatrices outPoly2Vertice.polytopicMatrices];
            end

            if nVertice >= 3 && numbPointsUniSpacedSub ~= 0
            %   Generate Matrixes distributed 3 vertices on the polytope
                verticesIndexs = 1:1:nVertice;
                combinationVertices = nchoosek(verticesIndexs,3);
                [m, ~] = size(combinationVertices);
                numPointsBy3Points4Each3Verts = round(numbPointsUniSpacedSub/m);
                if(numPointsBy3Points4Each3Verts<1)
                    numPointsBy3Points4Each3Verts = 5;
                    m = round(numPointsBy2Points/5);
                end
                for i=1:1:m
                    outPolyUniSpacedSub = genPolyUniSpacedSub(mat,combinationVertices(i,:),numPointsBy3Points4Each3Verts);

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