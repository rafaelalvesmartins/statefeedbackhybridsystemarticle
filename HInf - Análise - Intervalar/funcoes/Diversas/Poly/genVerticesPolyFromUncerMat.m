function polytopic = genVerticesPolyFromUncerMat(expDisc,prec)
    
    diffMat = abs(expDisc.inf - expDisc.sup);
    [variationRows,variationCols] = find(diffMat > prec); % Inputs nonzero
    
    uncertainMatrix = zeros(length(variationRows),2);
    uncertainMatrix(:,1) =  expDisc.inf(sub2ind(size(expDisc.inf),variationRows(:),variationCols(:)))';
    uncertainMatrix(:,2) =  expDisc.sup(sub2ind(size(expDisc.sup),variationRows(:),variationCols(:)))';
    
    text = 'combination = combvec(';
    for i=1:length(variationRows)
        if i ~= length(variationRows)
            text = [text  sprintf('uncertainMatrix(%d,:),',i)];
        else
            text = [text sprintf('uncertainMatrix(%d,:)',i)];
        end
    end
    text =[text ');'];
    
    eval(text);
    
    indexPolytopic = 1;
    polytopic{indexPolytopic} = expDisc.inf;
    if(length(combination)>=2) % Has at least one one uncertainty     
        for i=2:length(combination)
            indexPolytopic = indexPolytopic + 1;
            polytopic{indexPolytopic} = expDisc.inf;
            for j=1:length(variationRows)
                polytopic{indexPolytopic}(variationRows(j),variationCols(j)) = combination(j,i);
            end
        end  
    end
    
end