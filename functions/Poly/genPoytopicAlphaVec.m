function [alphaVec] = genPoytopicAlphaVec(nVert)

    alphaVec(1)=1-rand^(1/(nVert-1));
    
    for k=2:nVert-1
        alphaVec(k)=(1-sum(alphaVec(1:k-1)))*(1-rand^(1/(nVert-k)));
    end
    
    alphaVec(nVert)=1-sum(alphaVec(1:nVert-1));
    
end