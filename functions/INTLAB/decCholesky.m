function L=decCholesky(A)

    n = length( A ); 
    L = midrad(zeros(n,n),0);

    for i=1:n
        L(i, i) = sqrt( A(i, i) - L(i, :)*L(i, :)' );

        for j=(i + 1):n
            L(j, i) = ( A(j, i) - L(i, :)*L(j, :)' )/L(i, i);
        end
    end
    
    % ERROR      
    if(sum(sum(isnan(L)))>0)
        display('Error in Cholesky decompositon, there are NaN element: ');
        L
    end
end
 