function gerarSinLMIAlternadaDerivadaLMILab(A0,B0,E0,C0,D0,deltaA,deltaB,deltaE,deltaC,deltaD,alphas,betas,gammas,deltas,epsilons,nx,nu,nw,ny,sS,sZ,sL,S,Z,L,mu,i,c)
     
    %   S_Vec_nx
    S_Vec_nx = gerarVecRepInt(nx,sS{i});
    [S_Vec_nx,nVar,sSVec_nx] = lmivar(3,S_Vec_nx');
    
    %   L_TVec_nx
     L_TVec_nx = gerarVecRepInt(nx,sL{i});
    [L_TVec_nx,nVar,sL_TVec_nx] = lmivar(3,L_TVec_nx');
    
    %   L_Vec_nx
     L_Vec_nx = gerarVecRepInt(nx,sL{i}');
    [L_Vec_nx,nVar,sL_Vec_nx] = lmivar(3,L_Vec_nx');
    
    %   Z_Vec_nx
     Z_Vec_nx = gerarVecRepInt(nx,sZ{i});
    [Z_Vec_nx,nVar,sZVec_nx] = lmivar(3,Z_Vec_nx');
    
     %   SVec_ny
     S_Vec_ny = gerarVecRepInt(ny,sS{i});
    [S_Vec_ny,nVar,sSVec_ny] = lmivar(3,S_Vec_ny');
    
     %   L_TVec_ny
     L_TVec_ny = gerarVecRepInt(ny,sL{i});
    [L_TVec_ny,nVar,sL_TVec_ny] = lmivar(3,L_TVec_ny');
    
     %   L_Vec_ny
     L_Vec_ny = gerarVecRepInt(ny,sL{i}');
    [L_Vec_ny,nVar,sL_TVec_ny] = lmivar(3,L_Vec_ny');
    
    
     %   Z_Vec_ny
     Z_Vec_ny = gerarVecRepInt(ny,sZ{i});
    [Z_Vec_ny,nVar,sZ_Vec_ny] = lmivar(3,Z_Vec_ny');
    
    %   I_nwVec_nx
     I_nxVec_nxw = gerarVecRepInt(nw,eye(nx));
     I_nxVec_nxw = I_nxVec_nxw';
    


    lmiterm([c 1 1 S{i}],1,A0','s');
    lmiterm([c 1 1 L{i}],1,B0','s');
    for k=1:nx
       for j=1:nx
            ek = generatesEi(k,nx);
            lmiterm([c 1 1 alphas(k,j)],2*abs(deltaA(k,j))^2*(ek*ek'),1);
       end
    end
     for k=1:nx
       for j=1:nu
            ek = generatesEi(k,nx);
            lmiterm([c 1 1 betas(k,j)],2*abs(deltaB(k,j))^2*(ek*ek'),1);
       end
    end
    
    
    lmiterm([c 2 1 -L{i}],1,A0');
    lmiterm([c 2 1 -Z{i}],1,B0');
    
    
    lmiterm([c 3 1 0],E0');
    lmiterm([c 3 3 mu],-1,1);
    for k=1:nx
       for j=1:nw
             hj = generatesEi(j,nw);
            lmiterm([c 3 3 epsilons(k,j)],abs(deltaE(k,j))^2*(hj*hj'),1);
       end
    end
    
    
    lmiterm([c 4 1 S{i}],C0,1);
    lmiterm([c 4 1 -L{i}],D0,1);
    
    lmiterm([c 4 2 L{i}],C0,1);
    lmiterm([c 4 2 Z{i}],D0,1);
    
    lmiterm([c 4 4 0],-1);
    for k=1:ny
       for j=1:nx
            gk = generatesEi(k,ny);
            lmiterm([c 4 4 gammas(k,j)],2*abs(deltaC(k,j))^2*(gk*gk'),1);
       end
    end
    for k=1:ny
       for j=1:nu
             gk = generatesEi(k,ny);
            lmiterm([c 4 4 deltas(k,j)],2*abs(deltaD(k,j))^2*(gk*gk'),1);
       end
    end
    
    
    
    lmiterm([c 5 1 S_Vec_nx], 1, 1);
    for k=1:nx
       for j=1:nx
           ei = generatesEi(((k-1)*nx)+j,(nx*nx)); 
           lmiterm([c 5 5 alphas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 6 1 L_TVec_nx], 1, 1);
    for k=1:nx
       for j=1:nu
           ei = generatesEi(((k-1)*nu)+j,(nx*nu)); 
           lmiterm([c 6 6 betas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 7 2 L_Vec_nx], 1, 1);
    for k=1:nx
       for j=1:nx
           ei = generatesEi(((k-1)*nx)+j,(nx*nx)); 
           lmiterm([c 7 7 alphas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 8 2 Z_Vec_nx], 1, 1);
    for k=1:nx
       for j=1:nu
           ei = generatesEi(((k-1)*nu)+j,(nx*nu)); 
           lmiterm([c 8 8 betas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 9 1 S_Vec_ny], 1, 1);
    for k=1:ny
       for j=1:nx
           ei = generatesEi(((k-1)*nx)+j,(ny*nx)); 
           lmiterm([c 9 9 gammas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 10 1 L_TVec_ny], 1, 1);
    for k=1:ny
       for j=1:nu
           ei = generatesEi(((k-1)*nu)+j,(ny*nu)); 
           lmiterm([c 10 10 deltas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 11 2 L_Vec_ny], 1, 1);
    for k=1:ny
       for j=1:nx
           ei = generatesEi(((k-1)*nx)+j,(ny*nx)); 
           lmiterm([c 11 11 gammas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 12 2 Z_Vec_ny], 1, 1);
    for k=1:ny
       for j=1:nu
           ei = generatesEi(((k-1)*nu)+j,(ny*nu)); 
           lmiterm([c 12 12 deltas(k,j)],-1*(ei*ei'),1);
       end
    end
    
    lmiterm([c 13 1 0], I_nxVec_nxw);
    for k=1:nx
       for j=1:nw
           ei = generatesEi(((k-1)*nw)+j,(nx*nw)); 
           lmiterm([c 13 13 epsilons(k,j)],-1*(ei*ei'),1);
       end
    end
end