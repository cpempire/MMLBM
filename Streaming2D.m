function [fIn,Lattice] = Streaming2D(fIn,fOut,Lattice,LS2D,BCs,WhichBlock)
switch WhichBlock
    case { 'P', 'Q' }
        %% Streaming Step
        in = 1;  out = Lattice.Nx;       % in and out boundary x-direction indices
        col = 2:Lattice.Ny-1;            % in and out boundary y-direction indices
        up = Lattice.Ny; down = 1;       % up and down wall boundary y-direction indices
        row = 2:Lattice.Nx-1;            % up and down wall boundary x-direction indices

          % INNER PART
        for i = 1:9
            fIn(i,row+LS2D.cx(i),col+LS2D.cy(i)) = fOut(i,row,col); 
        end

          % INLET and OUTLET
        for i = [1,2,3,5,6,9]
            fIn(i,in+LS2D.cx(i),col+LS2D.cy(i)) = fOut(i,in,col);
            fIn(LS2D.opp(i),out+LS2D.cx(LS2D.opp(i)),col+LS2D.cy(LS2D.opp(i))) = fOut(LS2D.opp(i),out,col);
        end
          % UP WALL and DOWN WALL
        for i = [1,2,4,5,8,9]
            fIn(i,row+LS2D.cx(i),up+LS2D.cy(i)) = fOut(i,row,up);
            fIn(LS2D.opp(i),row+LS2D.cx(LS2D.opp(i)),down+LS2D.cy(LS2D.opp(i))) = fOut(LS2D.opp(i),row,down);
        end

          % INLET-UP and OUTLET-DOWN CORNER
        for i = [1,2,5,9]
            fIn(i,in+LS2D.cx(i),up+LS2D.cy(i)) = fOut(i,in,up);
            fIn(LS2D.opp(i),out+LS2D.cx(LS2D.opp(i)),down+LS2D.cy(LS2D.opp(i))) = fOut(LS2D.opp(i),out,down);
        end
          % INLET-DOWN and OUTLET-UP CORNER
        for i = [1,2,3,6]
            fIn(i,in+LS2D.cx(i),down+LS2D.cy(i)) = fOut(i,in,down);
            fIn(LS2D.opp(i),out+LS2D.cx(LS2D.opp(i)),up+LS2D.cy(LS2D.opp(i))) = fOut(LS2D.opp(i),out,up);
        end

    case 'C'
          % INNER PART
        for i = 1:9
            fIn(i,find(Lattice.InnerRegion==1)+LS2D.cx(i)+LS2D.cy(i)*Lattice.Nx) = fOut(i,find(Lattice.InnerRegion==1)); 
        end
        

%           % Boundary 
%         for k = 1:length(Lattice.Boundary)
%             for i = find(BCs.C(1:9,Lattice.Boundary(k))==1)'
%                 if (Lattice.Boundary(k)+LS2D.cx(LS2D.opp(i))+LS2D.cy(LS2D.opp(i))*Lattice.Nx<Lattice.Nx*Lattice.Ny...
%                         &&Lattice.Boundary(k)+LS2D.cx(LS2D.opp(i))+LS2D.cy(LS2D.opp(i))*Lattice.Nx>0)
%                         fIn(LS2D.opp(i), Lattice.Boundary(k)+LS2D.cx(LS2D.opp(i))+LS2D.cy(LS2D.opp(i))*Lattice.Nx) = fOut(LS2D.opp(i),Lattice.Boundary(k));
%                 end
%             end
%         end

        
        for k = Lattice.Boundary'
            for i = 1:9
                if k+LS2D.cx(i)+LS2D.cy(i)*Lattice.Nx<Lattice.Nx*Lattice.Ny&&...
                        k+LS2D.cx(i)+LS2D.cy(i)*Lattice.Nx>0 
                    fIn(i,k+LS2D.cx(i)+LS2D.cy(i)*Lattice.Nx) = fOut(i,k);
                end
            end
        end
                
        % obstacle as curved boundary 
        for i = Lattice.CurveBoundary'
            for k = 1:9
                if BCs.C(k,i)>0
                    fIn(LS2D.opp(k),i) = 1/(BCs.C(k,i)+1)*( (1-BCs.C(k,i))*fIn(k,i) + BCs.C(k,i)*fIn(k,i+LS2D.cx(k)+LS2D.cy(k)*Lattice.Nx) ) ...
                                            + BCs.C(k,i)/(BCs.C(k,i)+1)* fIn(LS2D.opp(k),i+LS2D.cx(LS2D.opp(k))+LS2D.cy(LS2D.opp(k))*Lattice.Nx );
                end
            end
        end
        
        fIn(:,Lattice.Obstacle) = 0;
        fEq(:,Lattice.Obstacle) = 0;
        
    otherwise
        disp('please specify which block')
end


Lattice.Rho = reshape(sum(fIn),1,Lattice.Nx,Lattice.Ny);
Lattice.Ux  = reshape(LS2D.cx*reshape(fIn,9,Lattice.Nx*Lattice.Ny),1,Lattice.Nx,Lattice.Ny)./Lattice.Rho;                             
Lattice.Uy  = reshape(LS2D.cy*reshape(fIn,9,Lattice.Nx*Lattice.Ny),1,Lattice.Nx,Lattice.Ny)./Lattice.Rho;