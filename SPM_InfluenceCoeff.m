function [I, J] = SPM_InfluenceCoeff(xi, yi, Xj, Yj, phi, S, numPan)
numPan = length(xi); %iteration variable
I = zeros(numPan, numPan);
J = zeros(numPan, numPan);

% xi & yi - vectors containing the x and y coordinates of the control
% points

%Defining Integral Terms:
for i = 1:numPan %iterating over the jth control point
    for j = 1:numPan %for each control point, iterate over j=1:n panels 
        if (j~=i) %avoid singularity
            %all of the coefficient terms determined in the solution of the
            %geometric integrals 
            A = -(xi(i)-Xj(j))*cos(phi(j))-(yi(i)-Yj(j))*sin(phi(j));
            B = (xi(i)-Xj(j))^2 +(yi(i)-Yj(j))^2;
            C = sin(phi(i)-phi(j));
            C_s = -cos(phi(i)-phi(j));
            D = -(xi(i)-Xj(j))*sin(phi(i))+(yi(i)-Yj(j))*cos(phi(i));
            D_s = (xi(i)-Xj(j))*cos(phi(i)) + (yi(i)-Yj(j))*sin(phi(i));
            E = sqrt(B-A^2);
            Sj = S(j);

            if ~isreal(E)
                E = 0;
            end
            %compute normal geometric integral for each panel at the ith
            %point
            I(i,j)= (C/2)*log((Sj^2 + 2*A*Sj + B)/B) + ((D - A*C)/E)*(atan2(Sj+A, E) ...
                -atan2(A, E)); 

            %compute tangent (surface) geometric interal for each panel at
            %the ith point
            % J(i,j) = (C_s/2)*log((Sj^2 + 2*A*Sj + B)/B) + ((D_s - A*C_s)/E)*(atan2(Sj+A, E) ...
            %     -atan2(A, E));
            J(i,j) = ((D-A*C)/(2*E))*log((Sj^2 + 2*A*Sj + B)/B) - C*(atan2(Sj+A, E) ...
                -atan2(A, E));
        end
    if (isnan(I(i,j)) || isinf(I(i,j)) || ~isreal(I(i,j)))
        I(i,j) = 0;
    end

 

    if (isnan(J(i,j)) || isinf(J(i,j)) || ~isreal(J(i,j)))
        J(i,j) = 0;
    end

    end

end
end