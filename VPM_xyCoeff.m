function [Nx , Ny] = VPM_xyCoeff(xi, yi, Xj, Yj, phi, S)
numPan = length(xi);

Nx = zeros(numPan,1);                                                       % Initialize Mx integral array
Ny = zeros(numPan,1);                                                       % Initialize My integral array

% Compute Mx and My
for i=1:numPan
    for j = 1:numPan
    % Loop over the j panels
    % Compute intermediate values
    A  = -(xi(i)-Xj(j))*cos(phi(j))-(yi(i)-Yj(j))*sin(phi(j));                    % A term
    B  = (xi(i)-Xj(j))^2+(yi(i)-Yj(j))^2;                                         % B term
    Cy = -cos(phi(j));                                                      % C term (X-direction)
    Dx = -(yi(i) - Yj(j));                                                        % D term (X-direction)
    Cx = sin(phi(j));                                                      % C term (Y-direction)
    Dy = xi(i) - Xj(j);                                                        % D term (Y-direction)
    E  = sqrt(B-A^2);                                                       % E term
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Nx, Ref [1]
    term1 = 0.5*Cx*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Nx equation
    term2 = ((Dx-A*Cx)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Nx equation
    Nx(j) = term1 + term2;                                                  % X-direction geometric integral
    
    % Compute Ny, Ref [1]
    term1 = 0.5*Cy*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Ny equation
    term2 = ((Dy-A*Cy)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Ny equation
    Ny(j) = term1 + term2;                                                  % Y-direction geometric integral
    
    % Zero out any NANs, INFs, or imaginary numbers
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
    end
end
end