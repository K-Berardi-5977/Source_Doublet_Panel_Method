
function Bp = naca4_airfoil(digits, chord, numPan, spacingType, closedTE)
%NACA4_AIRFOIL Generate boundary points for a NACA 4-digit airfoil.
%
% Bp = naca4_airfoil(digits, chord, numPan, spacingType, closedTE)
%
% digits      : '0012', '2412', etc.
% chord       : chord length
% numPan      : total number of body panels (must be even)
% spacingType : 'cosine' or 'uniform'
% closedTE    : true -> closed trailing edge, false -> standard NACA

    Nside = numPan/2;

    % ---- Parse NACA digits ----
    if isnumeric(digits)
        code = sprintf('%04d', digits);
    else
        code = char(string(digits));
    end
    
    % === NACA 4-Digit Airfoil Classification based on user input ===== %
    m = str2double(code(1))/100;
    p = str2double(code(2))/10;
    t = str2double(code(3:4))/100;

    % ---- Chordwise spacing ----
    switch lower(spacingType)
        case 'cosine'
            beta = linspace(0,pi,Nside+1).';
            x = 0.5*chord*(1 - cos(beta));
        case 'uniform'
            x = linspace(0,chord,Nside+1).';
        otherwise
            error('spacingType must be ''cosine'' or ''uniform''.');
    end

    xc = x/chord;

    % ---- Thickness distribution ----
    if closedTE
        a4 = -0.1036;
    else
        a4 = -0.1015;
    end

    yt = 5*t*chord*( ...
          0.2969*sqrt(xc) ...
        - 0.1260*xc ...
        - 0.3516*xc.^2 ...
        + 0.2843*xc.^3 ...
        + a4*xc.^4 );

    % ---- Camber line ----
    yc = zeros(size(x));
    dycdx = zeros(size(x));

    if m ~= 0
        i1 = xc <= p;
        i2 = xc > p;

        yc(i1) = chord*m/p^2 .* (2*p*xc(i1) - xc(i1).^2);
        yc(i2) = chord*m/(1-p)^2 .* ((1 - 2*p) + 2*p*xc(i2) - xc(i2).^2);

        dycdx(i1) = 2*m/p^2 .* (p - xc(i1));
        dycdx(i2) = 2*m/(1-p)^2 .* (p - xc(i2));
    end

    theta = atan(dycdx);

    % ---- Upper / lower surfaces ----
    xu = x - yt.*sin(theta);
    zu = yc + yt.*cos(theta);

    xl = x + yt.*sin(theta);
    zl = yc - yt.*cos(theta);

    % ---- Assemble CCW boundary ----
    lowerSurf = [flipud(xl), flipud(zl)];
    upperSurf = [xu(2:end), zu(2:end)];

    Bp = [lowerSurf; upperSurf];

    % Close the loop
    if norm(Bp(end,:) - Bp(1,:)) > 1e-12
        Bp = [Bp; Bp(1,:)];
    end
end