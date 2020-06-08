function [electricField, magneticField] = getAllFieldsFromRectCavity(a, b, c, epsilonR, frequency, modeWG, x, y, z)
    %getAllFieldsFromRectCavity - Returns the evaluation of the fields for a
    % rectangular Cavity of dimensions a x b x c, defined in [0,a]x[0,b]x[0,c]. In modeWG,
    % an array of 4x1 is needed: [m, n, p, TEflag]
    % TEflag might be either -1 (TE mode) or -2 (TM mode).
    % The returned values are symbolic expressions in x, y, z
    % muR = 1 is assumed.
    %
    % Syntax: [electricField, magneticField] = getAllFieldsFromRectCavity(a, b, c, epsilonR, modeWG)
    %

    if (nargin==6)
        xyzObj = math.AffineCoordinates('M', 3, false);
        xyz = xyzObj.coordinates;
        x = xyz(1);
        y = xyz(2);
        z = xyz(3);
    end
    m = modeWG(1);
    n = modeWG(2);
    p = modeWG(3);

    betaX = m*pi/a;
    betaY = n*pi/b;
    betaZ = p*pi/c;
    [epsilon0, mu0, c0] = settings.getVacuumConstants();
    epsilon = epsilonR*epsilon0;
    mu = mu0*1;
    omega = 2*pi*frequency;
    betaProp = omega/c0*sqrt(epsilonR);

    switch modeWG(4)
    % TE mode
    case -1
        Ex = betaY/epsilon*cos(betaX*x)*sin(betaY*y)*sin(betaZ*z);
        Ey = -betaX/epsilon*sin(betaX*x)*cos(betaY*y)*sin(betaZ*z);
        Ez = 0;
        Hx = 1i*betaX*betaZ/(omega*mu*epsilon)*sin(betaX*x)*cos(betaY*y)*cos(betaZ*z);
        Hy = 1i*betaY*betaZ/(omega*mu*epsilon)*cos(betaX*x)*sin(betaY*y)*cos(betaZ*z);
        Hz = -1i*(betaProp^2-betaZ^2)/(omega*mu*epsilon)*cos(betaX*x)*cos(betaY*y)*sin(betaZ*z);

    % TM mode
    case -2
        % Not from Balanis (but for analogy)
        Hx = betaY/mu0*sin(betaX*x)*cos(betaY*y)*cos(betaZ*z);
        Hy = -betaX/mu0*cos(betaX*x)*sin(betaY*y)*cos(betaZ*z);
        Hz = 0;
        Ex = 1i*betaX*betaZ/(omega*mu0*epsilon)*cos(betaX*x)*sin(betaY*y)*sin(betaZ*z);
        Ey = 1i*betaY*betaZ/(omega*mu0*epsilon)*sin(betaX*x)*cos(betaY*y)*sin(betaZ*z);
        Ez = -1i*(betaProp^2-betaZ^2)/(omega*mu0*epsilon)*sin(betaX*x)*sin(betaY*y)*cos(betaZ*z);

    otherwise
        error('Check sintax of modeWG')
    end

    electricField = [Ex, Ey, Ez];
    magneticField = [Hx, Hy, Hz];