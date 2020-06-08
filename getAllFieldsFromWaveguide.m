function [electricField, magneticField] = getAllFieldsFromWaveguide(waveguideType, a, b, epsilonR, frequency, modeWG, normFactor, x, y, z)
    %getAllFieldsFromWaveguide - Returns the evaluation of the fields for a
    % rectangular waveguide of dimensions a x b. In modeWG,
    % an array of 3x1 is needed: [m, n, TEflag] for
    % rectangular waveguide, 2x1 for parallel plate.
    % TEflag might be either -1 (TE mode) or -2 (TM mode).
    % For parallel plate, if TEFlag == 0, TEM mode is returned
    % The returned values are symbolic expressions in x, y.
    % These are thought to serve as patterns in the
    % transverse waveport, so z = 0.
    % muR = 1 is assumed.
    % waveguideType might be:
    % Rectangular
    % ParallelPlate
    %
    % Syntax: kc = getAllFieldsFromWaveguide(a, b, epsilonR, modeWG)
    %

    if (nargin == 6)
        normFactor = 1;
    end
    % If coordinates are not provided, just return the symbolic computation.
    if (nargin<=7)
        xyObj = math.AffineCoordinates('M', 2, false);
        xy = xyObj.coordinates;
        x = xy(1);
        y = xy(2);
    end
    m = modeWG(1);
    betaX = m*pi/a;
    [epsilon0, mu0, c0, eta0] = settings.getVacuumConstants();
    epsilon = epsilonR*epsilon0;
    omega = 2*pi*frequency;

    switch waveguideType

    case 'ParallelPlate'
        E0 = 1/normFactor;
        switch modeWG(3)
        % TEM mode
        case 0
            betaZ = omega/c0*sqrt(epsilonR);
            if (b>a)
                % This comes from a swapping in the geometrical mapping of the waveport.
                Ey = E0;
                Ex = 0;
                Ez = 0;
                Hy = 0;
                Hx = E0/eta0;
                Hz = 0;
            else
                Ex = E0;
                Ey = 0;
                Ez = 0;
                Hx = 0;
                Hy = E0/eta0;
                Hz = 0;
            end
        % TE mode
        case -1
            % It may happen that the axis are interchanged (when you have a swapping in build2DMeshFrom3DMesh)
            if (m == 0)
                n = modeWG(2);
                betaY = n*pi/b;
                betaCutoff = n*pi/b;
                betaZ = sqrt((omega/c0*sqrt(epsilonR))^2 - betaCutoff^2);

                Ey = 0;
                Ex = -E0*sin(betaY*y);
                Ez = 0;
                Hy = -E0*betaZ/(omega*mu0)*sin(betaY*y);
                Hx = 0;
                Hz = 1i*E0*betaY/(omega*mu0)*cos(betaY*y);
            else
                betaCutoff = betaX;
                betaZ = sqrt((omega/c0*sqrt(epsilonR))^2 - betaCutoff^2);

                Ex = 0;
                Ey = E0*sin(betaX*x);
                Ez = 0;
                Hx = -E0*betaZ/(omega*mu0)*sin(betaX*x);
                Hy = 0;
                Hz = 1i*E0*betaX/(omega*mu0)*cos(betaX*x);
            end
        % TM mode
        case -2
            % It may happen that the axis are interchanged.
            if (m == 0)
                n = modeWG(2);
                betaY = n*pi/b;
                betaCutoff = n*pi/b;
                betaZ = sqrt((omega/c0*sqrt(epsilonR))^2 - betaCutoff^2);

                Ey = E0*cos(betaY*y);
                Ex = 0;
                Ez = 1i*betaY/betaZ*E0*sin(betaY*y);
                Hy = 0;
                Hx = -omega*epsilon0*epsilonR/betaZ*E0*cos(betaY*y);
                Hz = 0;
            else
                betaCutoff = betaX;
                betaZ = sqrt((omega/c0*sqrt(epsilonR))^2 - betaCutoff^2);

                Ex = E0*cos(betaX*x);
                Ey = 0;
                Ez = 1i*betaX/betaZ*E0*sin(betaX*x);
                Hx = 0;
                Hy = omega*epsilon0*epsilonR/betaZ*E0*cos(betaX*x);
                Hz = 0;
            end
        otherwise
            error('Check sintax of modeWG')
        end

    case 'Rectangular'
        n = modeWG(2);
        betaY = n*pi/b;
        betaCutoff = sqrt(betaX^2 + betaY^2);
        betaZ = sqrt((omega/c0*sqrt(epsilonR))^2 - betaCutoff^2);

        switch modeWG(3)
        % TE mode
        case -1
            Amn = 1/normFactor;
            Ex = Amn*betaY/epsilon*cos(betaX*x)*sin(betaY*y);
            Ey = -Amn*betaX/epsilon*sin(betaX*x)*cos(betaY*y);
            Ez = 0;
            Hx = Amn*betaX*betaZ/(omega*mu0*epsilon)*sin(betaX*x)*cos(betaY*y);
            Hy = Amn*betaY*betaZ/(omega*mu0*epsilon)*cos(betaX*x)*sin(betaY*y);
            Hz = -1i*Amn*betaCutoff^2/(omega*mu0*epsilon)*cos(betaX*x)*cos(betaY*y);

            % A10 = sqrt(4*omega*mu0*epsilon^2/(betaCutoff^2*betaZ*a*b));
            % Ex = 0;
            % Ey = -A10/epsilon*pi/a*sin(betaX*x);
            % Ez = 0;
            % Hx = A10*betaZ/(omega*mu0*epsilon)*pi/a*sin(betaX*x);
            % Hy = 0;
            % Hz = -1i*A10/(omega*mu0*epsilon)*(pi/a)^2*cos(betaX*x);

        % TM mode
        case -2
            Bmn = 1/normFactor;
            Hx = Bmn*betaY/mu0*sin(betaX*x)*cos(betaY*y);
            Hy = -Bmn*betaX/mu0*cos(betaX*x)*sin(betaY*y);
            Hz = 0;
            Ex = -Bmn*betaX*betaZ/(omega*mu0*epsilon)*cos(betaX*x)*sin(betaY*y);
            Ey = -Bmn*betaY*betaZ/(omega*mu0*epsilon)*sin(betaX*x)*cos(betaY*y);
            Ez = -1i*Bmn*betaCutoff^2/(omega*mu0*epsilon)*sin(betaX*x)*sin(betaY*y);

        otherwise
            error('Check sintax of modeWG')
        end

    otherwise
        error ('Waveguide not supported yet.')
    end

    electricField = [Ex, Ey, Ez];
    magneticField = [Hx, Hy, Hz];
    % Add propagation constant
    if (nargin == 10)
        electricField = electricField*exp(-1i*betaZ*z);
        magneticField = magneticField*exp(-1i*betaZ*z);
    end