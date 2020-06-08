function [kc, indices] = getAllKcFromParallelWaveguide(b, number_kc, resolution)
%getAllKcFromParallelWaveguide - Returns a sorted list of kc and the respective mode in indices.
% Also, in indices, if 3rd index is -1 is a TE, if 3rd index is -2, is a TM mode.
% Arbitrary accuracy (specified by resolution) is provided thanks to the symbolic toolbox.
%
% Syntax: kc = getAllKcFromRectWAveguide(a,b,c,number_kc)

    % Default resolution
    if (nargin == 2)
        resolution = 16;
    end
    % Change the resolution of the symbolic package.
    digitsOld = digits(resolution);
    kc = sym(zeros(1, 2*(number_kc+1)^3));
    indices = zeros(3, 2*(number_kc+1)^3);

    %TEM mode
    kc(1) = 0;
    indices(:,1) = [0 0 0];

    counter = 2;

    for m = 1:number_kc
        % TE mode
        kc(counter) = m*pi/b;
        indices(:,counter) = [m, 0, -1];
        counter = counter + 1;
        % TM mode
        kc(counter) = m*pi/b;
        indices(:,counter) = [m, 0, -2];
        counter = counter + 1;
    end

    % To remove zeros.
    indices = indices(:,1:counter-1);
    kc = vpa(kc(:,1:counter-1),resolution);

    [kc, sorted_indices] = sort(kc);
    indices = indices(:,sorted_indices);

    kc = kc(1:number_kc);
    indices = indices(:,1:number_kc);

    % Restore previous resolution.
    digits(digitsOld)
