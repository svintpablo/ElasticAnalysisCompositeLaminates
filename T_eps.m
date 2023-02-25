function [sG,sL] = T_eps(sG,sL,theta)
% Function ideally takes theta, sG, and sL as epsilon (strain) global/local respectively. 3
% of the inputs to sL/sG should be defined, and three symbolic. Will solve for the
% remaining 3 values and return sG and sL.

% defining transformation matrix
%    theta = 45;
    s = sind(theta);
    c = cosd(theta);
    T_sigma = [ c^2 , s^2 ,  2*s*c;
                s^2 , c^2 , -2*s*c;
               -s*c , s*c , c^2-s^2];
    T_epsilon = T_sigma'^(-1);


    %example code 
%     sG = sym('s_g',[3 1]);
%     sL = sym('s_l',[3 1]);
% 
%     sL(1) = 3.5;
%     sG(2) = 2;
%     sL(3) = -0.5;

    eqns = T_epsilon*sG-sL;           % symbolic equations
    vars = symvar(eqns);    % get all unknown variables
    [A,B] = equationsToMatrix(eqns,vars); % creates unique system
    A = double(A);
    B = double(B);
    res = A\B;  %solves system for remaining variables
    
    % assigns values back to respective 'slots'
    for i = 1:3
        for j = 1:3
            if strcmp(char(vars(i)),char(sG(j)))
                sG(j) = res(i);
            end
        end
        for j = 1:3
            if strcmp(char(vars(i)),char(sL(j)))
                sL(j) = res(i);
            end
        end
    end
    % converts from symbolic array to double array
    sG = double(sG);
    sL = double(sL);
end
