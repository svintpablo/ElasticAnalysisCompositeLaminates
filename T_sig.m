function [sL] = T_sig(sG,theta)
% Function ideally takes theta, sG, and sL as sig global/local respectively. 3
% of the inputs to sL/sG should be defined, and three symbolic. Will solve for the
% remaining 3 values and return sG and sL.

% defining transformation matrix
%    theta = 45;
    s = sind(theta);
    c = cosd(theta);
    T_sigma = [ c^2 , s^2 ,  2*s*c;
                s^2 , c^2 , -2*s*c;
               -s*c , s*c , c^2-s^2];
    sL = T_sigma*sG;
end