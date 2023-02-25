function [Q_bar] = Q_bar(E1,E2,v12,G12,theta) %%% optimize just calc Q
% 
    S = [ 1/E1 , -v12/E1 ,  0;
       -v12/E1 ,    1/E2 ,  0;
           0   ,     0   , 1/G12];
    % Transform Matrix
    s = sind(theta);
    c = cosd(theta);
    T_sigma = [ c^2 , s^2 ,  2*s*c;
                s^2 , c^2 , -2*s*c;
               -s*c , s*c , c^2-s^2];
     
    % S_bar matrix
    S_bar = T_sigma'*S*T_sigma;
    
    Q = inv(S);
    Q_bar = inv(S_bar);
end