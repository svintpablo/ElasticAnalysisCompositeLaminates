%% ** Thermal loads **
% Pablo Tejeda
% 
%% 
%
function [N_T , M_T, alpha_k] = ThermalLoads(E1,E2,v12,G12, stack, t, alpha, Del_T, stackmap)
    n = numel(stack);
        h = 0;
    for i = 1:n
        mat = stackmap(i);
        h = h + t(mat);
    end

    z_0 = -h/2;

    N_T = [0 0 0]';
    M_T = [0 0 0]';
    alpha_k = zeros(3,n);



    for i = 1:n
        mat = stackmap(i);

        theta = stack(i);
        s = sind(theta);
        c = cosd(theta);

        z_top = z_0+t(mat)*i;
        z_bot = z_0+t(mat)*(i-1);
        
        Q_bar_k = Q_bar(E1(mat),E2(mat),v12(mat),G12(mat),theta); % consider collecting Qbar and alpha values

        alpha_k(1,i) = alpha(1, mat)*c^2 +alpha(2, mat)*s^2;
        alpha_k(2,i) = alpha(1, mat)*s^2 +alpha(2, mat)*c^2;
        alpha_k(3,i) = 2*(alpha(1, mat)-alpha(2, mat))*c*s;

        
        N_T = N_T + Del_T*Q_bar_k*alpha_k(:,i)*t(mat);
        M_T = M_T + .5*Del_T*Q_bar_k*alpha_k(:,i)*(z_top^2 - z_bot^2);


    end



end






