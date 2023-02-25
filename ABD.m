function [A, B, D, alph, beta, del] = ABD(E1,E2,v12,G12, stack, t, stackmap)
    A = zeros(3,3);
    B = zeros(3,3);
    D = zeros(3,3);
    n = numel(stack);
    
    h = 0;
    for i = 1:n
        mat = stackmap(i);
        h = h + t(mat);
    end


    z_0 = -h/2;
    for i = 1:n
        mat = stackmap(i);

        z_top = z_0+t(mat)*i;
        z_bot = z_0+t(mat)*(i-1);
        Q_bar_k = Q_bar(E1(mat),E2(mat),v12(mat),G12(mat),stack(i));
        A = A + Q_bar_k*(z_top - z_bot);
        B = B + .5*Q_bar_k*(z_top^2 - z_bot^2);
        D = D + Q_bar_k*(z_top^3 - z_bot^3)/3;
    end
    invA = inv(A);
    
    del = inv(D-B*invA*B);
    beta = -invA*B*del;
    alph = invA+invA*B*del*invA;
    

end
