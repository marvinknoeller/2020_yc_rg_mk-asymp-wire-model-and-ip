function [ XX,errN,errD,errU,errG,j ] = inv_test_maxwell_spline(Var,X,FF )
% inv_test returns all iterates of the newton-like algorithm and several
% error terms which appear in the function we aim to minimize.
%
% Input:
% Var - a structure which contains all the coefficients and parameters
% X - the initial guess for the curve
% FF - The far field we computed by means of Bem++
%
% Output:
% XX - a structure which contains all the iterates of the algorithm
% errN - the error of the far field
% errD - the error of the first constraint
% errU - the error of the second constraint
% errG - the full error
% j - it took j-1 steps to get to the final iterate
mu_rel = Var.mu_rel;
eps_rel = Var.eps_rel;
N=10;
N2=2*N;
th=linspace(pi/N,pi/N*(N-1),N-1);
Theta=ones(N2,1)*th;
th2=linspace(0,2*pi-2*pi/N2,N2)';
Phi=th2*ones(1,N-1);
% Cartesian coordinates of the sampling points and their number
Z = [sin(Theta(:)').*cos(Phi(:)'); sin(Theta(:)').*sin(Phi(:)'); cos(Theta(:)')];
NN = size(Z,2);
%% store the given parameters in local variables
%number of segments
n=Var.n;
%the far field computed with Bem++
E_infty_BEM=FF;

%% further setup
% weights of different measurements with respect to the L2 boundary norm
w = pi/N*sqrt(sin(Theta(:)));
% Number of integration points in each sub line segment
M=11;
% relative noise level. e.g. 0.1 corr. to 10% noise
if Var.noise==1
    err = (1*rand(N2*(N-1),3)-.5+1i*(1*rand(N2*(N-1),3)-.5));
    err = err/err_on_ff(w.',err.')*Var.noise_level*err_on_ff(w.',E_infty_BEM.');
    E_infty_BEM=E_infty_BEM+err;
    err_on_ff(w.',err.')
end
%%
norm_E_infty_BEM = err_on_ff(w.',E_infty_BEM.');
%%
%% Choose the regularization parameters (note that the n in front of Psi_1 and Psi_2 is included here)
% there could possibly exist different reg. par. for noise, however, we
% took always the same coefficients
%%
lambda1 = 2e-1;
lambda2 = 9e-1;
%% The first iterate
XX{1} = X;
%% The iterative reconstruction algorithm
%counts the iterates
j = 1;
% check index for finding out if we are at a local minimum
counter = 1;
% these parameters corr. to the extended algorithm.
jump_over = 1;  % just for skipping the first if condition
while j < 500 %the condition counter==0 has been replaced
    if counter >= 1 && jump_over == 0   % we moved
        %move
    elseif counter == 0 % we did not move
        if errN(j-2)<errD(j-2) && errN(j-2)<errU(j-2)
            lambda1 = lambda1*0.5;
            lambda2 = lambda2*0.5;
        elseif errN(j-2)<errD(j-2) && errN(j-2)>=errU(j-2)
            lambda1 = lambda1*0.5;
        elseif errN(j-2)>=errD(j-2) && errN(j-2)<errU(j-2)
            lambda2 = lambda2*0.5;
        else
            break
        end
    end
    jump_over = 0;
    %% Computation of the far fields, penalty terms and derivatives
    num_x = n+1;
    [X_in_between,ww,Pol_1,Pol_2] = SetupFarField(X,M,num_x,mu_rel,eps_rel);
    % far field and the Jacobian
    E_infty_Pert = FarField_Pert_Maxwell_E_spline(Var,X_in_between,Pol_1,Pol_2,ww,Z);
    [J,Jp1new,Jp2new] = Jacobian_FarField_maxwell_spline(Var,X,Z,M);
    % weighting corresponding to the L2-norm on the unit sphere
    gw = 1/norm_E_infty_BEM*(E_infty_Pert - E_infty_BEM).*w;
    gw = gw.';
    gw_vec = reshape(gw,3*NN,1);
    J = 1/norm_E_infty_BEM*w.*J;
    J_shape = reshape(permute(J,[2 1 3]),[3*NN,3*(n+1)]);
    % split it in real and imag part
    J1 = real(J_shape);
    J2 = imag(J_shape);
    gw1 = real(gw_vec);
    gw2 = imag(gw_vec);
    % first penalty term penalizes strong curvature
    p11 = Curv_Pen(X,M,num_x);
    % vector whose norm is the penalty term Psi_1
    p1 = p11(:);
    % Jacobian of the penalty term with respect to all points
    Jp1 = reshape(Jp1new,3*((num_x-1)*1*(M-1)+1),3*num_x);
    % second regularization term Psi_2 penalizes varying sublengths of
    % segments
    p22 = Seg_Pen_New(X,M,num_x);
    % vector whose norm is the penalty term Psi_2
    p2 = p22(:);
    % corresponding Jacobian
    Jp2 = Jp2new;
    %% Computation of the Gauss-Newton direction
    % Gauss-Newton iterate
    h = (J1'*J1 + J2'*J2  + lambda1^2*(Jp1'*Jp1) + lambda2^2*(Jp2'*Jp2))\...
        ([J1'  J2'  lambda1*Jp1' lambda2*Jp2']*[-gw1; -gw2; -lambda1*p1; -lambda2*p2]);
    % suggested movement of the end points
    H = reshape(h, 3, n+1);
    %% For plotting different penalty terms on logarithmic scale (for j=0 some terms equal zero)
    if j~=1
        errN(j-1)=sum(sum((abs(gw)).^2,1));
        errD(j-1)=lambda1^2*norm(p1)^2;
        errU(j-1)=lambda2^2*norm(p2)^2;
        errG(j-1)=sum(sum((abs(gw)).^2,1)) + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
    else
        errN(j)=sum(sum((abs(gw)).^2,1));
        errD(j)=lambda1^2*norm(p1)^2;
        errU(j)=lambda2^2*norm(p2)^2;
        errG(j)=sum(sum((abs(gw)).^2,1)) + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
    end
    %% Finding the minimum in the direction hh using golden ratio search
    golden = 2/(1 + sqrt(5));
    smax = 1;
    X1 = X;
    X4 = X + smax*H;
    X2 = X4 - golden*(X4 - X1);
    X3 = X1 + golden*(X4 - X1);
    [X2_in_between,ww2,Pol2_1,Pol2_2] = SetupFarField(X2,M,n,mu_rel,eps_rel);
    [X3_in_between,ww3,Pol3_1,Pol3_2] = SetupFarField(X3,M,n,mu_rel,eps_rel);
    E_infty_2 = FarField_Pert_Maxwell_E_spline(Var,X2_in_between,Pol2_1,Pol2_2,ww2,Z);
    E_infty_3 = FarField_Pert_Maxwell_E_spline(Var,X3_in_between,Pol3_1,Pol3_2,ww3,Z);
    % first penalty term
    p21 = Curv_Pen(X2,M,num_x);
    p21 = p21(:);
    % second penalty term
    p22 = Seg_Pen_New(X2,M,num_x);
    p22 = p22(:);
    % value of the functional to be minimized
    funct2 = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_2.')^2 + lambda1^2*norm(p21)^2 + lambda2^2*norm(p22)^2;
    % first penalty term
    p31 = Curv_Pen(X3,M,num_x);
    p31 = p31(:);
    % second penalty term
    p32 = Seg_Pen_New(X3,M,num_x);
    p32 = p32(:);
    % value of the functional to be minimized
    funct3 = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_3.')^2 + lambda1^2*norm(p31)^2 + lambda2^2*norm(p32)^2;
    counter = 0;
    % count how many times the farther point is chosen
    % ten divisions by the golden ratio search
    for k=1:10
        if funct2 <= funct3
            X4 = X3;
            X3 = X2;
            X2 = X4 - golden*(X4 - X1);
            Xnew = X2;
            funct3 = funct2;
            ind = 2;
        else
            X1 = X2;
            X2 = X3;
            X3 = X1 + golden*(X4 - X1);
            Xnew = X3;
            funct2 = funct3;
            ind = 3;
            counter = counter + 1;
        end
        [Xnew_in_between,wwnew,Polnew_1,Polnew_2] = SetupFarField(Xnew,M,n,mu_rel,eps_rel);
        E_infty_Pert = FarField_Pert_Maxwell_E_spline(Var,Xnew_in_between,Polnew_1,Polnew_2,wwnew,Z);
        % first penalty term
        p11 = Curv_Pen(Xnew,M,num_x);
        p1 = p11(:);
        % second penalty term
        p22 = Seg_Pen_New(Xnew,M,num_x);
        p2 = p22(:);
        % value of the functional to be minimized
        funkt = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_Pert.')^2 + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
        if ind==2
            funct2 = funkt;
        else
            funct3 = funkt;
        end
    end
    % conservative choices for the parameter values and sum of the movement
    total_movement = sum(sqrt(sum((X - X1).^2)));
    X = X1;
    % store the values
    XX{j+1} = X;
    %% Monitoring
    % monitor the progress
    X1 = splinepoints(X,10);
    plot3(X1(1,:),X1(2,:),X1(3,:), 'r-', 'LineWidth', 1)
    hold on
    plot3(X(1,:),X(2,:),X(3,:), 'r*', 'LineWidth', 1)
    hold off
    drawnow
    display([j, total_movement])
    j = j+1;
    
end

end