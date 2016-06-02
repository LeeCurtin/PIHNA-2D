function [R, J] = vasc_reaction(Q, ...
                                 rho, K, ...
                                 beta, gamma, alpha_h, ...
                                 delta_c, delta_h, K_M, q, lambda_a, omega, ...
                                 mu_v, ...
                                 alpha_n,...
                                 Ktrans_low,Ktrans_hi,dl,v0,a0)
% Given Q=[c; h; a; v; n; l], this function calculates the reaction
% R=R(Q), and the Jacobian J=R'(Q). (l added 8/1/2011 ahd)

%Ktrans = permeability coefficient
%dl = liquid removal rate

    c = Q(1);
    h = Q(2);
    a = Q(3);
    v = Q(4);
    n = Q(5);


%     T = (c+h+v+n)/K;
    T = (c+h+v+n);
   
%New Vascular Efficiency term 
    a_min = 1; %this has to be less than sqrt(2) to avoid a_fun <= 0
    a_fun = ((a - a0)/2 - a_min)*(a - a0) + 1; 
    %The above quadratic shrinks (a-a0<a_min) and stretches tanh (a-a0>a_min), Note:tanh curve is unchanged at a=a0+2*a_min
    ch_shift = 0.1; %The factor that determines how much stress the tumour puts on the vasculature by shifting v0 to the right
    save_space = 100*((v-v0*(1 + ch_shift*(c+h)))*(a_fun)); %Just to save space in typing this over and over again in Jacobian
    F = (tanh(save_space)+1)/2; %New Vascular Efficiency term!
    
    
%Calculate Ktrans for current values of a
    Ktrans = 0;
    %if(c>0.1) Ktrans = Ktrans_low; end
    
    Ktrans = Ktrans + a/(K_M+a)*Ktrans_hi;
    

% Production terms
    production_c = rho * c * (1-T);
    production_a = delta_c*c + delta_h*h;
    production_v = mu_v * a/ (1 + a) * v * (1-T);

% Conversion terms
    c_to_h = beta * c * (1 - F);
    c_to_n = alpha_n * (n/K) * c;

    h_to_c = gamma * h * F;
    h_to_n = alpha_h * h * (1 - F) + alpha_n*(n/K)*h;

    v_to_n = alpha_n*(n/K)*v; %alpha_n * n * v;

% Loss terms
    loss_a = lambda_a*a + q * production_v + omega*a*v;
   


% The reaction vector, R(Q)
    R = zeros(5,1);
    R(1) = production_c - c_to_h + h_to_c - c_to_n;
    R(2) = c_to_h - h_to_c - h_to_n;
    R(3) = production_a - loss_a;
    R(4) = production_v - v_to_n;
    R(5) = c_to_n + h_to_n + v_to_n;


% The Jacobian matrix, J = R'(Q)
    J = zeros(5);

% First row
    J(1,1) = rho*(1-T) - rho*c/K - beta*(1-F) - 50*beta*v0*c*a_fun*ch_shift*sech(save_space).^2 ... 
             - 50*gamma*v0*a_fun*h*ch_shift*sech(save_space).^2 - alpha_n*n/K;
    J(1,2) = -rho*c/K - 50*ch_shift*beta*c*v0*a_fun*sech(save_space).^2 + gamma*F - 50*gamma*h*v0*a_fun*ch_shift*sech(save_space).^2;
    J(1,3) = 50*gamma*h*c*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2 + 50*beta*c*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
    J(1,4) = -rho*c/K + 50*beta*c*a_fun*sech(save_space).^2 + 50*gamma*h*a_fun*sech(save_space).^2;
    J(1,5) = -(rho+alpha_n)*c/K;
    
% Second row
    J(2,1) = beta*(1-F) + 50*beta*v0*c*a_fun*ch_shift*sech(save_space).^2 + 50*gamma*v0*a_fun*h*ch_shift*sech(save_space).^2 ...
             - 50*alpha_h*h*v0*a_fun*ch_shift*sech(save_space).^2;
    J(2,2) = + 50*beta*c*v0*a_fun*ch_shift*sech(save_space).^2 - gamma*F + 50*gamma*v0*a_fun*h*ch_shift*sech(save_space).^2 ...
             - alpha_h*(1-F) - 50*alpha_h*h*v0*a_fun*ch_shift*sech(save_space).^2 - alpha_n*n/K;
    J(2,3) = - 50*gamma*h*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2 ...
             - 50*beta*c*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2 ...
             + 50*alpha_h*h*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
    J(2,4) = - 50*beta*c*a_fun*sech(save_space).^2 ...
             - 50*gamma*h*a_fun*sech(save_space).^2 + 50*alpha_h*h*a_fun*sech(save_space).^2;
    J(2,5) = -alpha_n*h/K;
    
% Third row
    J(3,1) = delta_c + q*mu_v*(a/(1+a))*v/K;
    J(3,2) = delta_h + q*mu_v*(a/(1+a))*v/K;
    J(3,3) = -lambda_a - q*mu_v*v*(1-T)/(1+a) + q*mu_v*a*v*(1-T)/(1+a)^2 - omega*v;
    J(3,4) = - q*mu_v*a/(1+a)*(1-T) + q*mu_v*a/(1+a)*(v/K) - omega*a;
    J(3,5) = + q*mu_v*a/(1+a)*(v/K);
    
% Fourth row
    J(4,1) = -mu_v*a/(1+a)*v/K;
    J(4,2) = -mu_v*a/(1+a)*v/K;
    J(4,3) = mu_v*v*(1-T)/(1+a) - mu_v*a*v*(1-T)/(1+a)^2;
    J(4,4) = mu_v*a*(1-T)/(1+a) - mu_v*a*v/(1+a) - alpha_n*n/K;
    J(4,5) = -mu_v*a/(1+a)*v/K - alpha_n*v/K;

% Fifth row
    J(5,1) = 50*alpha_h*h*v0*ch_shift*a_fun*sech(save_space).^2 + alpha_n*n/K;
    J(5,2) = alpha_h*(1-F) + 50*alpha_h*h*v0*a_fun*ch_shift*sech(save_space).^2 + alpha_n*n/K;
    J(5,3) = -50*alpha_h*h*(v-v0*(1+ch_shift*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
    J(5,4) = -50*alpha_h*h*a_fun*ch_shift*sech(save_space).^2 + alpha_n*n/K;
    J(5,5) = alpha_n*(c+h+v)/K;
% 
% function [R, J] = vasc_reaction(Q, ...
%                                  rho, K, ...
%                                  beta, gamma, alpha_h, ...
%                                  delta_c, delta_h, K_M, q, lambda_a, omega, ...
%                                  mu_v, ...
%                                  alpha_n,...
%                                  Ktrans_low,Ktrans_hi,dl,v0,a0)
% % Given Q=[c; h; a; v; n; l], this function calculates the reaction
% % R=R(Q), and the Jacobian J=R'(Q). (l added 8/1/2011 ahd)
% 
% %Ktrans = permeability coefficient
% %dl = liquid removal rate
% 
%     c = Q(1);
%     h = Q(2);
%     a = Q(3);
%     v = Q(4);
%     n = Q(5);
% 
% 
%     %T = (c+h+v+n)/K;
%     T = (c+h+v+n);
%    
% %New Vascular Efficieny term 
%     a_min = 1; %this has to be less than sqrt(2) to avoid a_fun <= 0
% 
%     a_fun = ((a - a0)/2 - a_min)*(a - a0) + 1; 
%     %The above quadratic shrinks (a-a0<a_min) and stretches tanh (a-a0>a_min), Note:tanh curve is unchanged at a=a0+2*a_min
%     save_space = 100*((v-v0*(1 + 0.01*(c+h)))*(a_fun)); %Just to save space in typing this over and over again in Jacobian
%     F = (tanh(save_space)+1)/2; %New Vascular Efficiency term!
%     
% %Calculate Ktrans for current values of a
%     Ktrans = 0;
%     %if(c>0.1) Ktrans = Ktrans_low; end
%     
%     Ktrans = Ktrans + a/(K_M+a)*Ktrans_hi;
%     
% 
% % Production terms
%     production_c = rho * c * (1-T);
%     production_a = delta_c*c + delta_h*h;
%     production_v = mu_v * a/ (K_M + a) * v * (1-T);
% 
% % Conversion terms
%     c_to_h = beta * c * (1 - F);
%     c_to_n = alpha_n*(n/K)*c;
% 
%     h_to_c = gamma * h * F;
%     h_to_n = alpha_h * h * (1 - F) + alpha_n*(n/K)*h;
% 
%     v_to_n = alpha_n*(n/K)*v; %alpha_n * n * v;
% 
% % Loss terms
%     loss_a = lambda_a*a + q * production_v + omega*a*v;
%    
% 
% 
% % The reaction vector, R(Q)
%     R = zeros(5,1);
%     R(1) = production_c - c_to_h + h_to_c - c_to_n;
%     R(2) = c_to_h - h_to_c - h_to_n;
%     R(3) = production_a - loss_a;
%     R(4) = production_v - v_to_n;
%     R(5) = c_to_n + h_to_n + v_to_n;
% 
% 
% % The Jacobian matrix, J = R'(Q)
%     J = zeros(5);
% 
% % First row
%     J(1,1) = rho*(1-T) - rho*c/K - beta*(1-F) - 0.5*beta*v0*c*a_fun*sech(save_space).^2 ... 
%              - 0.5*gamma*v0*a_fun*h*sech(save_space).^2 - alpha_n*n/K;
%     J(1,2) = -rho*c/K - 0.5*beta*c*v0*a_fun*sech(save_space).^2 + gamma*F - 0.5*gamma*h*v0*c*a_fun*sech(save_space).^2;
%     J(1,3) = 50*beta*c*(v-v0*(1+0.01*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
%     J(1,4) = -rho*c/K + 50*beta*c*a_fun*sech(save_space).^2 + 50*gamma*h*a_fun*sech(save_space).^2;
%     J(1,5) = -(rho+alpha_n)*c/K;
%     
% % Second row
%     J(2,1) = beta*(1-F) + 0.5*beta*v0*c*a_fun*sech(save_space).^2 + 0.5*gamma*v0*a_fun*h*sech(save_space).^2 ...
%              - alpha_h*h*v/(c+h+v)^2;
%     J(2,2) = + 0.5*beta*c*v0*a_fun*sech(save_space).^2 - gamma*F + 0.5*gamma*v0*a_fun*h*sech(save_space).^2 ...
%              - alpha_h*(1-T) + alpha_h*h*(1/K) - alpha_n*n;
%     J(2,3) = - 50*beta*c*(v-v0*(1+0.01*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
%     J(2,4) = - 50*beta*c*a_fun*sech(save_space).^2 ...
%              - 50*gamma*h*a_fun*sech(save_space).^2 + alpha_h*h*(1/K);
%     J(2,5) = -alpha_n*h;
%     
% % Third row
% %     J(3,1) = delta_c - q*mu_v*(a/(K_M+a))*v/K;
% %     J(3,2) = delta_h - q*mu_v*(a/(K_M+a))*v/K;
% %     J(3,3) = -lambda_a + q*mu_v*v*(1-T)/(K_M+a) - q*mu_v*a*v*(1-T)/(K_M+a)^2 - omega*v;
% %     J(3,4) = q*mu_v*a/(K_M+a)*(1-T) - q*mu_v*a/(K_M+a)*v/K - omega*a;
% %     J(3,5) = -q*mu_v*a/(K_M+a)*v/K;
% 
%     J(3,1) = delta_c + q*mu_v*(a/(K_M+a))*v/K;
%     J(3,2) = delta_h + q*mu_v*(a/(K_M+a))*v/K;
%     J(3,3) = -lambda_a - q*mu_v*v*(1-T)/(K_M+a) + q*mu_v*a*v*(1-T)/(K_M+a)^2 - omega*v;
%     J(3,4) = - q*mu_v*a/(K_M+a)*(1-T) + q*mu_v*a/(K_M+a)*v/K - omega*a;
%     J(3,5) = + q*mu_v*a/(K_M+a)*v/K;
%     
% % Fourth row
%     J(4,1) = -mu_v*a/(K_M+a)*v/K;
%     J(4,2) = -mu_v*a/(K_M+a)*v/K;
%     J(4,3) = mu_v*v*(1-T)/(K_M+a) - mu_v*a*v*(1-T)/(K_M+a)^2;
%     J(4,4) = mu_v*a*(1-T)/(K_M+a) - mu_v*a/(K_M+a)*v/K - alpha_n*n;
%     J(4,5) = -mu_v*a/(K_M+a)*v/K - alpha_n*v;
% 
% % Fifth row
%     J(5,1) = 0.5*alpha_h*v0*a_fun*sech(save_space).^2 + alpha_n*n;
%     J(5,2) = alpha_h*(1-F) + 0.5*alpha_h*h*v0*a_fun*sech(save_space).^2 + alpha_n*n;
%     J(5,3) = -50*alpha_h*(v-v0*(0.01*(c+h)))*(a-a_min-a0)*sech(save_space).^2;
%     J(5,4) = -50*alpha_h*a_fun*sech(save_space).^2 + alpha_n*n;
%     J(5,5) = alpha_n*(c+h+v);
