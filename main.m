%% 
clear;
clc;
%% system parameters
x_0 = [ 0.1 , 0 ;
        0.2 , 0 ;
        0.3 , 0 ;
        0.4 , 0 ;
        0.5 , 0 ; ];               % agent's initial state
yd = 0.1;                          % leader's state
A = [ 0 , 0 , 0 , 0 , 0;
      1 , 0 , 0 , 0 , 0 ;
      0 , 1.5 , 0 , 0 , 0 ;
      0 , 0 , 0 , 0 , 1 ;
      0 , 0 , 1 , 2 , 0 ;];        % adjacency matrix  aij
B = [ 1 ; 0 ; 0 ; 0.5 ; 0 ];       % Node connection weight b_i
C = [ 20 ; 95 ];                   % c_i1 c_i2
R = [ 0.4 ; 0.4 ];                 % r_i1 r_i2 
P = [ 0.61, 0.001; 
      0.61, 0.01 ; 
      0.61, 0.05 ; 
      0.61, 0.05 ; 
      0.62, 0.005; ];              % predefined accuracy p_i,j
D = cal_indegree_matrix(A);        % indegree matrix
%L = D - A;                        % laplacian matrix
ni = 2;                            % agent order
follower_num = 5;                  % number of agents
%% simulation parameters
T_start = 0;                                        % simulation starting time 
T_ending = 12;                                      % simulation ending time   
t = 0.001;                                          % simulation step time  
dt = 0.001;
N = (T_ending - T_start)/t;                         % simulation step  
%% simulation
x = zeros(5,6,N);                                   % x_i1  x_i2  x_dot_i1  x_dot_i2   yd   yd_dot
cauchy_i = 0.5 * ones(1,follower_num,N);            % the initial value of ¦Îi
cauchy_i_dot = zeros(1,follower_num,N);             % the value of ¦Îi_dot
cauchy_i_dot(1,:,1) = [1,1,1,1,1];
alpha_i1 = zeros(1,5,N);                            % ¦Ái1
alpha_i2 = zeros(1,5,N);                            % ui£¨¦Ái2£©
gamma_i = 1;                                        % the initial value of ¦£i   
leader_follow_error = zeros(1,5,N);                 % x_i1 - yd
x(:,:,1) = substitute_initial(x_0,yd,2,5,t);
e_i = zeros(follower_num,1,N);

for sim_period = 2 : 1 : N
    e_i(:,sim_period) = cal_syn_error(A,B,[x(1,1,sim_period-1);x(2,1,sim_period-1);
                             x(3,1,sim_period-1);x(4,1,sim_period-1);
                             x(5,1,sim_period-1)],x(1,5,sim_period-1),5);          % z_i,1
    for i = 1 : 1 : follower_num
        cauchy_i(1,i,sim_period) = cauchy_i(1,i,sim_period-1) + ( cauchy_i_dot(1,i,sim_period-1) * dt );
        z_i1 =  e_i(i,sim_period);
        state_1 = [x(1,1,sim_period-1),x(2,1,sim_period-1),x(3,1,sim_period-1),...
                   x(4,1,sim_period-1),x(5,1,sim_period-1),x(i,5,sim_period-1)];   % ¦Ö_i1
        fai_i1 = fai_ik(state_1,2,1,x(i,5,sim_period-1));                          % ¦Õi1
        w_i1 = sg_ik(z_i1,P(i,1),ni,1) * sqrt( (norm(fai_i1))^2 + R(1)^2 );
        alpha_i1(1,i,sim_period) = -( ( C(1)+0.5*(D(i,i)+B(i)) )*( (abs(z_i1)-P(i,1))^ni ) * sg_ik(z_i1,P(i,1),ni,1) +...
                         cauchy_i(1,i,sim_period)*w_i1 +(D(i,i)+B(i))*(P(i,2)+1)*sg_ik(z_i1,P(i,1),ni,1) )/ (D(i,i)+B(i));
        tao_i1 = (((abs(z_i1)-P(i,1)))^ni) * f_ik(z_i1,P(i,1)) * sg_ik(z_i1,P(i,1),ni,1) * w_i1;

        
        z_i2 =  x(i,2,sim_period-1) - alpha_i1(1,i,sim_period) ;
        state_2 = [x(1,1,sim_period-1),x(1,2,sim_period-1),x(2,1,sim_period-1),x(2,2,sim_period-1),...
                   x(3,1,sim_period-1),x(3,2,sim_period-1),x(4,1,sim_period-1),x(4,2,sim_period-1),...
                   x(5,1,sim_period-1),x(5,2,sim_period-1),x(i,5,sim_period-1)];   % ¦Ö_i2
        fai_i2 = fai_ik(state_2,2,2,x(i,5,sim_period-1));                          % ¦Õi2
        w_i2 = sg_ik(z_i2,P(i,2),ni,2) * sqrt( ((norm(fai_i2))^2) + R(2)^2 );
        tao_i2 = ( (abs(z_i2)-P(i,2))^(ni-1) ) * f_ik(z_i2,P(i,2)) * sg_ik(z_i2,P(i,2),ni,2) * w_i2 + tao_i1;
        alpha_i2(1,i,sim_period) = -( (C(2)+0.5)*(abs(z_i2)-P(i,2))*sg_ik(z_i2,P(i,2),ni,2) )-cauchy_i(1,i,sim_period)*w_i2+... 
                (-w_i1/(D(i,i)+B(i)))*tao_i2*gamma_i;
        
        cauchy_i_dot(1,i,sim_period) = gamma_i * tao_i2; 
        x(i,:,sim_period) = cal_follower_dynamics(x(i,:,sim_period-1),alpha_i2(1,i,sim_period),sim_period*t,dt);
        leader_follow_error(1,i,sim_period) = x(i,1,sim_period) - x(i,5,sim_period);
    end
end
%% Plotting
T = t * [1 : 1 : N];   % simulation time
y_d = zeros(1,N);      % trajectory of leader
for j = 1 : N
   y_d(j) = x(1,5,j); 
end

y_i = zeros(5,N);     % trajectory of agents y_i(i,:)-->i_th agent
for period = 1 : N
    for i = 1 : follower_num
           y_i(i,period) = x(i,1,period); 
    end
end
subplot(2,3,1);
plot(T,y_d,'r','linewidth',2);
hold on ;
for i = 1 : follower_num
   plot(T,y_i(i,:),'linewidth',1);
   hold on;
end
title('Output Trajectories');
xlabel('t(sec)');
ylabel('y_i');
legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');

% synchronization error
synchronization_error = zeros(follower_num,N);  
for c = 1 : N
    for d = 1 : follower_num
        synchronization_error(d,c) =e_i(d,1,c);
    end
end
subplot(2,3,2);
for j = 1 : follower_num
    plot(T,synchronization_error(j,:),'linewidth',1);
    hold on;
end
hold on;
title('synchronization error');
xlabel('t(sec)');
ylabel('syn error');
legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','NorthEast');

% error = x_i1 - yd
follow_error = zeros(follower_num,N);
for i = 1 : N
    for j = 1 : follower_num
        follow_error(j,i) = leader_follow_error(1,j,i);
    end
end
subplot(2,3,3);
for i = 1 : follower_num
   plot(T,follow_error(i,:),'linewidth',1);
   hold on;
end
hold on;
title('follow error(error = xi1 - yd)');
xlabel('t(sec)');
ylabel('follow error');
legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','NorthEast');

% ¦Îi
cauchy = zeros(follower_num,N);    
for a = 1 : N
    for b = 1 : follower_num
       cauchy(b,a) = cauchy_i(1,b,a); 
    end
end
subplot(2,3,4.5);
for j = 1 : follower_num
   plot(T,cauchy(j,:),'linewidth',1);
   hold on;
end
title('Parameter estimates ¦Îi');
xlabel('t(sec)');
ylabel('¦Î_i');
legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');

% ui£¨¦Ái2£©
ui = zeros(follower_num,N);     
for c = 1 : N
    for d = 1 : follower_num
        ui(d,c) = alpha_i2(1,d,c);
    end
end
subplot(2,3,5.5);
for j = 1 : follower_num
    plot(T,ui(j,:),'linewidth',1);
    hold on;
end
title('Controllers ui');
xlabel('t(sec)');
ylabel('u_i');
legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');




% %% plotting_test
% figure(1);
% for period = 1 : N
%     for i = 1 : follower_num
%            y_i(i,period) = x(i,1,period); 
%     end
% end
% plot(T,y_d,'r','linewidth',2);
% hold on ;
% for i = 1 : follower_num
%    plot(T,y_i(i,:),'linewidth',1);
%    hold on;
% end
% title('Output Trajectories');
% xlabel('t(sec)');
% ylabel('y_i');
% legend('Leader','Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');
% 
% figure(2);
% for c = 1 : N
%     for d = 1 : follower_num
%         synchronization_error(d,c) =e_i(d,1,c);
%     end
% end
% for j = 1 : follower_num
%     plot(T,synchronization_error(j,:),'linewidth',1);
%     hold on;
% end
% hold on;
% title('synchronization error');
% xlabel('t(sec)');
% ylabel('syn error');
% legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','NorthEast');
% % 
% figure(3);
% for i = 1 : N
%     for j = 1 : follower_num
%         follow_error(j,i) = leader_follow_error(1,j,i);
%     end
% end
% for i = 1 : follower_num
%    plot(T,follow_error(i,:),'linewidth',1);
%    hold on;
% end
% hold on;
% title('follow error(error = x_i - x_0)');
% xlabel('t(sec)');
% ylabel('follow error');
% legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','NorthEast');
% % 
% figure(4);
% for a = 1 : N
%     for b = 1 : follower_num
%        cauchy(b,a) = cauchy_i(1,b,a); 
%     end
% end
% for j = 1 : follower_num
%    plot(T,cauchy(j,:),'linewidth',1);
%    hold on;
% end
% title('Parameter estimates ¦Îi');
% xlabel('t(sec)');
% ylabel('¦Î_i');
% legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');
% % 
% figure(5);
% for c = 1 : N
%     for d = 1 : follower_num
%         ui(d,c) = alpha_i2(1,d,c);
%     end
% end
% for j = 1 : follower_num
%     plot(T,ui(j,:),'linewidth',1);
%     hold on;
% end
% title('Controllers ui');
% xlabel('t(sec)');
% ylabel('u_i');
% legend('Follower 1','Follower 2','Follower 3','Follower 4','Follower 5','Location','SouthEast');