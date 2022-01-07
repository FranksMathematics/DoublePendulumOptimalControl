function robot_control()

    % init
    clc;
    clear all;



    global N
    global beta
    global dest_x
    global dest_y



    beta = 1000000;

    dest_x = 2.0;
    dest_y = 2.0;


    global m1;
    global m2;
    global l1;
    global l2;
    global g;

    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 0.5;
    g = 9.81;

    % render video
    doVideo = 1;



    N = 40;

    % block structure of solution vector
    % theta_1
    % theta_2
    % omega_1
    % omega_3
    % u_1
    % u_2
    % repeat N times





    % initial data
    x0 = zeros(6*N,1);


    % Solve optimization problem
    qopt = optimset('MaxFunEvals', 500000, ...
                    'maxIter', 150, ...
                    'Display', 'iter');


    %q = fmincon(@fmin,x0,[],[],[],[],lb,ub,@mycon, qopt);
    q = fmincon(@fmin,x0,[],[],[],[],[],[],@mycon, qopt);


    % resolve data
    theta_1 = zeros(N,1);
    theta_2 = zeros(N,1);
    omega_1 = zeros(N,1);
    omega_2 = zeros(N,1);
    u_1 = zeros(N,1);
    u_2 = zeros(N,1);
    l = 1;
    for i=1:N
        theta_1(i) = q(l);
        theta_2(i) = q(l+1);
        omega_1(i) = q(l+2);
        omega_2(i) = q(l+3);
        u_1(i) = q(l+4);
        u_2(i) = q(l+5);

        l = l+6;
    end



    h = 1/(N-1);
    t = 0:h:1;
    figure(1);
    subplot(6,1,1) ; plot (t, theta_1) , grid on , title  ('Theta 1')
    subplot(6,1,2) ; plot (t, theta_2) , grid on , title  ('Theta 2')
    subplot(6,1,3) ; plot (t, omega_1) , grid on , title  ('Omega 1')
    subplot(6,1,4) ; plot (t, omega_2) , grid on , title  ('Omega 2')
    subplot(6,1,5) ; plot (t, u_1) , grid on , title  ('Control Theta 1')
    subplot(6,1,6) ; plot (t, u_2) , grid on , title  ('Control Theta 2')






    if doVideo == 1
        videofile = VideoWriter('robot_control.avi');
        % extend video to 10 seconds
        videofile.FrameRate = round(1/(10*h));
        open(videofile);
    end



    for i=1:N
        x0 = 0;
        y0 = 0;

        x1 = l1*sin( theta_1(i) );
        y1 = - l1*cos( theta_1(i) );

        x2 = x1 + l2*sin( theta_2(i) );
        y2 = y1 - l2*cos( theta_2(i) );


        f2 = figure(2);
        clf(f2);
        set(gcf, 'position', [0 0 1000 1000]);
        %clf(1);
        hold on
        plot(x0,y0,'o','markerfacecolor','b');
        hold on
        plot(x1,y1,'o','markerfacecolor','b');
        hold on
        plot(x2,y2,'o','markerfacecolor','b');
        hold on
        plot( [x0,x1], [y0,y1] );
        hold on
        plot( [x1,x2], [y1,y2] );
        hold off
        axis([-2 2 -2 2]);

        if doVideo == 1
            frame = getframe(gcf);
            writeVideo(videofile, frame);
        end
    end


    if doVideo == 1
        close(videofile);
    end


    end





function [c,ceq] = mycon(x)
    global N;

    % no inequalities
    c = [];



    h =  1  /(N-1);

    % use equality constraints

    % initial constraint
    % pendulum is at rest for t=0   (ceq = 0)
    ceq(1) = x(1);  % theta_1 = 0
    ceq(2) = x(2);  % theta_2 = 0
    ceq(3) = x(3);  % omega_1 = 0
    ceq(4) = x(4);  % omega_2 = 0



    % inegral approximation y_{k+1} - y_k - 0.5*h*(f_{k+1} + f_k) = 0
    j = 5; % number of equality constraint
    l = 1; % position inside vector x
    for i=1:N-1

        % extract y_i, u_i and y_{i+1}, u_{i+1} for easier calculation
        y_i = x(l: l+3);
        u_1_i = x(l+4);
        u_2_i = x(l+5);

        y_ii = x(l+6: l+9);
        u_1_ii = x(l+10);
        u_2_ii = x(l+11);

        % get RHS of ODE
        f_i  = ODE_RHS(y_i,  u_1_i,  u_2_i );
        f_ii = ODE_RHS(y_ii, u_1_ii, u_2_ii);

        % add equations
        ceq(j:j+3) = y_ii - y_i - 0.5*h*(  f_ii + f_i  );

        % increase numbers if equations and counter in vector
        j = j+4;  
        l = l+6;
    end


    % pendulum is at rest for t=T
    ceq(j+1) = x(6*N-3);  % omega_1 = 0
    ceq(j+2) = x(6*N-2);  % omega_2 = 0

end




function res = fmin(x)

    global N
    global beta
    global dest_x
    global dest_y

    h =1/(N-1);


    % Simple cost-term for the control

    l = 1; % position inside vector x
    res1 = 0;
    for i=1:N-1

        % extract y_i, u_i and y_{i+1}, u_{i+1} for easier calculation
        %y_i = x(l: l+4);
        u_1_i = x(l+4);
        u_2_i = x(l+5);
        u_1_ii = x(l+10);
        u_2_ii = x(l+11);

        res1 = res1 + 0.5*h*(u_1_i^2 + u_1_ii^2 + u_2_i^2 + u_2_ii^2);  

        % increase numbers if equations and counter in vector
        l = l+6;
    end



    % penalty term for position violation at the end
    y_end = x(6*(N-1)+1:6*N);
    [xp,yp] = AngularToCartesian(y_end);
    res2 = beta*(  (xp - dest_x)^2 + (yp - dest_y)^2   );


    res = res1 + res2;

end






function res = alpha1(t1, t2)
    global m1;
    global m2;
    global l1;
    global l2;
    global g;

    res = (l2/l1)*( m2/(m1+m2)  )*cos(t1-t2);
end
function res = alpha2(t1, t2)
    global m1;
    global m2;
    global l1;
    global l2;
    global g;

    res = (l1/l2)*cos(t1-t2);

end
function res = f1(t1, t2, w1, w2)
    global m1;
    global m2;
    global l1;
    global l2;
    global g;

    res = (-l2/l1)*(m2/(m1+m2))*w2*w2*sin(t1-t2) - (g/l1)*sin(t1);

end
function res = f2(t1, t2, w1, w2)
    global m1;
    global m2;
    global l1;
    global l2;
    global g;

    res = (l1/l2)*w1*w1*sin(t1-t2) - (g/l2)*sin(t2);

end
function res = g1(y)
 
    t1 = y(1);
    t2 = y(2);
    w1 = y(3);
    w2 = y(4);

    res = ( f1(t1,t2,w1,w2) - alpha1(t1,t2)*f2(t1,t2,w1,w2)  ) / (1 - alpha1(t1,t2)*alpha2(t1,t2));

end
function res = g2(y)

    t1 = y(1);
    t2 = y(2);
    w1 = y(3);
    w2 = y(4);

    res = ( -alpha2(t1,t2)*f1(t1,t2,w1,w2) + f2(t1,t2,w1,w2)  ) / (1 - alpha1(t1,t2)*alpha2(t1,t2));

end
function res = ODE_RHS(y,u1,u2)

    res = [ y(3); y(4); g1(y) + u1; g2(y) + u2];

end



function  [px,py] = AngularToCartesian(y)

   global l1;
   global l2;

   x1 = l1*sin( y(1) );
   y1 = - l1*cos( y(1) );
   
   px = x1 + l2*sin( y(2) );
   py = y1 - l2*cos( y(2) );

end





