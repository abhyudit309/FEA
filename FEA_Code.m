%Kishore Ram Sathia - ME18B085
%Abhyudit Singh Manhas - ME18B088
%AS5850 - Jan-May 2021 Project
%May 25, 2021
%2D C0-quadratic quadrilateral isoparametric element FEM Solution
clear;
clc;

s_EF                = 16.8E6; %sigma_nought along top edge
tau_EF              = 0;
s_AF                = 0;
tau_AF              = 0;
s_ED                = 0;
tau_ED              = 0;
s_AB                = 0;
tau_AB              = 0;
s_CD                = 0;
tau_CD              = 0;
tau_BC              = 0;

Lx                  = 1; %SI Units: metres
Ly                  = 1;
nodeB               = Lx/4;
nodeC               = 3*Lx/4;

n_AF                = 28;            %initial no of elements in the vertical direction (ignoring midpoints of edges)
n_AB                = (n_AF - 4)/4;  %initial no of elements in the horizontal direction between A and B (not including A, B, not including midpoints)
n_BC                = (n_AF - 4)/2;  %initial no of elements in the horizontal direction between B and C (not including B, C, not including midpoints)
n_CD                = (n_AF - 4)/4;  %initial no of elements in the horizontal direction between C and D (not including C, D, not including midpoints)
n = n_AF;

tolerance           = 0.018; %This allows a convergence with n=52. For any larger n, the program takes very long to solve, may even cause MATLAB to hang
iter_step           = 4;  %change this from 12 to a higher multiple of 4 for faster convergence
max_iters           = 10;
isConverged         = false;

%check for convergence by comparing derivative values of one iteration with
%next at 7x7 uniformly spaced points
xnodes_convcheck            = linspace(Lx/10, 9*Lx/10, 7);    %10% within the domain, so that not too close to the boundaries
ynodes_convcheck            = linspace(Ly/10, 9*Ly/10, 7);

prev_iteration_deriv_sx         = zeros(7,7); 
current_iteration_deriv_sx      = zeros(7,7);
prev_iteration_deriv_sy         = zeros(7,7); 
current_iteration_deriv_sy      = zeros(7,7);
prev_iteration_deriv_tauxy      = zeros(7,7); 
current_iteration_deriv_tauxy   = zeros(7,7);

%3 point Gauss-Legendre quadrature
GP                          = [-sqrt(3)/sqrt(5) 0 sqrt(3)/sqrt(5)]';
w                           = [5/9 8/9 5/9]';

%Material properties
E                           = 2E11; %Young's Modulus (Pa)
nu                          = 0.3;  %Poisson's ratio

C                           = E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 1-nu]; %3x3 matrix for plane stress

% %Calculate Smoothing matrix:
xi  = [GP; GP; GP];
eta = [];
for p = 1:length(GP)
    for q = 1:length(GP)
        eta = [eta GP(p)];
    end
end

Q = [];
for i = 1:length(eta)
    Q   = [Q [1 eta(i) xi(i)]'];
end
P = zeros(3);
for i = 1: length(xi)
   lP   = [1 xi(i) eta(i);xi(i) xi(i)^2 xi(i)*eta(i);eta(i) xi(i)*eta(i) xi(i)^2];
   P    = P + lP;
end
S = (P^-1)*Q; %The smoothing matrix

E_2     = [1 -1 -1;1 1 -1; 1 1 1; 1 -1 1; 1  0 -1; 1 1 0; 1 0 1; 1 -1 0]; %to avoid confusion with E
%Transformation matrix
TR      = E_2*S;

n_iters = 1;


while  ~isConverged
    x_AB        = linspace(0, nodeB, n_AB+2);
    x_BC        = linspace(nodeB, nodeC, n_BC+2);
    x_CD        = linspace(nodeC, Lx, n_CD+2);
    x           = [x_AB(1:end-1) x_BC(1:end-1) x_CD];
    y           = linspace(0, Ly, n_AF);
    
    x_mid       = conv(x, [0.5 0.5], 'valid');
    x_withmid   = [x;[x_mid 0]];
    x_withmid   = x_withmid(:)';
    x_withmid   = x_withmid(1:end-1);
    
    K           = zeros(2*(3*n^2 - 2*n));    %Assembled K matrix  n^2 + 2*n*(n-1)
    
    f           = zeros(2*(3*n^2 - 2*n), 1); %Assembled f matrix
    
    for i = 1:(length(y) - 1)
        for j = 1:(length(x) - 1) 
            
           x_ele    = [x(j); x(j+1); x(j+1); x(j)]; %corner points of element counter-clockwise from bottom left
           y_ele    = [y(i); y(i); y(i+1); y(i+1)]; %corner points of element counter-clockwise from bottom left
           K_ele    = zeros(16);    %for a quadratic element, local K is (2*8)x(2*8)
           f_ele    = zeros(16, 1);  
           
           hlen_ele = x_ele(3) - x_ele(1);
           vlen_ele = y_ele(3) - y_ele(1);
           
           
           %iterate over Gauss Points:
           for p = 1:length(GP)
               for q = 1:length(GP)
                   B        = [];
                   xi       = GP(p);
                   eta      = GP(q);
                   J_ele    = Jacobian(x_ele, y_ele, xi, eta);
                   for k = 1:8 %8 trial functions
                    B       = [B diag(J_ele^-1*delPhiK_delIso(k, xi, eta))];
                   end
                   B3       = [B(2,2:2:end);B(1,1:2:end)];
                   B(3,:)   = (B3(:))';
                   K_ele    = K_ele + B'*C*B*det(J_ele)*w(p)*w(q);
               end
           end
           
           %connectivity matrix:
           nodeNumbers_ele  = [(i-1)*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 2; i*(3*n_AF-1)+(2*j-1)+2; i*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 1; (i-1)*(3*n_AF-1)+2*n_AF-1+j+1; i*(3*n_AF-1)+(2*j-1) + 1;(i-1)*(3*n_AF-1)+2*n_AF-1+j];    %global node numbers used by the element
           %see comment at bottom for node-numbering strategy 
           
           %Assemble into global stiffness matrix - Transfer from local K to global K based on node numbers used by the element
           %Reference: https://caendkoelsch.wordpress.com/2017/12/03/how-are-stiffness-matrices-assembled-in-fem/
           for p = 1:8
               for q = 1:8
                   K(2*nodeNumbers_ele(p)-1, 2*nodeNumbers_ele(q)-1) = K(2*nodeNumbers_ele(p)-1, 2*nodeNumbers_ele(q)-1) + K_ele(2*p-1, 2*q-1);
                   K(2*nodeNumbers_ele(p), 2*nodeNumbers_ele(q)) = K(2*nodeNumbers_ele(p), 2*nodeNumbers_ele(q)) + K_ele(2*p, 2*q);
                   K(2*nodeNumbers_ele(p)-1, 2*nodeNumbers_ele(q)) = K(2*nodeNumbers_ele(p)-1, 2*nodeNumbers_ele(q)) + K_ele(2*p-1, 2*q);
                   K(2*nodeNumbers_ele(p), 2*nodeNumbers_ele(q)-1) = K(2*nodeNumbers_ele(p), 2*nodeNumbers_ele(q)-1) + K_ele(2*p, 2*q-1);
               end
           end
           
           %Generate local f matrix:
           if y_ele(1) == 0 && x_ele(1)<nodeB       %Bottom edge
               f_ele(1)     = hlen_ele*tau_AB/6;
               f_ele(2)     = hlen_ele*s_AB/6;
               f_ele(3)     = 4*hlen_ele*tau_AB/6;
               f_ele(4)     = 4*hlen_ele*s_AB/6;
               f_ele(9)     = hlen_ele*tau_AB/6;
               f_ele(10)    = hlen_ele*s_AB/6;
           elseif y_ele(1) == 0 && x_ele(1)<nodeC
               f_ele(1)     = hlen_ele*tau_BC/6;
               f_ele(3)     = 4*hlen_ele*tau_BC/6;
               f_ele(9)     = hlen_ele*tau_BC/6;
               %Specify nothing for ty. Anyway this part has v = 0
           elseif y_ele(1) == 0 && x_ele(1)<Lx
               f_ele(1)     = hlen_ele*tau_CD/6;
               f_ele(2)     = hlen_ele*s_CD/6;
               f_ele(3)     = 4*hlen_ele*tau_CD/6;
               f_ele(4)     = 4*hlen_ele*s_CD/6;
               f_ele(9)     = hlen_ele*tau_CD/6;
               f_ele(10)    = hlen_ele*s_CD/6;
           end

           if x_ele(1) ==  0 %Left Edge
               f_ele(1)     = f_ele(1) + vlen_ele*s_AF/6;
               f_ele(2)     = f_ele(2) + vlen_ele*tau_AF/6;
               f_ele(7)     = f_ele(7) + vlen_ele*s_AF/6;
               f_ele(8)     = f_ele(8) + vlen_ele*tau_AF/6;
               f_ele(15)    = f_ele(15) + 4*vlen_ele*s_AF/6;
               f_ele(16)    = f_ele(16) + 4*vlen_ele*tau_AF/6;
           elseif x_ele(2) == Lx %Right Edge
               f_ele(3)     = f_ele(3) + vlen_ele*s_ED/6;
               f_ele(4)     = f_ele(4) + vlen_ele*tau_ED/6;
               f_ele(5)     = f_ele(5) + vlen_ele*s_ED/6;
               f_ele(6)     = f_ele(6) + vlen_ele*tau_ED/6;
               f_ele(11)    = f_ele(11) + 4*vlen_ele*s_ED/6;
               f_ele(12)    = f_ele(12) + 4*vlen_ele*tau_ED/6;
           end

           if y_ele(3) == Ly %Top Edge
               f_ele(5)     = f_ele(5) + hlen_ele*tau_EF/6;
               f_ele(6)     = f_ele(6) + hlen_ele*s_EF/6;
               f_ele(7)     = f_ele(7) + hlen_ele*tau_EF/6;
               f_ele(8)     = f_ele(8) + hlen_ele*s_EF/6;
               f_ele(13)    = f_ele(13) + 4*hlen_ele*tau_EF/6;
               f_ele(14)    = f_ele(14) + 4*hlen_ele*s_EF/6;
           end
       
           %Transfer from local K to global K based on node numbers used by the element
           for p = 1:8
              f(nodeNumbers_ele(p)*2-1) = f(nodeNumbers_ele(p)*2-1) + f_ele(2*p-1);
              f(nodeNumbers_ele(p)*2)   = f(nodeNumbers_ele(p)*2) + f_ele(2*p);
           end
  
       end
    end
    
    %small piece of code to check banded nature of K:
    %    A = full(K);
    %    A = ~~A;
    %    figure; hAxes = gca; imagesc(A);
    
    
    %Get reduced matrix
    roller_nodes_rows         = 2*(2*n_AB+1+2) : 2: 2*(2*n_AB+1 + 2 + 2*(n_BC+1));
    
    %Since without specifying u in atleast one place we would still have a
    %rank-deficient system, we set u at mid point of bottom edge BC to be 0 (since symmatric)
    %if the following line is commented/removed, then K_reduced will be non-invertible
    roller_nodes_rows         = [roller_nodes_rows 2*n_AB+1+2+2*n_AB+1 + 2 + 2*(n_BC+1)-1];
   
    
    K_reduced                       = K;
    K_reduced(roller_nodes_rows,:)  = [];
    K_reduced(:,roller_nodes_rows)  = [];
    f_reduced                       = f;
    f_reduced(roller_nodes_rows)    = [];
    
    
    %Solve:
    a_reduced                       = K_reduced\f_reduced;
    
    %Here we need to put the known values of displacement back in a_reduced
        alpha = a_reduced(1:4*(n_AB+ 1));
        for i = 4*(n_AB+ 1) + 1:2*n_AB+1+2+2*n_AB+1 + 2 + 2*(n_BC+1)-1 - 1
            beta    = [a_reduced(i) ; 0];
            alpha   = [alpha ; beta];
        end
        b           = [alpha ; a_reduced(2*n_AB+1+2+2*n_AB+1 + 2 + 2*(n_BC+1)-1:2*(2*n-1)-2*n_BC-4)];
    a_final         = [b(1:2*n_AB+1+2+2*n_AB+1 + 2 + 2*(n_BC+1)-1-1); 0 ; 0 ; b(2*n_AB+1+2+2*n_AB+1 + 2 + 2*(n_BC+1)-1:end); a_reduced(2*(2*n-1)-2*n_BC-4+1:end)];
    
    %Flux (stress) calculation
    s_global_x          = zeros(3*n^2 - 2*n, 1);
    s_global_y          = zeros(3*n^2 - 2*n, 1);
    tau_global_xy       = zeros(3*n^2 - 2*n, 1);
    for i = 1:(length(y) - 1)
        for j = 1:(length(x) - 1) 
            s_nodal_x           = zeros(8, 1);    %initialise
            s_gauss_x           = zeros(9, 1);
            s_nodal_y           = zeros(8, 1);    %initialise
            s_gauss_y           = zeros(9, 1);
            tau_nodal_xy        = zeros(8, 1); %initialise
            tau_gauss_xy        = zeros(9, 1);
            
            x_ele               = [x(j); x(j+1); x(j+1); x(j)]; %corner points of element counter-clockwise from bottom left
            y_ele               = [y(i); y(i); y(i+1); y(i+1)]; %corner points of element counter-clockwise from bottom left
           
            %connectivity matrix:
            nodeNumbers_ele     = [(i-1)*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 2; i*(3*n_AF-1)+(2*j-1)+2; i*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 1; (i-1)*(3*n_AF-1)+2*n_AF-1+j+1; i*(3*n_AF-1)+(2*j-1) + 1;(i-1)*(3*n_AF-1)+2*n_AF-1+j];    %global node numbers used by the element
            %see comment at bottom for node-numbering strategy 
            
            a_ele               = zeros(16, 1);
            for p = 1:8
                a_ele(2*p-1)    = a_final(2*nodeNumbers_ele(p)-1);
                a_ele(2*p)      = a_final(2*nodeNumbers_ele(p));
            end
            
            for p = 1:length(GP)
               for q = 1:length(GP)
                   B        = [];
                   xi       = GP(q); %specifically written as GP(q) not GP(p), since order of Gauss Points important here
                   eta      = GP(p);
                   J_ele    = Jacobian(x_ele, y_ele, xi, eta);
                   for k = 1:8 %8 trial functions
                    B       = [B diag(J_ele^-1*delPhiK_delIso(k, xi, eta))];
                   end
                   B3           = [B(2,2:2:end);B(1,1:2:end)];
                   B(3,:)       = (B3(:))';
                   tau_gauss    = C*B*a_ele;
                   s_gauss_x((p-1)*length(GP) + q)      = tau_gauss(1);
                   s_gauss_y((p-1)*length(GP) + q)      = tau_gauss(2);
                   tau_gauss_xy((p-1)*length(GP) + q)   = tau_gauss(3);
               end
            end
            
            %Nodal flux values for element:
            s_nodal_x           = TR*s_gauss_x;
            s_nodal_y           = TR*s_gauss_y;
            tau_nodal_xy        = TR*tau_gauss_xy;
            
            %Flux averaging part 1 - assemble to global flux vector 
            for p = 1:8
              s_global_x(nodeNumbers_ele(p))        = s_global_x(nodeNumbers_ele(p)) + s_nodal_x(p);
              s_global_y(nodeNumbers_ele(p))        = s_global_y(nodeNumbers_ele(p)) + s_nodal_y(p);
              tau_global_xy(nodeNumbers_ele(p))     = tau_global_xy(nodeNumbers_ele(p)) + tau_nodal_xy(p);
            end
            
        end
    end
    
    %create a variable in which we store the x and y value of each node
    coords          = zeros(3*n^2 - 2*n, 2);
    for p = 1:3*n^2 - 2*n
        row         = ceil(p/(3*n-1));
        within_row  = p - (row-1)*(3*n-1);
        
        %x coordinate
        if within_row > 2*n-1
           coords(p, 1) = x(within_row - (2*n-1));
        else
           coords(p, 1) = x_withmid(within_row);
        end
        
        %y coordinate
        toadd = 0;
        if within_row > 2*n-1
            toadd = Ly/(n-1)/2;
        end
        coords(p, 2) = (row-1)*Ly/(n-1)+ toadd;
    end
    
    %Flux averaging part 2 - divide by 4 for the interior nodes, by 2 for
    %the edge nodes and by leave the corner nodes alone
    for p = 1:3*n^2 - 2*n
        row         = ceil(p/(3*n-1));
        within_row  = p - (row-1)*(3*n-1);
        
        if (within_row > 2*n-1) || (within_row < 2*n-1 && mod(within_row, 2)==0)
            if ~(coords(p, 1) == 0 || coords(p, 1) == 1)
               if ~(coords(p, 2) == 0 || coords(p, 2)==1)
                   s_global_x(p)            = s_global_x(p)/2;
                   s_global_y(p)            = s_global_y(p)/2;
                   tau_global_xy(p)         = tau_global_xy(p)/2;
               end
            end
        else
            if coords(p, 1) == 0 || coords(p, 1) == 1
               if ~(coords(p, 2) == 0 || coords(p, 2)==1)
                    s_global_x(p)           = s_global_x(p)/2;
                    s_global_y(p)           = s_global_y(p)/2;
                    tau_global_xy(p)        = tau_global_xy(p)/2;
               end
            else
               if (coords(p, 2) == 0 || coords(p, 2)==1)
                    s_global_x(p)           = s_global_x(p)/2;
                    s_global_y(p)           = s_global_y(p)/2;
                    tau_global_xy(p)        = tau_global_xy(p)/2;
               else
                   s_global_x(p)            = s_global_x(p)/4;
                   s_global_y(p)            = s_global_y(p)/4;
                   tau_global_xy(p)         = tau_global_xy(p)/4;
               end
            end
        end
    end
    
    %Find the flux at the 7x7 test points
    for p = 1:length(xnodes_convcheck)
       for q = 1:length(ynodes_convcheck)
           testnode_x = xnodes_convcheck(p);
           testnode_y = ynodes_convcheck(q);
           
           %find the element in which the test point lies
           row = intersect(find(coords(:, 1) == x(find(x > testnode_x, 1, 'first'))),  find(coords(:, 2) == y(find(y > testnode_y, 1, 'first'))));
           i = ceil(row/(3*n-1));
           j = row - (i-1)*(3*n-1);
           j = (j-1)/2;
           nodeNumbers_ele = [(i-1)*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 2; i*(3*n_AF-1)+(2*j-1)+2; i*(3*n_AF-1)+(2*j-1); (i-1)*(3*n_AF-1)+(2*j-1) + 1; (i-1)*(3*n_AF-1)+2*n_AF-1+j+1; i*(3*n_AF-1)+(2*j-1) + 1;(i-1)*(3*n_AF-1)+2*n_AF-1+j];    %global node numbers used by the element
           
           %do a 2D interpolation using the flux at the nodes:
           x_ele = coords(:,1);
           y_ele = coords(:,2);
           x_ele = x_ele(nodeNumbers_ele);
           y_ele = y_ele(nodeNumbers_ele);
           
           F = scatteredInterpolant(x_ele, y_ele, s_global_x(nodeNumbers_ele));
           current_iteration_deriv_sx(p, q) = F(testnode_x, testnode_y);
           F = scatteredInterpolant(x_ele, y_ele, s_global_y(nodeNumbers_ele));
           current_iteration_deriv_sy(p, q) = F(testnode_x, testnode_y);
           F = scatteredInterpolant(x_ele, y_ele, tau_global_xy(nodeNumbers_ele));
           current_iteration_deriv_tauxy(p, q) = F(testnode_x, testnode_y);
           
       end
    end
    
    isConverged                     = norm(current_iteration_deriv_sx-prev_iteration_deriv_sx)/norm(current_iteration_deriv_sx)< tolerance && ...
        norm(current_iteration_deriv_sy-prev_iteration_deriv_sy)/norm(current_iteration_deriv_sy)< tolerance && ...
        norm(current_iteration_deriv_tauxy-prev_iteration_deriv_tauxy)/norm(current_iteration_deriv_tauxy)< tolerance;

    prev_iteration_deriv_sx            = current_iteration_deriv_sx;
    prev_iteration_deriv_sy            = current_iteration_deriv_sy;
    prev_iteration_deriv_tauxy         = current_iteration_deriv_tauxy;
    
    %update n values for next iteration (h convergence)
    n_AF                = n_AF +iter_step;        
    
    n_AB                = (n_AF - 4)/4;  
    n_BC                = (n_AF - 4)/2;  
    n_CD                = (n_AF - 4)/4; 
    n                   = n_AF;
    
    n_iters = n_iters + 1;
    
    if n_iters>max_iters %so that the program doesn't enter an infinite loop in case it doesn't converge fast enough
        break
    end
end
    n_iters             = n_iters - 1;
    n_AF                = n_AF - iter_step;        
    
    n_AB                = (n_AF - 4)/4;  
    n_BC                = (n_AF - 4)/2;  
    n_CD                = (n_AF - 4)/4; 
    n                   = n_AF;


%Plotting:
scatter(coords(:,1), coords(:,2),'k.')
xlim([-0.2 1.2]);
ylim([0 1.4]);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
displaced_x = coords(:,1) + a_final(1:2:end);
displaced_y = coords(:,2) + a_final(2:2:end);
figure();
scatter(displaced_x, displaced_y, 'r.');
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
% 
figure;
plot3(coords(:,1), coords(:,2), a_final(1:2:end), 'k.');
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('u (m)','FontSize',14);
grid on;
pbaspect([1 1 1]);
figure;
plot3(coords(:,1), coords(:,2), a_final(2:2:end), 'k.');
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('v (m)','FontSize',14);
grid on;
pbaspect([1 1 1]);
% 
[Y,X] = meshgrid(coords(:,2), coords(:,1));
figure;
vq = griddata(coords(:,1), coords(:,2), s_global_x, X, Y);
surf(X, Y,vq);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
shading interp;
colormap jet;
colorbar;
view(2);
pbaspect([1 1 1]);
% 
[Y,X] = meshgrid(coords(:,2), coords(:,1));
figure;
vq = griddata(coords(:,1), coords(:,2), s_global_y, X, Y);
surf(X, Y,vq);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
shading interp;
colormap jet;
colorbar;
view(2);
pbaspect([1 1 1]);
% 
colormap jet;
cmap = jet(256);
v = rescale(tau_global_xy, 1, 256); 
numValues = length(tau_global_xy);
markerColors = zeros(numValues, 3);
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end
% Create the scatter plot.
figure;
scatter(coords(:,1), coords(:,2), 40, markerColors, 'filled' );
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
%colorbar;
view(2);
pbaspect([1 1 1]);

% Relative error plotting:
% s_x = [1 0.077477074 0.039802251 0.024057073 0.021725196 0.019342744 0.017647205];
% s_y = [1 0.011023764 0.006873848 0.005474206 0.004885956 0.003423156 0.002966996];
% tau_xy = [1 0.064989227 0.052245499 0.029071743 0.018169344 0.015023063 0.008479651];
% N = [28 32 36 40 44 48 52];
% semilogy(N, s_x, 'b.-', 'LineWidth', 2);
% hold on;
% semilogy(N, s_y, 'r.-', 'LineWidth', 2);
% hold on;
% semilogy(N, tau_xy, 'k.-', 'LineWidth', 2);
% xlabel("N", 'FontSize', 12)
% ylabel("Relative error", 'FontSize', 12)
% yline(0.018,'-','\epsilon = 0.018', 'LineWidth', 1.5)
% legend('e_1','e_2','e_3','', 'FontSize', 12)
% grid on;
% pbaspect([1 1 1])



function toReturn = delPhiK_delIso(k, xi, eta)
    switch k
        case 1
            delPhidelXi     = 0.25*(1-eta)*(2*xi+eta);
            delPhidelEta    = 0.25*(1-xi)*(2*eta+xi);
        case 2
            delPhidelXi     = 0.25*(1-eta)*(2*xi-eta);
            delPhidelEta    = 0.25*(1+xi)*(2*eta-xi);
        case 3
            delPhidelXi     = 0.25*(1+eta)*(2*xi+eta);
            delPhidelEta    = 0.25*(1+xi)*(2*eta+xi);
        case 4
            delPhidelXi     = 0.25*(1+eta)*(2*xi-eta);
            delPhidelEta    = 0.25*(1-xi)*(2*eta-xi);
        case 5
            delPhidelXi     = -xi*(1-eta);
            delPhidelEta    = -0.5*(1-xi^2);
        case 6
            delPhidelXi     = 0.5*(1-eta^2);
            delPhidelEta    = -(1+xi)*eta;
        case 7
            delPhidelXi     = -xi*(1+eta);
            delPhidelEta    = 0.5*(1-xi*xi);
        case 8
            delPhidelXi     = -0.5*(1-eta^2);
            delPhidelEta    = -(1-xi)*eta;
    end
    toReturn = [delPhidelXi delPhidelEta]';
end
function toReturn = Jacobian(x, y, xi, eta) %assume 4 points at corners, 4 at centres of edges
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    y1 = y(1);
    y2 = y(2);
    y3 = y(3);
    y4 = y(4);
    delXdelXi   = 0.25*((x2 - x1)*(1 - eta) + (x3 - x4)*(1 + eta));
    delXdelEta  = 0.25*((x4 - x1)*(1 - xi) + (x3 - x2)*(1 + xi));
    delYdelXi   = 0.25*((y2 - y1)*(1 - eta) + (y3 - y4)*(1 + eta));
    delYdelEta  = 0.25*((y4 - y1)*(1 - xi) + (y3 - y2)*(1 + xi));
    
    toReturn    = [delXdelXi delYdelXi;delXdelEta delYdelEta];
end

% Node numbering strategy followed:
% Left to right bottom last row, left to right second last row ...
% 
% 34 35 36 37 38 39 40
% 30    31    32    33
% 23 24 25 26 27 28 29
% 19    20    21    22   row with midpoints of edges 
% 12 13 14 15 16 17 18
%  8     9    10    11   row with midpoints of edges
%  1  2  3  4  5  6  7

