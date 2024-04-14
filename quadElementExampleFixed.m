%Robert Liu, 260981372
%I, Robert Liu, would like to exercise my right to use the late penalty waiver for this assignment. 
%Please note that all the answers to the questions are defined by comments with the question number, such as `%QUESTION 1` or `%QUESTION 5`. 

%Question 5, one line that took a line time was the stiffness calculation line which took 2.4 seconds (line 135)
%After using num2cell and cellfun for the function defined in step 1 and Gfun, 
%the time it takes now is 1.8 seconds (line 139)
%the time it takes to compute the forces as well improved from 0.49 seconds to 0.19 seconds (directly above the stiffness matrix calculations)
%ml divide has an improved performance of 0.005 seconds
%Thus showing an increase in performance after using num2cell, cellfun and sparsing the matrices
%All changes can be found in comments with QUESTION 5 as the header


clear
% suppose we have a quad (or quads) in plane with bi-linear interpolation
gridN = 10;
dx = 2/(gridN-1); % size of elements (needed for gradient)
[X1,X2] = meshgrid(linspace(-1,1,gridN));
P = [reshape(X1,1,[]); reshape(X2,1,[]) ];
P0 = P; % initial state
Pdot = zeros(size(P));

% See P with linear indexing as a vector of alternating x y coords
ix = 1:numel(X1);
ixg = reshape( ix, size(X1) );
% define elements on this grid, CW, starting upper left in matrix, which is
% actually lower left in the plot of the grid
i = reshape(ixg(1:gridN-1,1:gridN-1),[],1);
el = [ i, i+1, i+gridN+1, i+gridN ];
N = size(el,1); % number of elements

ph = drawElements( el, P, 0 ); % keep the plot handle for faster plotting

% compute deformation gradient at a location within an element
% can start with the deformation function, given location (a,b) in [0,1]^2
% within the element.  Could also see this as big X in R^2.
Xs = sym('Xs',[2,1]); % call it Xs with s meaning symbolic
phi1 = (1-Xs(1))*(1-Xs(2)); 
phi2 = (1-Xs(1))*Xs(2); 
phi3 = Xs(1)*Xs(2); 
phi4 = Xs(1)*(1-Xs(2));

% deformaiton of point (X1,X2) in element j maps to sum phi(i) P(el(j, i))
% so:  x = phi = phi1(x)*p1 + phi2(x)*p2 + phi3(x)*p3 + phi4(x)*p4;
% and: F = [dx1/dX dx1/dX2; dx2/dX dx2/dX2]
% note we'll fix for element scaling (defined top as dx) later

% deformation gradient needs grad of 4 shape functions wrt a and b
% these have grad wrt xs(1) in first coord and grad wrt xs(2) in second
dphidxs = [ gradient( phi1, Xs ) gradient( phi2, Xs ) gradient( phi3, Xs ) gradient( phi4, Xs ) ];
dphidxsfun = @(a,b) eval(subs( dphidxs, Xs, [a;b] ));

% evaluation of shape function gradients at different quadrature points
Bp5p5 = computeB(dphidxsfun(.5,.5)); % Quadrature point right in the middle! won't be enough
Bq1to4 = [ computeB(dphidxsfun(1/3,1/3));
           computeB(dphidxsfun(1/3,2/3));
           computeB(dphidxsfun(2/3,1/3));
           computeB(dphidxsfun(2/3,2/3)) ];
Bel = Bq1to4; % choose Bp5p5 for one point, or Bq1to4 for 4 point quadrature

% fill out B matrix for a given element to compute deformaiton gradient when multiplied by P
B = zeros(N*size(Bel,1),numel(P)); % deformation gradient F has 4 entries
for j = 1:N
    eldof = [el(j,:)*2-1;el(j,:)*2];
    B((j-1)*size(Bel,1)+(1:size(Bel,1)),eldof) = Bel;
end
B = B / dx; % scale by element size as X space assumed 1x1 domain up to now


Fs = sym('Fs', [2,2],'real');
Es = 1/2*(Fs.'*Fs-eye(2));
J = det(Fs);
nu = 0.4; % Poisson's raio
E = 100; % Young's modulus
C = Fs'*Fs; % right Cauchy-Green tensor
[mu, lambda] = toLame( nu, E ); % mu (Shear modulus) and lambda (Lamé's first parameter)

neoHookeanEnergy2D = true;

%QUESTION 4
if neoHookeanEnergy2D
    % Neo-Hookean energy density
    psi = mu/2*(trace(C) - 2) - mu*log(J) + lambda/2*(log(J))^2;
else
    % St. Venant-Kirchhoff energy density
    psi = 1/2 * lambda * (trace(Es))^2 + mu * trace( Es'*Es );
end

    
G = gradient(psi,Fs(:));
Gfun = matlabFunction(G);
% CAREFUL with order of the parameters to Gfun !!!!
%COMMENTED OUT DUE TO QUESTION 5
% myGfun = @(F) Gfun(F(1),F(3),F(2),F(4));

% QUESTION 1
q1 = hessian(psi, Fs(:));
q1Func = matlabFunction(q1);
%COMMENTED OUT DUE TO QUESTION 5
% myq1Func = @(F) q1Func(F(1),F(3),F(2),F(4));
%C = hessian(psi,F(:)); % we only need this for backward Euler.  Easy to add!

% do some time stepping

% compute lumped mass
density = 25;
mdiag = zeros(numel(P),1);
for i=1:N
    ix = el(i,:);
    mdiag( ix*2-1 ) = mdiag( ix*2-1 ) + 1/4*dx^2*density;
    mdiag( ix*2 ) = mdiag( ix*2 ) + 1/4*dx^2*density;
end
Minv = diag(1./mdiag);

dt = 0.01;


useSymplecticEuler = false;

for t = 1:500
    F = B*P(:); % compute deformation gradients for current state
    forces = zeros( size(P) );
    %QUESTION 2
    stiffness = zeros(numel(P),numel(P));

    %QUESTION 5 FIRST CHANGE
    Fcell = num2cell(reshape(F, 4, []), 1);
    G = cellfun(@(F) Gfun(F(1),F(3),F(2),F(4)), Fcell, 'UniformOutput', false);
    Q = cellfun(@(F) q1Func(F(1),F(3),F(2),F(4)), Fcell, 'UniformOutput', false);

    for j=1:size(B,1)/4 % go through all the quadrature points of all elements 
        ix = j*4-3:j*4; % indicies for accessing the jth deformation gradient

        %QUESTION 2
        % forces(:) = forces(:) - B(ix,:)'*myGfun(F(ix));
        % stiffness = stiffness - B(ix,:)' * myq1Func(F(ix)) * B(ix,:);

        %QUESTION 5 SECOND CHANGE
        forces(:) = forces(:) - B(ix,:)'*G{j};
        stiffness = stiffness - B(ix,:)' * Q{j} * B(ix,:);
    end

    % add gravity
    forces(2,:) = forces(2,:) - 9.8 * ones(1,size(P,2));

    if useSymplecticEuler
        % FE integration
        Pdot(:) = Pdot(:) + dt * Minv * forces(:);
        P(:) = P(:) + dt * Pdot(:); 
    else
        %QUESTION 3
        freeIndices = setdiff(1:size(P, 2), gridN:gridN:gridN*gridN);
        freeIndicesDof = sort([2*freeIndices-1, 2*freeIndices]);

        % mdiagMatrix = diag(mdiag);

        %QUESTION 5 THIRD CHANGE
        mdiagMatrix = spdiags(mdiag, 0, numel(P), numel(P)); % Initialize as sparse matrix
        a = mdiagMatrix - dt^2 * stiffness;
        b = dt * forces(:) + dt^2 * stiffness * Pdot(:);
        a = a(freeIndicesDof, freeIndicesDof);
        a = sparse(a);
        b = b(freeIndicesDof);

        newPdot = mldivide(a,b);
        newPdotMatrix = reshape(newPdot, [2, numel(newPdot)/2]);

        Pdot(:, freeIndices) = Pdot(:, freeIndices) + newPdotMatrix;
        Pdot(:,gridN:gridN:gridN*gridN) = 0;
        P(:) = P(:) + dt * Pdot(:); 
    end

    
    % reset pinned particles to P0
    P(:,gridN:gridN:gridN*gridN) = P0(:,gridN:gridN:gridN*gridN);
    disp('Positions after resetting pinned particles:');
    disp(P);
    %P(:,[gridN,gridN*gridN]) = P0(:,[gridN,gridN*gridN]);
    drawElements( el, P, ph ); %QUESTION 5 TAKES A LOT OF TIME
end

function B = computeB( g )
    % assemble B submatrix for one element's set of dofs at given gradient
    B = zeros(4,8);
    B(1,1:2:end) = g(1,:);
    B(2,2:2:end) = g(1,:);
    B(3,1:2:end) = g(2,:);
    B(4,2:2:end) = g(2,:);
end

function ph = drawElements( el, P, ph )
    % draw the elements, or just up the data in the plot handle provided (faster)
    Pm = [P [nan; nan]];
    % create index list for plotting, where we loop through element
    % indicies, repeating the first, and then include a reference to the
    % nan point added above to break the line
    elm = [el, el(:,1), size(Pm,2)*ones(size(el(:,1))) ];
    elmr = reshape( elm', 1, [] );
    x = Pm(1,elmr);
    y = Pm(2,elmr);
    if ( ph == 0 )
        ph = plot(x,y,'-b');
        axis equal;
        axis( [-2,2,-2,2] );
    else
        ph.XData = x;
        ph.YData = y;
    end
    getframe;
end

function [mu, lambda] = toLame(nu, E)
    % TOLAME Converts Poisson ratio and nu and Young's modulus E to Lamé parameters.
    mu = E / 2 / (1 + nu);
    lambda = E * nu / (1 + nu) / (1 - 2 * nu);
end


% some handy testing code... change the initial conditions!
%  P = P - [1.5;1.5];
%  theta = -90;
%  c = cosd(theta); s = sind(theta);
%  R = [c,-s;s,c];
%  S = eye(2);
%  % S = eye(2)*1.1;
%  S = [0.5 0;0 1];
%  P = S*R*P;
% more handy testing code... change the initial conditions!
%  Pdot(:,1) = [-2;-2];
%  Pdot(:,4) = [2;2];
%  Pdot(1,:) = [-1 -1 1 1];
%  Pdot(2,:) = [-2 2 -2 2];
%  Pdot = R * Pdot;
%  Pdot = [ P(2,:); -P(1,:) ];
    % draw the deformation gradient for first element.
    % These are indeed correct as far as I can tell...
%     F = B*P(:); % compute deformation gradients for current state
%     p1 = P(:,el(1,1));
%     p2 = P(:,el(1,2));
%     p3 = P(:,el(1,3));
%     p4 = P(:,el(1,4));
%     pq1 = eval(phiFun( [1/3;1/3],p1,p2,p3,p4));
%     pq2 = eval(phiFun( [1/3;2/3],p1,p2,p3,p4));
%     pq3 = eval(phiFun( [2/3;1/3],p1,p2,p3,p4));
%     pq4 = eval(phiFun( [2/3;2/3],p1,p2,p3,p4));
%     disp([pq1 pq2 pq3 pq4]);
%     F1 = reshape(F(1:4),2,2);
%     F2 = reshape(F(4+(1:4)),2,2);
%     F3 = reshape(F(8+(1:4)),2,2);
%     F4 = reshape(F(12+(1:4)),2,2);
%     
%     hold on;
%     pts = [ pq1+F1(:,1), pq1, pq1 + F1(:,2)];
%     plot( pts(1,:), pts(2,:) );
%     pts = [ pq2+F2(:,1), pq2, pq2 + F2(:,2)];
%     plot( pts(1,:), pts(2,:) );
%     pts = [ pq3+F3(:,1), pq3, pq3 + F3(:,2)];
%     plot( pts(1,:), pts(2,:) );
%     pts = [ pq4+F4(:,1), pq4, pq4 + F4(:,2)];
%     plot( pts(1,:), pts(2,:) );
%     getframe

