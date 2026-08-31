clear all
close all
clc
a  = 1e-2;
er = 1;
e0 = 8.854e-12;
m0 = pi*4e-7;
c0 = 3e8;
f  = 35.5e9;

pnm  = [2.405 3.832 5.135 5.520 6.380 7.016];

gd = [1;
      0; 
      0;
      a];
d1 = decsg(gd);

[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);

%pdeplot(p,e,t); axis equal ; axis tight;

Nm      = 6;         %Nmodes
Nn      = size(p,2); %Nnodes 
Ne      = size(t,2); %Nelements
Nd      = size(e,2); %Nedges
node_id = ones(Nn,1);
X0      = zeros(Nn,1);

%loop se oles tis akmes gia na orisoume synthikes dirichlet ston kymatodigo
for id  = 1:Nd
    n1  = e(1,id);
    n2  = e(2,id);
    r1  = e(6,id);
    r2  = e(7,id);
    x1  = p(1,n1);
    y1  = p(2,n1);
    x2  = p(1,n2);
    y2  = p(2,n2);

    if(r1==0 ||r2==0) 
        node_id(n1) = 0; 
        node_id(n2) = 0;
    end;
end;

counter = 0;
index   = zeros(Nn,1);
for in  = 1:Nn
    if(node_id(in) == 1)
        counter     = counter + 1;
        index(in)   = counter;
    end;
end;
Nf = counter  %Number of unknowns

S  = spalloc(Nf,Nf,7*Nf);
T  = spalloc(Nf,Nf,7*Nf);

% Loop pou sarwnei ola ta stoixeia, gia tin teliki synathroisi :
for ie = 1:Ne
    n(1:3) = t(1:3,ie); % Global aroithmisi kathe komvou
    x(1:3) = p(1,n(1:3));
    y(1:3) = p(2,n(1:3));
    region = t(4,ie);
    De     = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    b(1)   = (y(2)-y(3))/De;
    b(2)   = (y(3)-y(1))/De;
    b(3)   = (y(1)-y(2))/De;
    c(1)   = (x(3)-x(2))/De;
    c(2)   = (x(1)-x(3))/De;
    c(3)   = (x(2)-x(1))/De;
    Ae     = abs(De/2); %To embadon tou trigwnou

    %Ypologismos topikou pinaka
    for i  = 1:3 
        for j = 1:3 
            Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;
            if(i==j)
                Te(i,j) = Ae/6;
            else
                Te(i,j) = Ae/12;
            end;
            if(node_id(n(i)) == 1 )
                if(node_id(n(j)) == 1 )
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j); %Προσθέτουμε στον ολικό πίνακα , την συνεισφορά του τοπικού.
                    T(index(n(i)),index(n(j))) = T(index(n(i)),index(n(j))) + Te(i,j); %Ολικός Πίνακας Μάζας
                end;
            end;
        end;
    end;
end;

[Ez,kc2] = eigs(S,T,12,'sm','MaxIterations',350);
kc2 = kc2(:, [1 3 4 6 8 10]);
Ez  = Ez(:, [1 3 4 6 8 10]);

for i = 1:12
        nonzero_values = kc2(kc2 ~= 0); 
        D = diag(nonzero_values);
    end;

for im = 1: Nm
    for in = 1:Nn
        if (node_id(in)==1)
            X0(in) = Ez(index(in),im);
        end;
    end;  
   fc_FEM(im) = c0*sqrt(D(im,im))/(2*pi)
   fc(im)     = (c0*pnm(im))/(2*pi*a)
   error_per_cent(im) = 100*(fc_FEM-fc)/fc
   figure;
   pdeplot(p,e,t,'xydata',X0,'colorbar','on'); colormap jet; axis equal; axis tight;
   title(sprintf('TM Mode %d', im));

   [EZx,EZy] = pdegrad(p,t,X0);
   Ex      = (-1j*(sqrt((2*pi*f)^2*m0*e0 - D(im,im))*EZx))/D(im,im);
   Ey      = (-1j*(sqrt((2*pi*f)^2*m0*e0 - D(im,im))*EZy))/D(im,im);
   E       = sqrt(abs(Ex).^2 + abs(Ey).^2); %Εγκάρσιο Ε στην διατομή
   figure()
   pdeplot(p,e,t,'xydata',E,'colorbar','on'); colormap jet; axis equal; axis tight;
   title(sprintf('TM Mode %d (Electric Field)', im))
end;
    

   