clear all
clc
a  = 0.76e-3;
bi  = 1.75e-3;
er = 1;
e0 = 8.854e-12;
V  = 1;
tol= 10e-5;
max= 1e4;

gd = [1 1;
      0 0;
      0 0;
      a bi];

ns = [82 82 ;
      49 50];
sf = 'R2-R1';
d1=decsg(gd,sf,ns);

[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);

Nn      = size(p,2); %Nnodes 
Ne      = size(t,2);  %Nelements
Nd      = size(e,2); %Nedges
node_id = ones(Nn,1);
X0      = zeros(Nn,1);

pdeplot(p,e,t); axis equal ; axis tight;

%loop se oles tis akmes gia na orisoume synthikes dirichlet stous agwgous
%tou omoaxonikou kalodiou
for id  = 1:Nd
    n1  = e(1,id);
    n2  = e(2,id);
    r1  = e(6,id);
    r2  = e(7,id);
    x1  = p(1,n1);
    y1  = p(2,n1);
    x2  = p(1,n2);
    y2  = p(2,n2);

    if((r1==0 ||r2==0) && (sqrt(x1^2+y1^2)<1.1*a))
        node_id(n1) = 0; 
        node_id(n2) = 0;
        X0(n1)      = V;
        X0(n2)      = V;
    end;
    if((r1==0 ||r2==0) && (sqrt(x1^2+y1^2)>0.9*bi))
        node_id(n1) = 0; 
        node_id(n2) = 0;
        X0(n1)      = 0;
        X0(n2)      = 0;
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
Nf   = counter  %Number of unknowns

S = spalloc(Nf,Nf,7*Nf);
B = zeros(Nf,1);
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
    Ae     = abs(De/2); %to embadon tou trigwnou

    %Ypologismos topikou pinaka
    for i  = 1:3 
        for j = 1:3 
            Se(i,j) = e0*er*(b(i)*b(j)+c(i)*c(j))*Ae;
            if(node_id(n(i)) == 1 )
                if(node_id(n(j)) == 1 )
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j); %Προσθέτουμε στον ολικό πίνακα , την συνεισφορά του τοπικού.
                else 
                    B(index(n(i))) = B(index(n(i))) -Se(i,j)*X0(n(j)); %j γνωστοί
                end;
            end;
        end;
    end;
end;

 tic
 X=S\B;
% X=bicg(S,B,tol,max); %Βρίσκω τα Χ (το δυναμικό) σε όλους τους κομβους στους οποίους ήταν άγνωστοι
 toc

for in = 1:Nn
    if(node_id(in) == 1)
        X0(in)=X(index(in));
    end;
end;

figure(1)
pdeplot(p,e,t,'xydata',X0,'contour','on'); axis equal; axis tight; colormap jet;
title('Μεταβολή του δυναμικού Φ:', 'FontSize',16,'FontWeight','bold');

figure(2)
pdegplot(d1);
hold on;
[Ex,Ey] = pdegrad(p,t,X0);
pdeplot(p,e,t,'FlowData',-[Ex;Ey]); axis equal; axis tight;
title('Διανυσματικό Πεδίο Ε:', 'FontSize',16,'FontWeight','bold');

Spq = spalloc(Nn,Nn,7*Nn);
We_grad=0;
for ie = 1:Ne
    n(1:3) = t(1:3,ie); % Global aroithmisi kathe komvou
    x(1:3) = p(1,n(1:3));
    y(1:3) = p(2,n(1:3));
    region = t(4,ie);
    er1    = er(region);
    De     = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    b(1)   = (y(2)-y(3))/De;
    b(2)   = (y(3)-y(1))/De;
    b(3)   = (y(1)-y(2))/De;
    c(1)   = (x(3)-x(2))/De;
    c(2)   = (x(1)-x(3))/De;
    c(3)   = (x(2)-x(1))/De;
    Ae     = abs(De/2); %to embadon tou trigwnou

    E2     = Ex(ie)^2 + Ey(ie)^2;


    for i  = 1:3 
        for j = 1:3 
            Se(i,j) = e0*er1*(b(i)*b(j)+c(i)*c(j))*Ae;
            Spq(n(i),n(j)) = X0(n(i))*Se(i,j)*X0(n(j)) + Spq(n(i),n(j));
        end;
    end;
end;

We     = (1/2)*Spq;

We_tot = sum(We(:))

C_Ni     = (2*We_tot)/V^2 

C_exact  = (2*pi*e0*er)/(log(bi/a))

error    = C_Ni - C_exact


