clear all
clc
w = 4e-2;
d = 1e-2;
t = 2e-3;
A = 7*w;
B = 7*w;
V = 100;
er=[1,2.2];
e0 = 8.854e-12;

gd=[3    3      3      3;
    4    4      4      4;
   -A/2  w/2    w/2    w/2;
    A/2  w/2    w/2    w/2; 
    A/2 -w/2   -w/2   -w/2;
   -A/2 -w/2   -w/2   -w/2; 
   -B/2  d/2   -d/2-t -d/2;
   -B/2  d/2+t -d/2    d/2;
    B/2  d/2+t -d/2    d/2;
    B/2  d/2   -d/2-t -d/2];

ns=[82 82 82 82;
    49 50 51 52];
sf='R1-R2-R3+R4';
d1=decsg(gd,sf,ns);

[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
% [p,e,t] = refinemesh(d1,p,e,t);
% [p,e,t] = refinemesh(d1,p,e,t);


Nn = size(p,2); %Nnodes 
Ne = size(t,2);  %Nelements
Nd = size(e,2); %Nedges

node_id = ones(Nn,1);
X0      = zeros(Nn,1);
%loop se oles tis akmes gia na orisoume synthikes dirichlet stis plakes tou
%piknoti
for id  = 1:Nd
    n1  = e(1,id);
    n2  = e(2,id);
    r1  = e(6,id);
    r2  = e(7,id);
    x1  = p(1,n1);
    y1  = p(2,n1);
    x2  = p(1,n2);
    y2  = p(2,n2);

    if((r1==0 ||r2==0) && (y1>0 && y1<B/4 && abs(x1)<1.2*w/2 && y2>0 && y2<B/4 && abs(x2)<1.2*w/2))
        node_id(n1) = 0; 
        node_id(n2) = 0;
        X0(n1)      = V/2;
        X0(n2)      = V/2;
    end;
    if((r1==0 ||r2==0) && (y1<0 && y1>-B/4 && abs(x1)<1.2*w/2 && y2<0 && y2>-B/4 && abs(x2)<1.2*w/2))
        node_id(n1) = 0; 
        node_id(n2) = 0;
        X0(n1)      =-V/2;
        X0(n2)      =-V/2;
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

pdeplot(p,e,t); axis equal ; axis tight;
% for in = 1:Nn
%     x  = p(1,in); 
%     y  = p(2,in);
%     text(x,y,num2str(X0(in)));
% end;

% for ie = 1:Ne
%     n(1:3) = t(1:3,ie); %oi komvoi tou kathe stoixiou
%     x(1:3) = p(1,n(1:3));
%     y(1:3) = p(2,n(1:3));
%     region = t(4,ie);
%     %Συντεταγμένες Βαρύτκεντρου 
%     xcm  = sum(x)/3; 
%     ycm  = sum(y)/3;
%     text(xcm,ycm,num2str(region));
% end;

S = spalloc(Nf,Nf,7*Nf);
B = zeros(Nf,1);
% Loop pou sarwnei ola ta stoixeia, gia tin teliki synathroisi :
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

    %Ypologismos topikou pinaka
    for i  = 1:3 
        for j = 1:3 
            Se(i,j) = e0*er1*(b(i)*b(j)+c(i)*c(j))*Ae;
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

X=S\B; %Βρίσκω τα Χ (το δυναμικό) σε όλους τους κομβους στους οποίους ήταν άγνωστο

for in = 1:Nn
    if(node_id(in) == 1)
        X0(in)=X(index(in));
    end;
end;

figure(1)
pdeplot(p,e,t,'xydata',X0,'contour','on'); axis equal; axis tight; colormap jet;
title('Μεταβολή του δυναμικού Φ:', 'FontSize',16,'FontWeight','bold')

figure(2)
pdegplot(d1);
hold on;

[Ex,Ey] = pdegrad(p,t,X0);
pdeplot(p,e,t,'FlowData',-[Ex;Ey]); axis equal; axis tight;
title('Διανυσματικό Πεδίο Ε:', 'FontSize',16,'FontWeight','bold')
%Fu=pdeInterpolant(p,t,X0)

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
    We_grad= We_grad + 0.5*er1*e0*E2*Ae;


    for i  = 1:3 
        for j = 1:3 
            Se(i,j) = e0*er1*(b(i)*b(j)+c(i)*c(j))*Ae;
            Spq(n(i),n(j)) = X0(n(i))*Se(i,j)*X0(n(j)) + Spq(n(i),n(j));
        end;
    end;
end;

We     = (1/2)*Spq;

We_tot = sum(We(:));

C_exact   = (e0*er(2)*w)/d  

C_FEM   = (2*We_tot)/V^2 

C_grad = (2*We_grad)/V^2;

relative_error_percent = 100*abs(C_FEM - C_exact)/C_exact 

                




