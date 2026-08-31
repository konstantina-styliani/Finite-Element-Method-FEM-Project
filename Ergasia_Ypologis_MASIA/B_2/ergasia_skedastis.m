clear all
clc
l     = 1;
a_vals = [l/2, 2*l, 5*l];
%R     = l +a;
%R     = 2*l +a;
e0    = 8.854e-12;
m0    = pi*4e-7;
c0    = 3e8;
f     = 300e6;
k     = (2*pi)/l;
omega = 2*pi*f;
E0    = 1

fig_count = 0; 
for ia = 1:3
    a  = a_vals(ia)
Rvals = [(l/2 + a),(l +a),(2*l +a)]; %3 R gia kathe a

for m=1:3
    R=Rvals(m);


aR = -1i*k - 1/(2*R); %σταθερός συντελεστής 1ης τάξης ABC


gd = [1 1;
      0 0;
      0 0;
      a R];

ns = [82 82 ;
      49 50];
sf = 'R2-R1';
d1=decsg(gd,sf,ns);
[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);

Nn      = size(p,2); %Nnodes 
Ne      = size(t,2);  %Nelements
Nd      = size(e,2); %Nedges
node_id = ones(Nn,1);
node_C2 = ones(Nn,1);

EzI   = E0 * exp(-1i*k*p(1,:));
EzS      = zeros(Nn,1);
%loop se oles tis akmes gia na orisoume synthikes dirichlet ston teleia
%agwgimo kylindro , kai na vroume poioi komvoi einai sto exwteriko orio C2
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
        EzS(n1)     =-EzI(n1);
        EzS(n2)     =-EzI(n2);
    end;
    if((r1==0 ||r2==0) && (sqrt(x1^2+y1^2)>0.9*R))
        node_C2(n1) = 0; 
        node_C2(n2) = 0;
    end;
end;

% for in = 1:Nn
%     x  = p(1,in); 
%     y  = p(2,in);
%     text(x,y,num2str(node_C2(in)));
% end;

counter = 0;
index   = zeros(Nn,1);
for in  = 1:Nn
    if(node_id(in) == 1)
        counter     = counter + 1;
        index(in)   = counter;
    end;
end;
Nf = counter  %Number of unknowns

A   = spalloc(Nf,Nf,7*Nf);
B   = zeros(Nf,1);


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
            Se(i,j) = (1/m0)*(b(i)*b(j)+c(i)*c(j))*Ae;
            if(i==j)
                Te(i,j) = e0*(Ae/6);
            else
                Te(i,j) = e0*(Ae/12);
            end;
            AE(i,j) = Se(i,j) - (omega^2)*Te(i,j);
            if(node_id(n(i)) == 1 )
                if(node_id(n(j)) == 1 )
                    A(index(n(i)),index(n(j))) = A(index(n(i)),index(n(j))) + AE(i,j);
                else
                    B(index(n(i)))             = B(index(n(i))) + AE(i,j)*EzI(n(j));

                end;
            end;
        end;
    end;
end;

coun = 0;
for in  = 1:Nn
    if(node_C2(in) == 0)
        coun     = coun + 1;
    end;
end;
coun=coun %komvoi sto synoro C2

%ABC 1ης τάξης:
for id = 1:Nd
    n1  = e(1,id);
    n2  = e(2,id);
    x1  = p(1,n1);
    y1  = p(2,n1);
    x2  = p(1,n2);
    y2  = p(2,n2);
    for i  = 1:2 
        for j = 1:2 
            if(node_C2(n1) == 0 && node_C2(n2) == 0 && node_id(n1)==1 && node_id(n2)==1) %An einai sto orio C2 kai oi dyo komvoi
                    Le      = sqrt((x1 - x2)^2 + (y1 - y2)^2);
                    if(i==j)
                        Tc2(i,j) = (1/m0)*(Le/3);
                    else
                        Tc2(i,j) = (1/m0)*(Le/6);
                    end;
                A(index(n1),index(n2)) = A(index(n1),index(n2)) - aR*Tc2(i,j);
                end;
                 
            end;
            
        end;
    end;

% K = (A + (-1i*k - (1/(2*R)))*Tc2 );

Ez =A\B

for in = 1:Nn
    if(node_id(in) == 1)
        EzS(in)=Ez(index(in));
    end;
end;

Eol =zeros(Nn,1)
for in = 1:Nn
    Eol(in)=abs(EzI(in)+EzS(in))
end;

fig_count = fig_count + 1;

figure(fig_count)
pdeplot(p,e,t,'xydata',Eol); axis equal; axis tight; colormap jet;
title(sprintf('Συνολικό Ηλεκτρικό Πεδίο E, για R = %.1f και a = %.1f', R, a));
end;
end;

