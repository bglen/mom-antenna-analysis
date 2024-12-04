%       mom1.m
%
%   method of momemnts calculation of
%   WORKS WITH all examples and problems
%   need to clean up and add explanations
%   two unequal plates - it is slow because of the iterative solution
% maybe implement the direct solution which solves for V based on equating
% charge on plates
%   NEGATIVE POTENTIALS, NEGATIVE CAPACITANCE
%     the following data are used:2
%     r,  vector of charge densities on subdomains
%     n1x= number of subdomains in the x direction
%     n1y= number of subdomains in the y direction
%     ax, the x-dimension of the plate
%     ay, the y-dimension of the plate
%     vp, potential on the plate
%
clear
        disp('in the data entry that follows use a space to separate input values')
        disp('when multiple values are required they must be entered as an array, including the brackets')
ne=[0 0];
vp=[0 0];
a=[0 0];
b=[0 0];
na=[0 0];
nb=[0 0];
    npl = input('enter number of plates (1 or 2)  ---> ');
        if npl==1;
            disp('enter coordinates of the four corners of the plate in meters');
            [c1]=input('enter [x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4] ---> ');
            x1(1)=c1(1);
            y1(1)=c1(2);
            z1(1)=c1(3);
            x2(1)=c1(4);
            y2(1)=c1(5);
            z2(1)=c1(6);
            x3(1)=c1(7);
            y3(1)=c1(8);
            z3(1)=c1(9);
            x4(1)=c1(10);
            y4(1)=c1(11);
            z4(1)=c1(12);
            disp('enter no. of subdomains on sides 1 and 2 of the plate, separated by a space')
            disp('note: side 1 is the side between corners 1 and 2')
            disp('note: side 2 is the side between corners 2 and 3')
            [n1]=input(' [N1, N2] ---> ');
            na(1)=n1(1);
            nb(1)=n1(2);
        else
             disp('enter coordinates of the four corners of plate no 1')
            [c1]=input('enter [x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4] ---> ')
            disp('enter no. of subdomains on sides 1 and 2 of the plate no. 1, separated by a space')
            disp('note: side 1 is the side between corners 1 and 2')
            disp('note: side 2 is the side between corners 2 and 3')
            [n1]=input(' [N1, N2] ---> ')
                        na(1)=n1(1);
                        nb(1)=n1(2);
            disp('enter coordinates of the four corners of plate no 2')
            [c2]=input('enter [x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4] ---> ')
            disp('enter no. of subdomains on sides 1 and 2 of the plate no. 2, separated by a space')
            disp('note: side 1 is the side between corners 1 and 2')
            disp('note: side 2 is the side between corners 2 and 3')
            [n1]=input(' [N1, N2] ---> ')
                        na(2)=n1(1);
                        nb(2)=n1(2);
            x1(1)=c1(1);
            y1(1)=c1(2);
            z1(1)=c1(3);
            x2(1)=c1(4);
            y2(1)=c1(5);
            z2(1)=c1(6);
            x3(1)=c1(7);
            y3(1)=c1(8);
            z3(1)=c1(9);
            x4(1)=c1(10);
            y4(1)=c1(11);
            z4(1)=c1(12);
            x1(2)=c2(1);
            y1(2)=c2(2);
            z1(2)=c2(3);
            x2(2)=c2(4);
            y2(2)=c2(5);
            z2(2)=c2(6);
            x3(2)=c2(7);
            y3(2)=c2(8);
            z3(2)=c2(9);
            x4(2)=c2(10);
            y4(2)=c2(11);
            z4(2)=c2(12);
        end
        vol=input('enter potential difference between plates or potntial on single plate  ---> ');
%        disp('notes: it is assumed that if a single plate is specified the potential difference is between the plate and infinity')
%        disp('if the two plates have different areas, only the potential on plate no. 1 is used'\n)
%        disp('the potential on the second plate is calculated assuming v1 to be the potential difference between plates'\n)
%        disp('the forgoing does not affect computation of capacitance but
%        affects computation of charge density')
kk=0;
kk1=na(1)*nb(1);
kk2=na(2)*nb(2);
kkp=kk1+kk2
for i=1:kkp
    r(i)=0.0;
    v(i)=0.0;
    x(i)=0.0;
    y(i)=0.0;
    z(i)=0.0;
    for j=1:kkp
        ak(i,j)=0.0;
    end
end
 for k=1:npl
        ne(k)=0;
        dx1=(x2(k)-x1(k))/na(k);
        dy1=(y2(k)-y1(k))/na(k);
        dz1=(z2(k)-z1(k))/na(k);
        dx2=(x3(k)-x2(k))/nb(k);
        dy2=(y3(k)-y2(k))/nb(k);
        dz2=(z3(k)-z2(k))/nb(k);
        dx3=(x3(k)-x4(k))/na(k);
        dy3=(y3(k)-y4(k))/na(k);        
        dz3=(z3(k)-z4(k))/na(k);
        dx4=(x4(k)-x1(k))/nb(k);
        dy4=(y4(k)-y1(k))/nb(k);        
        dz4=(z4(k)-z1(k))/nb(k); 
        a(k)=sqrt(dx1^2+dy1^2+dz1^2);
        b(k)=sqrt(dx2^2+dy2^2+dz2^2);

 
    for j=1:nb(k)
        for i=1:na(k)
            ne(k)=ne(k)+1;
            kk=kk+1;
            
            xa=x1(k)+(i-1)*dx1+(j-1)*dx4;
            xb=xa+dx1;
            xc=xa+dx4;
            xd=xb+dx2;
            x(kk)=(xa+xb+xc+xd)/4;
            
            xa=y1(k)+(i-1)*dy1+(j-1)*dy4;
            xb=xa+dy1;
            xc=xa+dy4;
            xd=xb+dy2;
            y(kk)=(xa+xb+xc+xd)/4;
            
            xa=z1(k)+(i-1)*dz1+(j-1)*dz4;
            xb=xa+dz1;
            xc=xa+dz4;
            xd=xb+dz2;
            z(kk)=(xa+xb+xc+xd)/4;
        end
    end
 end
ne;
kk;
a;
b;
x;
y;
z;
 %  The following defines potentials on plates and an index iv used to
 %  treat the potentials
     if npl==1
        vp(1)=vol;
        iv=0;
     else
        if abs(ne(1)*a(1)*b(1)-ne(2)*a(2)*b(2))<(0.01*ne(1)*a(1)*b(1));
%if ne(1)*a(1)*b(1)==ne(2)*a(2)*b(2)
            vp(1)=vol/2;
            vp(2)=-vol/2;
            iv=0;
        else
            vp(1)=vol;
            vp(2)=0.0;
            iv=1;
        end
     end

     err=0.0;
     
   if iv==1
       err=input('enter accuracy required (0.01=1%, 0.001=0.1%, etc.). This is only used with two unequal plates  --->  ')
   end
vp;
      	%pi=3.141592;
      	e0=8.854e-12;
        %e0=1.0;
%
% generate the matrix of coefficients
%
        swi=1;
        q1=0;
%      	qerr=1e-19; %careful - this needs to be small!
      	dv=(vp(1)-vp(2))/1000;


 while swi==1
           	n=ne(1)+ne(2);
            n;
       for i=1:n
        	if i<=ne(1)
               v(i)=vp(1);
            end

            if i>ne(1)
               v(i)=vp(2);
            end

        for j=1:n
        	
                if j<=ne(1);
                     k=1;
                else
                k=2;
                end

                if j==i
                    t1=a(k)*log((b(k)+sqrt(a(k)^2+b(k)^2))/a(k));
                    t2=b(k)*log((a(k)+sqrt(a(k)^2+b(k)^2))/b(k));
                    ak(i,j)=(1/(2*pi*e0))*(t1+t2);
                else
                    ak(i,j)=a(k)*b(k)/(4*pi*e0*sqrt((x(j)-x(i))^2+(y(j)-y(i))^2+(z(j)-z(i))^2));
                end
                end
        end
v;
ak;
    %r=ak\v';
    %r=v/ak;
    
    nn=n-1;
    for i=1:nn
        for m=i+1:n
            fact=ak(m,i)/ak(i,i);
            v(m)=v(m)-v(i)*fact;
            for j=i+1:n
                ak(m,j)=ak(m,j)-ak(i,j)*fact;
            end
        end
    end
    r(n)=v(n)/ak(n,n);
    for i=1:nn
        sum=0.0;
        kkk=n-i;
        for j=kkk+1:n
            sum=sum+ak(kkk,j)*r(j);
        end
        r(kkk)=(v(kkk)-sum)/ak(kkk,kkk);
    end
    r;
    q1=0.0;
    q2=0.0;
    if npl==1
        for i=1:ne(1)
            q1=q1+a(1)*b(1)*r(i);
            swi=0;
        end
    else
            for i=1:ne(1)
                 q1=q1+a(1)*b(1)*r(i);
                 qerr=abs(q1*err);
            end
            for i=ne(1)+1:ne(1)+ne(2)
                q2=q2+a(2)*b(2)*r(i);
            end
            swi=0;
    end
        if iv==1
            dq=q1+q2;
                if abs(dq)<=qerr;
                    swi=0;
                else
                vp(1)=vp(1)-dv;
                vp(2)=vp(2)-dv;
                swi=1;
                end
        end

 end
  vp;
  disp('The following are the charge densities on the individual elements')
  disp('The elements are listed starting at the lower left corner of the plate')
  disp('and ending at the upper right corner of the plate if a single plate is used')
  disp('If two plates are used, then the second plate elements follow the first in the same sequence')
r  
q1
q2
disp(' Capacitance is (in Farads) : ')
 capacitance=q1/vol
