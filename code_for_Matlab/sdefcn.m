function dydt = sdefcn(t,y,A,B,C)     
global beta1 beta2 beta3 mu11 mu12 mu13 mu14 mu21 mu22 mu23 mu24 mu31 mu32 mu33 mu34 
global ktt11 ktt12 ktt13 ktt14 ktt21 ktt22 ktt23 ktt24 ktt31 ktt32 ktt33 ktt34

gama1=0.1*rand;gama2=0.1*rand;gama3=0.1*rand;gama4=0.1*rand; % degradation rate of TGs

dydt = zeros(7,1);
dydt(1)=A-beta1.*y(1);
dydt(2)=B-beta2.*y(2);
dydt(3)=C-beta3.*y(3);
        
dydt(4)=mu11*y(1)./(ktt11+y(1))+mu21*y(2)./(ktt21+y(2))+mu31*y(3)./(ktt31+y(3))-gama1*y(4);
dydt(5)=mu12*y(1)./(ktt12+y(1))+mu22*y(2)./(ktt22+y(2))+mu32*y(3)./(ktt32+y(3))-gama2*y(5);
dydt(6)=mu13*y(1)./(ktt13+y(1))+mu23*y(2)./(ktt23+y(2))+mu33*y(3)./(ktt33+y(3))-gama3*y(6);
dydt(7)=mu14*y(1)./(ktt14+y(1))+mu24*y(2)./(ktt24+y(2))+mu34*y(3)./(ktt34+y(3))-gama4*y(7);

