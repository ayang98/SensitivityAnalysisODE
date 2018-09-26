function FDS2
%[t1,xa] = ode15s(@VDP2,(0:.1:20),[2 0]);
%disp(xa);
[a,b] = VDP3(0.1,.2);
end

function dydt = VDP1(t,y) %function for initial u 
global u; %IMPORTANT- in the thesis presentation, the initial u used is 1
dydt= [y(2); u*(1-y(1)^2)*y(2)-y(1)];
end


function dydt = VDP2(t,y) %function for initial u plus delta u
u = 1.5;
dydt= [y(2); u*(1-y(1)^2)*y(2)-y(1)];
end





function [t1,final_y] = VDP3(delta_t,delta_u)
tic
global u;
u=1;
[t1,xa] = ode15s(@VDP1,(0:delta_t:20),[2 0]);
list = zeros(length(t1),2);
values_u = [1, 1.2];
for j=1:length(values_u)
    u = values_u(j);
    [t1,xa] = ode15s(@VDP1,(0:delta_t:20),[2 0]);
    %y_list = [y_list xa(i,1)];
    list(:,j) =  xa(:,1);
end


final_y = (list(:,2)-list(:,1))/delta_u; 

plot(t1, list(:,1));
hold on 
plot(t1, list(:,2));
title('y(p+delta) and y(p)')
ylabel('y')
xlabel('time(t)')
legend('y(p)','y(p+delta)')

figure();
plot(t1, final_y);
title('Sensitivity With Respect to u');
ylabel('$\frac{dy1}{du}$','Interpreter','latex','FontSize',14);
set(get(gca,'ylabel'),'rotation',0)
xlabel('time(t)')

toc

%disp(length(t1))
%disp(length(list(:,1)))
end
