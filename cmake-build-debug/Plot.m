clear all
clc
close all
M = readmatrix('Integrazione.txt');
C = readmatrix('Control_vector.txt');
X = readmatrix('Snake.txt');

figure()
plot3(M(:,13),M(:,12),M(:,11))
% xlim([-500 500]);
% ylim([-500 500]);
 zlim([50 200]);
 xlabel('E[m]')
 ylabel('N[m]')
 zlabel('h[m]')
 grid on

figure()
plot(M(:,13),M(:,12))
% axis equal
xlabel('E[m]')
ylabel('N[m]')
% xlim([-5000 5000]);
% ylim([-500 500]);
grid on

figure()
plot(M(:,1),M(:,11))
axis equal
% xlim([-5000 5000]);
% ylim([-500 500]);
grid on

figure()
subplot(2,2,1)
plot(C(:,1),C(:,2))
legend('De[°]')

subplot(2,2,2)
plot(C(:,1),C(:,3))
legend('Da[°]')

subplot(2,2,3)
plot(C(:,1),C(:,4))
legend('Dth[%]')

figure()
plot(M(:,1),M(:,2))
hold on
plot(M(:,1),M(:,3))
plot(M(:,1),M(:,4))
legend('u[m/s]','v[m/s]','w[m/s]')

figure()
plot(M(:,1),M(:,5))
hold on
plot(M(:,1),M(:,6))
plot(M(:,1),M(:,7))
legend('p[rad/s]','q[rad/s]','r[rad/s]')

figure()
plot(M(:,1),M(:,8))
hold on
plot(M(:,1),M(:,9))
plot(M(:,1),M(:,10))
plot(X(:,1)-20.5,X(:,2))
legend('phi[rad]','theta[rad]','psi[rad]','psi_{ref} [rad]')

figure()
subplot(2,2,1)
plot(M(:,1),M(:,14))
legend('V[rad]')

subplot(2,2,2)
plot(M(:,1),M(:,15))
legend('alpha[rad]')

subplot(2,2,3)
plot(M(:,1),M(:,16))
legend('beta[rad]')
