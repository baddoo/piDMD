%% Pendulum example
%  This example applies energy-preserving piDMD to reconstruct the trajectory of
%  a double pendulum from noisy measurements.
addpath('../src')
rng(1); % Set random seed

% Define parameters of problem
l1=1; l2=1.5; % Lengths of rods
m1=1 ; m2=1.5; g=9.81; % Masses and gravity
params = [l1, l2, m1, m2, g]; % Concatenate parameters

% Construct linearised energy inner product
W(1,1) = (m1/2+m2/2)*g*l1;
W(2,2) = m2/2*g*l2;
W(3,3) = (m1/2+m2/2)*l1^2;
W(4,3) = m2/2*l1*l2;
W(3,4) = W(4,3);
W(4,4) = m2/2*l2^2;
C = chol(W); % Calculate inner product

% Set number of samples and span
tend = 30; nt = 1000;
tspan= linspace(0,tend,nt);

% Set initial conditions
theta1= 0.4; theta1_prime=0;
theta2= 0.7; theta2_prime=0;
y0=[theta1 theta1_prime theta2 theta2_prime];

% Solve ODE
[t,y]=ode45(@(t,y) pendulum(t,y,params), tspan,y0);

% Extract data
th1 = y(:,1); th2 = y(:,3); th1dt = y(:,2); th2dt = y(:,4);
x = [th1'; th2'; th1dt'; th2dt'];
xn = x + 1e-1*std(x,[],2).*randn(size(x)); % Add noise
data = C*xn; % Rescale measurements into energy norm
nTrain = nt-1;
X = data(:,1:nTrain); Y = data(:,2:nTrain+1);

% Train models
[piA,piVals] = piDMD(X,Y,'orthogonal');
[exA,exVals] = piDMD(X,Y,'exact');

% Perform reconstructions
piRec = zeros(4,nt); piRec(:,1) = data(:,1);
exRec = zeros(4,nt); exRec(:,1) = data(:,1);
for j = 2:nt
piRec(:,j) = piA(piRec(:,j-1));
exRec(:,j) = exA(exRec(:,j-1));
end

% Rescale reconstructions back into physical norm
fpiRec = C\piRec; fexRec = C\exRec;

%% Plot results
figure(1); LW = 'LineWidth'; IN = 'Interpreter'; LT = 'Latex'; FS = 'FontSize';
subplot(3,1,1)
c1 = .8*[1 1 1]; c2 = .8*[1 1 1];
plot(tspan,xn(1,:),LW,2,'Color',c1)
hold on
plot(tspan,xn(2,:),LW,2,'Color',c2)
hold off; trajPlot('measurements'); xticklabels([])
subplot(3,1,2)
plot(tspan,x(1,:),LW,3,'Color', c1)
hold on
plot(tspan,fpiRec(1,:),'b--',LW,2)
plot(tspan,fexRec(1,:),'r--',LW,2)
ylabel('$\theta_1$',IN,LT)
hold off; trajPlot('$\theta_1$'); xticklabels([])
subplot(3,1,3)
l1=plot(t,x(2,:),LW,3,'Color', c2);
hold on
l2=plot(t,fpiRec(2,:),'b--',LW,2);
l3=plot(t,fexRec(2,:),'r--',LW,2);
hold off; xlabel('time',FS,20,IN,LT)
trajPlot('$\theta_2$')
legend([l1,l2,l3],{'truth','piDMD','exact DMD'},IN,LT)
function yp = pendulum(~, y, params)

l1=params(1);  l2=params(2); 
m1=params(3);  m2=params(4); 
g=params(5);

a = (m1+m2)*l1 ;
b = m2*l2*cos(y(1)-y(3)) ;
c = m2*l1*cos(y(1)-y(3)) ;
d = m2*l2 ;
e = -m2*l2*y(4)* y(4)*sin(y(1)-y(3))-g*(m1+m2)*sin(y(1)) ;
f = m2*l1*y(2)*y(2)*sin(y(1)-y(3))-m2*g*sin(y(3)) ;
yp=zeros(4,1);
yp(1) = y(2);
yp(3)= y(4) ;
yp(2)= (e*d-b*f)/(a*d-c*b) ;
yp(4)= (a*f-c*e)/(a*d-c*b) ;
end

function trajPlot(j) % Nice plot of trajectories
yticks([-pi/4,0,pi/4]); yticklabels([{'$-\pi/4$'},{'0'},{'$\pi/4$'}])
set(gca,'TickLabelInterpreter','Latex','FontSize',20);grid on
ylim([-1,1])
ylabel(j,'Interpreter','latex','FontSize',20)
end