% Download the cylinder data from http://dmdbook.com/
load CYLINDER_ALL.mat
p=1;
%%  Compute DMD (Phi are eigenvectors)
VORTALLms = VORTALL - 0*mean(VORTALL,2);
X = VORTALLms(:,1:end-1);
X2 = VORTALLms(:,2:end);
[U,S,V] = svd(X,'econ');

r = 12;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Xproj = U'*X;
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);

eigs = diag(eigs); [~,idx] = sort(abs((eigs)-1));
eigs = eigs(idx); W = W(:,idx);

Phi = X2*V*inv(S)*W;

%% Produce noisy data from many samples
nNoise = 4; % Do not put this higher that 5
noiseMag = 1.25e0;
Xn = zeros(size(VORTALL,1),150*nNoise); Xn2 = Xn;
for jj = 1:nNoise 
%VORTALLN = VORTALLms + noiseMag*randn(size(VORTALL));
VORTALLN = VORTALLms + noiseMag*sqrt(var(VORTALL,[],2)).*randn(size(VORTALL));
Xn(:,150*(jj-1)+(1:150)) = VORTALLN(:,1:end-1);
Xn2(:,150*(jj-1)+(1:150)) = VORTALLN(:,2:end);
end

% Xn = Xn(:,1:5);
% Xn2 = Xn2(:,1:5);

%%
[Un,Sn,Vn] = svd(Xn,'econ');
Un = Un(:,1:r);
Sn = Sn(1:r,1:r);
Vn = Vn(:,1:r);
AtildeN = Un'*Xn2*Vn*inv(Sn);
%AtildeN = U'*Xn2*V*inv(S);
[Wn,eigsNoisy] = eig(AtildeN);

eigsNoisy = diag(eigsNoisy); [~,idx] = sort(abs((eigsNoisy)-1));
eigsNoisy = eigsNoisy(idx); Wn = Wn(:,idx);

Phin = Xn2*Vn/Sn*Wn;

%% Procrustes

%        [Ux,~,~] = svd(X,0); Ux = Ux(:,1:r);
        YprojP = Un'*Xn2; XprojP = Un'*Xn; % Project X and Y onto principal components
    [Uyx, ~, Vyx] = svd(YprojP*XprojP',0);
    AprojP = Uyx*Vyx';    
            [redVecsProc,eigsProc] = eig(AprojP);
        eigsProc = diag(eigsProc);
        PhiProc = Un(:,1:r)*redVecsProc;

%[AtildeProc, eigsProc, PhiProc] = piDMD(Xn,Xn2,'orthogonal',r);

% Procrustes projected
%AtildeProc = piDMD(X1POD,X2POD,'orthogonal');
%[Wproc,eigsProc] = eig(AtildeProc);
%PhiProc = Xn2*Vn*inv(Sn)*Wn;
%magProc = vecnorm(Un'*PhiProc,2,1);
% [~,idx] = sort(abs(eigsProc-1));
% eigsProc = eigsProc(idx); PhiProc = PhiProc(:,idx);
% jp = 5;
f7 = figure(7); LW = 'LineWidth';
plotCylinder(reshape(X(:,1),nx,ny));
set(gca,'position',[0.025 0.025 .95 .95]); % axis off;
set(gca,'TickLength',[0 0],LW,2);yticklabels([]);xticklabels([])
f7.Position(3:4) = 50*[9,4]; f7.Color = 'w';
if p; print('../../images/cylinderPlots/data','-dpng'); end
f8 = figure(8);
plotCylinder(reshape(Xn(:,100),nx,ny)); 
set(gca,'position',[0.025 0.025 .95 .95]); % axis off;
set(gca,'TickLength',[0 0],LW,1);yticklabels([]);xticklabels([])
f8.Position(3:4) = 50*[9,4]; f8.Color = 'w';
if p; print('../../images/cylinderPlots/noisyData','-dpng'); end
% figure(9)
% plotCylinder(reshape(abs(Phin(:,jp))-abs(PhiProc(:,jp)),nx,ny),[],[]); colorbar
%%  Plot  spectrum
LW = 'LineWidth'; FS = 'FontSize';
f1 = figure(19);
theta = (0:1:100)*2*pi/100;
set(gcf, 'renderer', 'painters')
%plot(1i*(theta-pi),'k--')
scatter(imag(log(diag([1; eigs]))),real(log(diag([1;eigs]))),130,'o',LW,1.5,'MarkerEdgeColor',[1 1 1]*.5)
hold on,
scatter(imag(log(diag([1;eigsNoisy]))),real(log(diag([1;eigsNoisy]))),130,'^r',LW,1.5)
scatter(imag(log(diag([1;eigsProc]))),real(log(diag([1;eigsProc]))),130,'xb',LW,1.5)
xline(0); yline(0)
hold off
grid on; ylabel('$\Re[\lambda]$','Interpreter','Latex',FS,25)
xlabel('$\Im[\lambda]$','Interpreter','Latex',FS,25)
set(gca,'TickLabelInterpreter','Latex',FS,15)
axis([-1.2 1.2 -.15 .02]); f1.Color = 'w';
box on;
%set(gca,'position',[0 0 .995 1]); 
f1.Position(3:4) = 3*[200,90];
%set(gca,'TickLength',[0 0],LW,2); box on
if p; print('../../images/cylinderPlots/spectrum','-dpng'); end


% LW = 'LineWidth';
% f1 = figure(1);
% theta = (0:1:100)*2*pi/100;
% set(gcf, 'renderer', 'painters')
% plot(cos(theta),sin(theta),'k--')
% hold on,
% scatter(real(diag(eigs)),imag(diag(eigs)),40,'o',LW,1,'MarkerEdgeColor',[1 1 1]*.5)
% scatter(real(diag(eigsNoisy)),imag(diag(eigsNoisy)),40,'^r',LW,1)
% xline(0); yline(0)
% hold off
%  axis equal; axis(1.3*[-.5 1.3 -.6 .6]);
% set(gca,'position',[0 0 .995 1]); f1.Position(3:4) = [150,100];
% set(gca,'TickLength',[0 0],LW,2); box on
% if p; print('../images/fullExSpecCyli','-dpng'); end
% 
% 
% f2 = figure(2);
% set(gcf, 'renderer', 'painters')
% plot(cos(theta),sin(theta),'k--')
% hold on;
% scatter(real(diag(eigs)),imag(diag(eigs)),40,'o',LW,1,'MarkerEdgeColor',[1 1 1]*.5)
% scatter(real((eigsProc)),imag((eigsProc)),40,'bx',LW,1)
% xline(0); yline(0)
% hold off
%  axis equal; axis(1.3*[-.5 1.3 -.6 .6]);
%  set(gca,'position',[0 0 .995 1]); f2.Position(3:4) = [150,100];
% set(gca,'TickLength',[0 0],LW,2); box on
% if p; print('../images/fullPiSpecCyli','-dpng'); end

%return
%% test reduced order model
m = 199;
procRec = zeros(r,m); Eproc = zeros(1,m);
exactRec = zeros(r,m); Eexact = zeros(1,m);
sg = sign(diag(Xproj/XprojP(:,1:150))); % permute POD modes so they match the noisy ones
Xproj = sg.*Xproj;
procRec(:,1) = Xproj(:,1); Eproc(1) = norm(procRec(:,1))^2;
exactRec(:,1) = Xproj(:,1); Eexact(1) = norm(exactRec(:,1))^2;
for j = 2:m
procRec(:,j) = AprojP*procRec(:,j-1);
exactRec(:,j) = AtildeN*exactRec(:,j-1);
Eproc(j) = norm(procRec(:,j))^2;
Eexact(j) = norm(exactRec(:,j))^2;
end
Etrue = vecnorm(U'*X,2,1).^2;

% figure(21)
% plot(Etrue);
% hold on
% plot(Eproc)
% plot(Eexact)
% hold off

%%
j = 0;
for np = 2:2:10
    j = j+1;
    f = figure(10+j);
    scl = 1;
    sclp= max(abs(Xproj(np,:)));
plot(0:150-1,Xproj(np,:)'/scl,LW,4,'Color',.8*[1 1 1])
hold on
plot(0:m-1,procRec(np,:)'/scl,'b--',LW,2.5)
plot(0:m-1,exactRec(np,:)'/scl,'r--',LW,2.5)
hold off
grid on
xlim([0,70]); ylim(1.1*[-sclp,sclp])%; yticks([-1:.5:1])
%xl = xlim;
%yl = ylim;%
%plot(xl,ones(1,2)*yl(1), '-k',  ones(1,2)*xl(1), yl,'-k', 'LineWidth',1.5)  % Left & Lower Axes
%plot(xl,ones(1,2)*yl(2), '-k',  ones(1,2)*xl(2), yl,'-k', 'LineWidth',1.5)  % Right & Upper Axes
set(gca,'position',[0.1 0.15 .8 .7]); % axis off;
%set(gca,'TickLength',[0 0],LW,.5);
set(gca,'TickLabelInterpreter','Latex',FS,10)
%yticklabels([]);xticklabels([])
f.Position(3:4) = 30*[12,4]; f.Color = 'w';
if p; print(['../../images/cylinderPlots/pod-' num2str(j)],'-dpng'); end
end
%%
% 
% f4 = figure(4);
% imagesc(AtildeN)
% imagesc(AtildeN'*AtildeN)
% for j = 1:r; xline(.5+j); yline(.5+j); end
% caxis(.5*[-1,1]); axis equal; axis tight
% colormap redblue
% set(gca,'position',[0 0 1 1]); f4.Position(3:4) = [200,200];
% set(gca,'TickLength',[0 0],LW,3)
% if p; print('../images/fullPexLinopCyli','-dpng'); end
% 
% f5 = figure(5);
% %imagesc(AtildeProc)
% %imagesc(AtildeProc'*AtildeProc)
% for j = 1:r; xline(.5+j); yline(.5+j); end
% caxis([-1,1]); axis equal; axis tight
% colormap redblue
% set(gca,'position',[0 0 1 1]); f5.Position(3:4) = [200,200];
% set(gca,'TickLength',[0 0],LW,3)
% if p; print('../images/fullPiLinopCyli','-dpng'); end

function plotCylinder(VORT)

%f1 = figure
vortmin = -10;  % only plot what is in -5 to 5 range
vortmax = 10;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

imagesc(VORT); % plot vorticity field
colormap(redblue);
%cmocean('thermal');

% clean up axes
set(gca,'XTick',[1 50 100 150 200 250 300 350 400 449],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8'})
set(gca,'YTick',[1 50 100 150 199],'YTickLabel',{'2','1','0','-1','-2'});
%set(gcf,'Position',[100 100 300 130])
hold on

% add contour lines (positive = solid, negative = dotted)
%contour(VORT,[-5.5:.5:-.5 -.25 -.125],':k','LineWidth',1.2)
%contour(VORT,[.125 .25 .5:.5:5.5],'-k','LineWidth',1.2)
axis equal; axis tight
caxis([-2.5,2.5])
theta = (1:100)/100'*2*pi;
x = 49+25*sin(theta);
y = 99+25*cos(theta);
fill(x,y,[.3 .3 .3])  % place cylinder
plot(x,y,'k','LineWidth',1.2) % cylinder boundary
hold off
end