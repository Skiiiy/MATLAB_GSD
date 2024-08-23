close all
clear
clc

MT=readmatrix("NEWtable8.3.xlsx","Sheet","Sheet1","Range",'mach');        % 马赫数
AT=readmatrix("NEWtable8.3.xlsx","Sheet","Sheet1","Range",'A');           % A面积函数

% 初始波阵面定义  波阵面分辨率：a0<=0.01;CFL条件：ndt<0.1，/最大马赫数
jmax=500;imax=201;ndt=0.1;
xmin=-0.01;xmax=10;
ymin=1;ymax=6;
nx0=1;ny0=0;
gamma=1.4;pkz=11;
% m0=sqrt(1+(gamma+1)./(2*gamma).*pkz);
m0=2.9326;
c0=1;
a0=(ymax-ymin)/(imax-1);dt=ndt*a0;T=(jmax-1)*dt;

% 壁面定义 凹theta+，凸theta-
nnn=1000;nn=10000;n=20;
yuan=0.0001;

r=0.5;
xA=linspace(-1,0.5-r-yuan*a0,nn);
yA=linspace(1,1,nn);
AxB=linspace(0.5-r+yuan*a0,0.5+r-yuan*a0,nn);
AyB=sqrt(r.^2-(0.5-AxB).^2)+1;
Bx=linspace(0.5+r+yuan*a0,7,nn*10);
By=linspace(1,1,nn*10);
ox=[xA,AxB,Bx];
oy=[yA,AyB,By];

% theta=-89.9;
% xA=linspace(-1,1-yuan*a0,nn);
% yA=linspace(1,1,nn);
% Ax=linspace(1+yuan*a0*cosd(theta),5,nn);
% Ay=linspace(1-yuan*a0*sind(theta),1+(5-1)*tand(theta),nn);
% ox=[xA,Ax];
% oy=[yA,Ay];

[ox,idox]=unique(ox);
oy=oy(idox);
oxx=linspace(min(ox),max(ox),nn);
YYY=spline(ox,oy);
oyy=ppval(YYY,oxx);

% XYY=[ox;oy];
% tt=1:length(XYY);
% yyy=spline(tt,[[ox(2)-ox(1);oy(2)-oy(1)],XYY,[ox(end)-ox(end-1);oy(end)-oy(end-1)]]);
% pxxx=ppval(yyy,linspace(tt(1),tt(end),nn));

% figure()
% subplot(1,1,1),
% % plot(pxxx(1,:),pxxx(2,:),'linewidth',1)
% plot(oxx,oyy) % title('波阵面位置')
% xlim([-1 7]),ylim([-3 5]),daspect([1 1 1])
% hold on

% 定义各节点坐标矩阵
px=zeros(imax,jmax);py=zeros(imax,jmax);
% 定义各节点速度方向单位向量矩阵
nx=zeros(imax,jmax);ny=zeros(imax,jmax);
% 截面积
A=zeros(imax,jmax);
% 马赫数&实际速度
M=zeros(imax,jmax);V=c0.*M;
% 离散弧长
S=zeros(imax,jmax);
% 样条函数
x0=zeros(1,jmax);y0=zeros(1,jmax);
xn=zeros(1,jmax);yn=zeros(1,jmax);
% 各管常数值
k=zeros(imax,jmax);k2=zeros(imax,jmax);
% 加减点判断数
sigma=zeros(imax,jmax);

%% j=1,2

A(:,1)=a0;M(:,1)=m0;
px(:,1)=xmin;py(:,1)=(ymin:a0:ymax).';
nx(:,1)=nx0;ny(:,1)=ny0;
% plot(px(:,1),py(:,1),'linewidth',1)
% drawnow;
% plot(px(:,1),py(:,1),'.')

px(:,2)=px(:,1)+M(:,1).*nx(:,1).*dt;
py(:,2)=py(:,1)+M(:,1).*ny(:,1).*dt;

a_min=globalsearch(px(2,2),py(2,2));
px(1,2)=a_min;
py(1,2)=ppval(YYY,a_min);
% plot(px(:,2),py(:,2),'linewidth',1)
% plot(px(:,2),py(:,2),'.')

for i = 2:imax
    S(i,1:2)=S(i-1,1:2)+sqrt((px(i,1:2)-px(i-1,1:2)).^2+((py(i,1:2)-py(i-1,1:2)).^2));
    sigma(i,1:2)=(S(i,1:2)-S(i-1,1:2))./a0;
end

% 拟合样条函数用于法向量计算
y0(1,2)=(py(2,2)-py(1,2))./(S(2,2)-S(1,2));
x0(1,2)=(px(2,2)-px(1,2))./(S(2,2)-S(1,2));

yn(1,2)=(py(imax,2)-py(imax-1,2))./(S(imax,2)-S(imax-1,2));
xn(1,2)=(px(imax,2)-px(imax-1,2))./(S(imax,2)-S(imax-1,2));

Y=spline(S(:,2).',[y0(1,2) py(:,2).' yn(1,2)]);
X=spline(S(:,2).',[x0(1,2) px(:,2).' xn(1,2)]);

DY(:,2)=((ppval(Y,S(:,2)+0.01.*a0))-(ppval(Y,S(:,2)-0.01.*a0)))./(0.02.*a0);
DX(:,2)=((ppval(X,S(:,2)+0.01.*a0))-(ppval(X,S(:,2)-0.01.*a0)))./(0.02.*a0);
D(:,2)=sqrt((DX(:,2)).^2+(DY(:,2)).^2);

nx(:,2)=DY(:,2)./D(:,2);
ny(:,2)=-DX(:,2)./D(:,2);

for i = 2:imax
    if i>imax-1
        A(i,1:2)=0.5.*(S(i,1:2)-S(i-1,1:2));
    else
        A(i,1:2)=0.5.*((S(i+1,1:2)-S(i-1,1:2)));
    end
end
A(1,1:2)=0.5.*(S(2,1:2)-S(1,1:2));

% GSD
for i = 1:imax
    k(i,1)=interp1(MT,AT,M(i,2-1),"linear")./A(i,2-1);
    aaa=k(i,1).*A(i,2);
    M(i,2)=interp1(AT,MT,aaa,"linear");
    V(i,2)=M(i,2).*c0;
end

%% j>=3 蛙跳格式或重新计算 ttt=1 点数变化之后；ttt=2 两步蛙跳，点数不变
ttt=0;
for j = 3:jmax
    DDD=0;JJJ=0;mmm=0;


    if ttt==1  % 重新按照初始计算
        px(:,j)=px(:,j-1)+M(:,j-1).*nx(:,j-1).*dt;
        py(:,j)=py(:,j-1)+M(:,j-1).*ny(:,j-1).*dt;

        a_min=globalsearch(px(2,j),py(2,j));
        px(1,j)=a_min;
        py(1,j)=ppval(YYY,a_min);

        for i = 2:imax
            S(i,j)=S(i-1,j)+sqrt((px(i,j)-px(i-1,j)).^2+((py(i,j)-py(i-1,j)).^2));
            sigma(i,j)=(S(i,j)-S(i-1,j))./a0;
        end

        for i = 2:imax
            if i>imax-1
                A(i,j)=0.5.*(S(i,j)-S(i-1,j));
            else
                A(i,j)=0.5.*((S(i+1,j)-S(i-1,j)));
            end
        end
        A(1,j)=0.5.*(S(2,j)-S(1,j));

        % GSD
        for i = 1:imax
            k(i,j-1)=interp1(MT,AT,M(i,j-1),"linear")./A(i,j-1);
            aaa=k(i,j-1).*A(i,j);
            M(i,j)=interp1(AT,MT,aaa,"linear");
            V(i,j)=M(i,j).*c0;
        end

    else    % 按照两步蛙跳计算
        px(:,j)=px(:,j-2)+2.*M(:,j-1).*nx(:,j-1).*dt;
        py(:,j)=py(:,j-2)+2.*M(:,j-1).*ny(:,j-1).*dt;

        a_min=globalsearch(px(2,j),py(2,j));
        px(1,j)=a_min;
        py(1,j)=ppval(YYY,a_min);

        for i = 2:imax
            S(i,j)=S(i-1,j)+sqrt((px(i,j)-px(i-1,j)).^2+((py(i,j)-py(i-1,j)).^2));
            sigma(i,j)=(S(i,j)-S(i-1,j))./a0;
        end

        for i = 2:imax
            if i>imax-1
                A(i,j)=0.5.*(S(i,j)-S(i-1,j));
            else
                A(i,j)=0.5.*((S(i+1,j)-S(i-1,j)));
            end
        end
        A(1,j)=0.5.*(S(2,j)-S(1,j));

        % GSD
        for i = 1:imax
            k(i,j-2)=interp1(MT,AT,M(i,j-2),"linear")./A(i,j-2);
            aaa=k(i,j-2).*A(i,j);
            M(i,j)=interp1(AT,MT,aaa,"linear");
            V(i,j)=M(i,j).*c0;
        end

    end

    % 拟合样条函数用于加点插值
    SX=[S(:,j).';px(:,j).'];
    sx=spline(S(:,j),[[S(2,j)-S(1,j);px(2,j)-px(1,j)],SX,[0;1]]);

    SY=[S(:,j).';py(:,j).'];
    sy=spline(S(:,j),[[S(2,j)-S(1,j);py(2,j)-py(1,j)],SY,[1;1]]);

    % 加点
    for i = 2:imax
        if sigma(i,j)>1.5
            DDD=DDD+1;
            inputS=0.5*(S(i,j)+S(i-1,j));
            inputxx=ppval(sx,inputS);
            inputyy=ppval(sy,inputS);
            px(:,j)=[px(1:i-1+(DDD-1),j);inputxx(2,1);px(i+(DDD-1):imax-1,j)];
            py(:,j)=[py(1:i-1+(DDD-1),j);inputyy(2,1);py(i+(DDD-1):imax-1,j)];
            M(:,j)=[M(1:i-1+(DDD-1),j);M(i-1+(DDD-1),j);M(i+(DDD-1):imax-1,j)];
        end
    end

    % 减点情况1：波阵面点运动到壁面范围内
    for i = 1:imax
        if py(i,j)<ppval(YYY,px(i,j)) && ppval(YYY,px(i,j))<py(i+1,j)
            JJJ=JJJ+1;
            px(:,j)=[px(i+1:imax,j);px(imax,j).*ones(i,1)];
            py(:,j)=[py(i+1:imax,j);(py(imax,j).*ones(i,1)+a0.*linspace(1,i,i).')];
            M(:,j)=[M(i+1:imax,j);M(imax,j).*ones(i,1)];
            for u=2:imax
                S(u,j)=S(u-1,j)+sqrt((px(u,j)-px(u-1,j)).^2+((py(u,j)-py(u-1,j)).^2));
                sigma(u,j)=(S(u,j)-S(u-1,j))./a0;
            end
        end
    end

    % 减点情况2：压缩波导致单个空间步长过小   
    for i = 2:imax
        if sigma(i,j)<0.5
            JJJ=JJJ+1;
            px(:,j)=[px(1:i-1,j);px(i+1:imax,j);px(imax,j)];
            py(:,j)=[py(1:i-1,j);py(i+1:imax,j);py(imax,j)+a0];
            M(:,j)=[M(1:i-1,j);M(i+1:imax,j);M(imax,j)];
            for u=2:imax
                S(u,j)=S(u-1,j)+sqrt((px(u,j)-px(u-1,j)).^2+((py(u,j)-py(u-1,j)).^2));
                sigma(u,j)=(S(u,j)-S(u-1,j))./a0;
            end
        end
    end

    % 平滑程序（分奇、偶点）
    if mod(j,n)==0
        mmm=mmm+1;
        for i = 2:imax-1
            if mod(i,2)==0
                px(i,j)=0.5*(px(i+1,j)+px(i-1,j));
                py(i,j)=0.5*(py(i+1,j)+py(i-1,j));
                M(i,j)=0.5*(M(i+1,j)+M(i-1,j));
            end
        end
        for i = 2:imax-1
            if mod(i+1,2)==0
                px(i,j)=0.5*(px(i+1,j)+px(i-1,j));
                py(i,j)=0.5*(py(i+1,j)+py(i-1,j));
                M(i,j)=0.5*(M(i+1,j)+M(i-1,j));
            end
        end
    end

    a_min=globalsearch(px(2,j),py(2,j));
    px(1,j)=a_min;
    py(1,j)=ppval(YYY,a_min);

    for i = 2:imax
        S(i,j)=S(i-1,j)+sqrt((px(i,j)-px(i-1,j)).^2+((py(i,j)-py(i-1,j)).^2));
    end

    for i = 2:imax
        if i>imax-1
            A(i,j)=0.5.*(S(i,j)-S(i-1,j));
        else
            A(i,j)=0.5.*((S(i+1,j)-S(i-1,j)));
        end
    end
    A(1,j)=0.5.*(S(2,j)-S(1,j));

    % 重新拟合样条函数用于法向量计算
    y0(1,j)=(py(2,j)-py(1,j))./(S(2,j)-S(1,j));
    x0(1,j)=(px(2,j)-px(1,j))./(S(2,j)-S(1,j));

    yn(1,j)=(py(imax,j)-py(imax-1,j))./(S(imax,j)-S(imax-1,j));
    xn(1,j)=(px(imax,j)-px(imax-1,j))./(S(imax,j)-S(imax-1,j));

    Y=spline(S(:,j).',[y0(1,j) py(:,j).' yn(1,j)]);
    X=spline(S(:,j).',[x0(1,j) px(:,j).' xn(1,j)]);

    DY(:,j)=((ppval(Y,S(:,j)+0.01.*a0))-(ppval(Y,S(:,j)-0.01.*a0)))./(0.02.*a0);
    DX(:,j)=((ppval(X,S(:,j)+0.01.*a0))-(ppval(X,S(:,j)-0.01.*a0)))./(0.02.*a0);
    D(:,j)=sqrt((DX(:,j)).^2+(DY(:,j)).^2);

    nx(:,j)=DY(:,j)./D(:,j);
    ny(:,j)=-DX(:,j)./D(:,j);

    if DDD+JJJ+mmm>0
        ttt=1;
    else
        ttt=2;
    end

    outputj=j
end

plotj=280;

figure()
subplot(1,1,1),
% plot(pxxx(1,:),pxxx(2,:),'linewidth',1)
plot3(oxx,oyy,ones(1,nn).*10,"k",'linewidth',1) % title('波阵面位置')
xlim([-0.5 5]),ylim([0.5 5]),daspect([1 1 1]),box on
hold on
% plot3(px(:,1),py(:,1),ones(imax,1).*10,"w",'linewidth',1)
% subplot(1,2,2)
s=surf(px(:,1:plotj),py(:,1:plotj),M(:,1:plotj),'FaceAlpha',1,'FaceColor','interp'),s.EdgeColor = 'none';
% xlim([-0.5 5]),ylim([-0.5 2]),daspect([1 1 1]),
view(2),colormap jet

for j=1:plotj
    XY=[px(:,j).';py(:,j).'];
    t=1:length(XY);
    yy=spline(t,[[px(2,j)-px(1,j);py(2,j)-py(1,j)],XY,[px(imax,j)-px(imax-1,j);py(imax,j)-py(imax-1,j)]]);
    pxx=ppval(yy,linspace(t(1),t(end),nn));
    if mod(j,20)==0
        plot3(pxx(1,:),pxx(2,:),ones(1,nn).*10,"w",'linewidth',1)
    end
end

f = gcf;
exportgraphics(f,'123456789.png','Resolution',600)