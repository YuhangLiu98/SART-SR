function u = TV_SB_denoising_1D_Lp_Dal(f,mu,phi,p)
if phi == 0
    u = f;
    return;
end
[rows,cols,height]  = size(f);

% Reserve memory for the auxillary variables
u           = f;

x           = zeros(rows,cols,height);
bx          = zeros(rows,cols,height);

sizeI2D = [rows,cols];
fd = [1 0;0 -1]; % Differential operator in diagonal direction
otfFd = psf2otf(fd,sizeI2D);
Denormin2 = abs(otfFd).^2 ;
if height>1
    Denormin2 = repmat(Denormin2,[1,1,height]); % 各层滤波 2D->3D
end
lambda = 2*phi;
while lambda < 1e5
    if p==-inf    % L0范数
        Denormin   = mu + lambda*Denormin2;%分母
        % h-v subproblem2
        h = D1d(u);%[diff(u,1,2), u(:,1,:) - u(:,end,:)]; %列方向一阶差分（后向差分）

        t = h.^2<phi/lambda;
        h(t)=0; 

        % S subproblem1
        Normin2 = D1dt(h);%[h(:,end,:) - h(:, 1,:), -diff(h,1,2)];%列方向一阶差分（前向差分）
        FS = fft2(mu*f + lambda*(Normin2))./Denormin;  %分子偏x偏y
        u = real(ifft2(FS));
        lambda = lambda*2.0;   % β更新
    elseif p == 1    % L1范数   
        dx       = D1d(u);
        x        = shrink1(dx+bx, phi/lambda);

        Denormin = mu + lambda*Denormin2;
        Fu       = fft2(mu*f+lambda*D1dt(x-bx))./Denormin;
        u        = real(ifft2(Fu));      
        bx       = bx+dx-x;

        lambda = lambda*2.0;
    else     % Lp范数 
        dx       = D1d(u);
        x        = shrinkP(dx+bx, phi/lambda,p);

        Denormin = mu + lambda*Denormin2;
        Fu       = fft2(mu*f+lambda*D1dt(x-bx))./Denormin;
        u        = real(ifft2(Fu));      
        bx       = bx+dx-x;
        lambda = lambda*2.0;
    end
end

return

%左上角
function d = D1d(u)
[rows,cols,height] = size(u);
d = zeros(rows,cols,height);
d(2:rows,2:cols,:) = u(2:rows,2:cols,:)-u(1:rows-1,1:cols-1,:);
d(1,1,:) = u(1,1,:)-u(rows,cols,:);
d(1,2:cols,:) = u(1,2:cols,:)-u(rows,1:cols-1,:);
d(2:rows,1,:) = u(2:rows,1,:)-u(1:rows-1,cols,:);
return

function d = D1dt(u)
[rows,cols,height] = size(u);
d = zeros(rows,cols,height);
d(1:rows-1,1:cols-1,:) = u(1:rows-1,1:cols-1,:) - u(2:rows,2:cols,:);
d(rows,cols,:) = u(rows,cols,:) - u(1,1,:);
d(rows,1:cols-1,:) = u(rows,1:cols-1,:) - u(1,2:cols,:);
d(1:rows-1,1,:) = u(1:rows-1,cols,:) - u(2:rows,1,:);
return

function xs = shrink1(x,lambda)

s = abs(x);
xs = sign(x).*max(s-lambda,0);

return;

function xs = shrinkP(x,lambda,p)

s = abs(x);
xs = sign(x).*max(s-(lambda^(2-p)).*(s.^(p-1)),0);
return;

%  % 知识点:
%  % 以x方向的二阶差分为例，fft处理一维赋值成二维，等同于fft2直接处理二维，三维同理；fft、fft2、fftn与otfFx不同
%  uker1 = zeros(1,31);
%  uker1(1,1) = 2; uker1(1,2) = -1; uker1(1,31) = -1; 
%  uker1 = fft(uker1) ;
%  uker1 = repmat(uker1,[31,1]); 
% 
%  uker2 = zeros(31,31);
%  uker2(1,1) = 2; uker2(1,2) = -1; uker2(1,31) = -1; 
%  uker2 = fft2(uker2) ;
%  
%  uker3 = zeros(31,31,31);
%  uker3(1,1,1) = 2; uker3(1,2,1) = -1; uker3(1,31,1) = -1; 
%  uker3    = fftn(uker3);

%  fx = [1, -1]; %Differential operator in x direction
%  otfFx = psf2otf(fx,[31,31]);
%  otfFx = abs(otfFx).^2 ;
%  
%  surf(abs(uker3(:,:,1))-abs(otfFx)); title('corresponding |OTF|');
%  axis square; axis tight
