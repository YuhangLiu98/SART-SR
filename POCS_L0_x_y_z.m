function [ f,errorL2,qualMeasOut ]  = POCS_L0_x_y_z(head,proj,geo,angles,maxiter,smooth_lambda,smooth_normType,u,varargin)
%   'lambdared':   Reduction of lambda.Every iteration
%                  lambda=lambdared*lambda. Default is 0.99
%
%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20
%
%   'alpha':       Defines the TV hyperparameter. default is 0.002
%
%   'alpha_red':   Defines the reduction rate of the TV hyperparameter
%
%   'Ratio':       The maximum allowed image/TV update ration. If the TV
%                  update changes the image more than this, the parameter
%                  will be reduced.default is 0.95
%   'maxL2err'     Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%
% 'OrderStrategy'  Chooses the subset ordering strategy. Options are
%                  'ordered' :uses them in the input order, but divided
%                  'random'  : orders them randomply
%                  'angularDistance': chooses the next subset with the
%                                     biggest angular distance with the ones used.
%--------------------------------------------------------------------------

%% parse inputs
blocksize=1;
[beta,beta_red,ng,verbose,alpha,alpha_red,rmax,epsilon,OrderStrategy,QualMeasOpts,nonneg]=parse_inputs(proj,geo,angles,varargin);
measurequality=~isempty(QualMeasOpts);

[alphablocks,orig_index]=order_subsets(angles,blocksize,OrderStrategy);

angles_reorder=cell2mat(alphablocks);
index_angles=cell2mat(orig_index);

% does detector rotation exists?
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end
%% Create weigthing matrices for the SART step
% the reason we do this, instead of calling the SART fucntion is not to
% recompute the weigths every ASD-POCS iteration, thus effectively doubling
% the computational time
% Projection weigth, W


geoaux=geo;
geoaux.sVoxel([1 2])=geo.sVoxel([1 2])*1.1; % a Bit bigger, to avoid numerical division by zero (small number)
geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'ray-voxel');  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;


% Back-Projection weigth, V
V=computeV(geo,angles,alphablocks,orig_index);


% initialize image.
f=zeros(geo.nVoxel','single');

%%
stop_criteria=0;
iter=0;
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
rotDetector=geo.rotDetector;
DSD=geo.DSD;
DSO=geo.DSO;
errorL2 = [];
while ~stop_criteria %POCS
    f0=f;
    if (iter==0 && verbose==1)
        tic;
    end
    iter=iter+1;
    
    for jj=1:size(angles,2)
        if size(offOrigin,2)==size(angles,2)
            geo.offOrigin=offOrigin(:,index_angles(:,jj));
        end
        if size(offDetector,2)==size(angles,2)
            geo.offDetector=offDetector(:,index_angles(:,jj));
        end
        if size(rotDetector,2)==size(angles,2)
            geo.rotDetector=rotDetector(:,index_angles(:,jj));
        end
        if size(DSD,2)==size(angles,2)
            geo.DSD=DSD(jj);
        end
        if size(DSO,2)==size(angles,2)
            geo.DSO=DSO(jj);
        end
        %         proj_err=proj(:,:,jj)-Ax(f,geo,angles(:,jj));          %                                 (b-Ax)
        %         weighted_err=W(:,:,jj).*proj_err;                   %                          W^-1 * (b-Ax)
        %         backprj=Atb(weighted_err,geo,angles(:,jj));            %                     At * W^-1 * (b-Ax)
        %         weigth_backprj=bsxfun(@times,1./V(:,:,jj),backprj); %                 V * At * W^-1 * (b-Ax)
        %         f=f+beta*weigth_backprj;                          % x= x + lambda * V * At * W^-1 * (b-Ax)
        % Enforce positivity
        f=f+beta* bsxfun(@times,1./V(:,:,index_angles(:,jj)),Atb(W(:,:,index_angles(:,jj)).*(proj(:,:,index_angles(:,jj))-Ax(f,geo,angles_reorder(:,jj))),geo,angles_reorder(:,jj)));
        % non-negativity constrain
        if nonneg
            f=max(f,0);
        end
    end
    
    geo.offDetector=offDetector;
    geo.offOrigin=offOrigin;
    geo.DSD=DSD;
    geo.DSO=DSO;
    geo.rotDetector=rotDetector;
    % Save copy of image.

    % compute L2 error of actual image. Ax-b
    dd=im3Dnorm(Ax(f,geo,angles)-proj,'L2');
    % compute change in the image after last SART iteration
    dp_vec=(f-f0);
    dp=im3Dnorm(dp_vec,'L2');
    
    if iter==1
        dtvg=alpha*dp;
        %Convert the steepest-descent step-size from a fraction of a
        %step-size to an absolute image distance on the first iteration.
    end
    f0=f;
    
    % =========================================================================
    f0 = permute(f0,[1 3 2]); % [x,z,y]->[x,y,z]
    f1 = TV_SB_denoising_1D_Lp(f0,1,smooth_lambda(1),smooth_normType(1));     % x-direction
    f1 = permute(f1,[2 1 3]); % [x,y,z]->[y,x,z]
    f2 = TV_SB_denoising_1D_Lp(f1,1,smooth_lambda(2),smooth_normType(2));     % y-direction
    f2 = permute(f2,[2 1 3]); % [y,x,z]->[x,y,z]
    f3 = TV_SB_denoising_1D_Lp_Dal(f2,1,smooth_lambda(3),smooth_normType(3)); % left diagonal
    f3 = rot90(f3); % 逆时针旋转90°，以便在另一个对角线操作
    f4 = TV_SB_denoising_1D_Lp_Dal(f3,1,smooth_lambda(4),smooth_normType(4)); % right diagonal
    f4 = rot90(f4,3); %再次逆时针旋转270°，恢复成原样
    f4 = permute(f4,[1 3 2]); % [x,y,z]->[x,z,y] 
    f = TV_SB_denoising_1D_Lp(f4,1,smooth_lambda(5),smooth_normType(5));      % z-direction
    f0 = permute(f0,[1 3 2]);
    f = single(f);
    
    f=minimizeTAwTV(f,dtvg,ng,u); %   TAwTV, CUDA is used
%     f=minimizeTV(f,dtvg,ng);    %   TV,    CUDA is used
    % ==========================================================================
    % compute change by TAwTV/TV min
    dg_vec=(f-f0);
    dg=im3Dnorm(dg_vec,'L2');
    % if change in TAwTV/TV is bigger than the change in SART AND image error is still bigger than acceptable
    if dg>rmax*dp && dd>epsilon
        dtvg=dtvg*alpha_red;
    end
    % reduce SART step
    beta=beta*beta_red;
    % Check convergence criteria
    % ==========================================================================
    
    %Define c_alpha as in equation 21 in the journal
    c=dot(dg_vec(:),dp_vec(:))/(norm(dg_vec(:),2)*norm(dp_vec(:),2));
    %This c is examined to see if it is close to -1.0
    
    if iter>=maxiter
%         if verbose
%             disp(['Stopping criteria met']);
%             disp(['   c    = ' num2str(c),      '(Desired: c<-0.99)']);
%             disp(['   beta = ' num2str(beta),   '(Desired: beta<0.005)']);
%             disp(['   iter = ' num2str(iter), ]);
%         end
        stop_criteria=true;
    end

    if (iter==1 && verbose==1)
        expected_time=toc*maxiter;
        disp('ADS-POCS');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    if measurequality
        qualMeasOut(:,iter)=Measure_Quality(head,f,QualMeasOpts);
    end
    errornow = norm(f(:)-head(:))/norm(head(:));
    errorL2=[errorL2 errornow];
end



end

function [beta,beta_red,ng,verbose,alpha,alpha_red,rmax,epsilon,OrderStrategy,QualMeasOpts,nonneg]=parse_inputs(proj,geo,angles,argin)

opts=     {'lambda','lambda_red','tviter','verbose','alpha','alpha_red','ratio','maxl2err','orderstrategy','qualmeas','nonneg'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:ASD_POCS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:ASD_POCS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extract value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('TIGRE:ASD_POCS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
        end
        val=argin{jj};
    end
    % parse inputs
    switch opt
        % Verbose
        %  =========================================================================
        case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE:Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
            % Lambda
            %  =========================================================================
            % Its called beta in ASD-POCS
        case 'lambda'
            if default
                beta=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('CBCT:ASD_POCS:InvalidInput','Invalid lambda')
                end
                beta=val;
            end
            % Lambda reduction
            %  =========================================================================
        case 'lambda_red'
            if default
                beta_red=0.99;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid lambda')
                end
                beta_red=val;
            end
            % Number of iterations of TV
            %  =========================================================================
        case 'tviter'
            if default
                ng=20;
            else
                ng=val;
            end
            %  TV hyperparameter
            %  =========================================================================
        case 'alpha'
            if default
                alpha=0.002; % 0.2
            else
                alpha=val;
            end
            %  TV hyperparameter redution
            %  =========================================================================
        case 'alpha_red'
            if default
                alpha_red=0.95;
            else
                alpha_red=val;
            end
            %  Maximum update ratio
            %  =========================================================================
        case 'ratio'
            if default
                rmax=0.95;
            else
                rmax=val;
            end
            %  Maximum L2 error to have a "good image"
            %  =========================================================================
        case 'maxl2err'
            if default
                epsilon=im3Dnorm(FDK(proj,geo,angles),'L2')*0.2; %heuristic
            else
                epsilon=val;
            end
        case 'orderstrategy'
            if default
                OrderStrategy='random';
            else
                OrderStrategy=val;
            end
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:ASD_POCS:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'nonneg'
            if default
                nonneg=true;
            else
                nonneg=val;
            end
        otherwise
            error('TIGRE:ASD_POCS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in ASD_POCS()']);
            
    end
end

end






