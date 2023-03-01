function [y,param] = eyeArtifact(x,param)


% find blink indicies
bIdx = diff([0 find(sum(abs(x'))==0)]);
param.numberOfBlinks = sum(bIdx~=1);

% define previous window
samples = param.window*param.fs;

% preallocate quality measure
param.quality = zeros(1,size(x,2));

% preallocate output
y = zeros(size(x));

% calculate mean and std dev. of nonzero elements
[ii,~,v] = find(x');
mu = accumarray(ii,v,[],@mean)';
sigma = accumarray(ii,v,[],@std)';

% define lower and upper bound to find invalied samples 
lt = mu-param.factor*sigma;
ht = mu+param.factor*sigma;

% set artifacts to 0 (if less or more than 3 times of mean value it is invalied)
x(x~=0 & x<lt | x>ht) = 0;  

% fill start and ending using mean if zero
x(1:2,x(1,:)==0) = repmat(mu(x(1,:)==0),[2 1]);
x(end-1:end,x(end,:)==0) = repmat(mu(x(end,:)==0),[2 1]);

% loop over columns
for col = 1:size(x,2)
    idx = find(x(:,col)==0);
    if size(idx,1)==1
        bsIdx = [];
    else
        [~,bsIdx] = findpeaks(abs([0; diff(idx)]));
    end
    bs = idx([1; bsIdx])-samples;
    bs(bs<1) = 1;
    be = idx([bsIdx-1; size(idx,1)])+samples;
    be(be>size(x,1)) = size(x,1);
    param.quality(col) = 1-sum(be-bs,1)/size(x,1);
    x(bs:be,col) = 0;
end

% fill start and ending using mean if zero
x(1:2,x(1,:)==0) = repmat(mu(x(1,:)==0),[2 1]);
x(end-1:end,x(end,:)==0) = repmat(mu(x(end,:)==0),[2 1]);

% interpolate missing values
for col = 1:size(y,2)
    b = x(:,col)~=0;
    Y = cumsum(b-diff([1;b])/2);
    y(:,col) = interp1(1:nnz(b),x(b,col),Y,'linear');
end

% smooth eye data
y = medfilt1(y,param.order,[],1);
% y = movmean(y,param.order,1);

% replace first value
y(1,:) = mu;

% if param.debug == 1
%     figure(1);
%     clf;
%     plot(param.time,x(:,1),'LineWidth',2);
%     hold on; grid on;
%     plot(param.time,y(:,1),'--','LineWidth',2);
%     plot(param.time,repelem(ht(1),size(param.time,2)),'k','LineWidth',2);
%     plot(param.time,repelem(lt(1),size(param.time,2)),'k','LineWidth',2);
%     xlabel('Time [sec]','Interpreter','latex');
%     ylabel(param.label(1),'Interpreter','latex');
%     axis tight;
%     ax = gca;
%     ax.TickLabelInterpreter = 'latex';
%     ax.FontSize = 14;
%     ax.YLim = [20 40];
%     ax.LineWidth = 2;
%     legend('Raw','Cleaned','Threshold','Interpreter','latex','Location','north','Orientation','horizontal');
%     
%     print('-f1','C:\Programming\Matlab\adrl-eeg-saeb\rawdata\code\results\figures\fig_eye_artifact','-dpng');
%     
%     figure(2);
%     clf;
%     bar(param.quality);
%     ax = gca;
%     ax.XTickLabel = param.label;
%     ax.XTickLabelRotation = 60;
%     ax.TickLabelInterpreter = 'latex';
%     ax.FontSize = 14;
%     ylabel('Quality index','Interpreter','latex');
%     
% %     print('-f2','C:\Programming\Matlab\adrl-eeg-saeb\rawdata\code\results\figures\fig_eye_quality','-dpng');
% end

end

% eof