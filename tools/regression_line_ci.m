function regression_line_ci(x1, y, xrange, varargin)

if length(varargin) > 0
    marker_type = varargin{1};
    marker_color = varargin{2};
    patch_color = varargin{3};
else
    marker_type = 'square';
    marker_color = 'black';
    patch_color = 'k';
end


X = [ones(size(x1)) x1];
[b,bint] = regress(y,X);
disp('Interval slope:')
disp(bint(2,2)-b(2))
disp('Interval intercept:')
disp(bint(1,2)-b(1))
xval = xrange(1):0.01:xrange(2);
yhat = b(1)+b(2)*xval;
ylow = bint(1,1)+bint(2,1)*xval;
yupp = bint(1,2)+bint(2,2)*xval;
plot(x1,y, marker_type, 'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', marker_color);
hold on;
p6=plot(xval,yhat,patch_color,'linewidth',3);
p6.Color(4) = 0.5;
fontSize = 12;
hold on

xpara = [min(xval) max(xval) max(xval) min(xval)]; ypara = [min(ylow) max(ylow)  max(yupp) min(yupp)];
patch(xpara, ypara, patch_color,'EdgeColor','white')
alpha(0.1)
leg=legend('Observed values','Regression line','95% C.I');
set(leg,'location','best')
xlabel('Observed', 'FontSize', fontSize);
ylabel('Predicted', 'FontSize', fontSize);
set(gcf,'color','white')
bias=sum(y-x1)/length(y);
tbl = table(y , x1)
mdl = fitlm(tbl,'linear')
icc_rep = icc(x1,y);
str=[    ' N = ',sprintf('%d',mdl.NumObservations),...
sprintf('\n y = %.2f x + %.2f s^{-1}',b(2),b(1)),...    
sprintf('\n R^2 = %.2f',mdl.Rsquared.Ordinary),...
sprintf('\n ICC = %.3f', icc_rep)
]
annotation('textbox',[.15 0.89 0 0],'string',str,'FitBoxToText','on','EdgeColor','none','fontname', 'Arial', 'FontSize',12, 'FaceAlpha',.0)     
legend1 = legend('boxoff');
set(legend1,...
    'Position',[0.51111111893422 0.159027781503068 0.36874999217689 0.135416662941376]);