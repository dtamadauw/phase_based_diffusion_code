function r = icc(x,y)
% Intraclass Correlation Coefficient
%
% Provided to complement the madicc.m function.  Differs from usual
% Pearson's correlation coefficient only in that the mean and standard
% deviation are pooled (assumed equal) over the two variables.
%
% 2014-07-08
% Thomas Nichols http://warwick.ac.uk/tenichols

I=find(all(~isnan([x(:) y(:)]),2));
if isempty(I)
  r=NaN;
else
  mx    = mean(x(I));
  my    = mean(y(I));
  sx    = std(x(I));
  sy    = std(y(I));
  mn    = (mx+mx)/2;
  sxy   = sqrt((sx^2+sy^2)/2);

  if sxy==0
    r = NaN;
  else
    r = sum((x(I)-mn).*(y(I)-mn)) / sxy^2 / (length(I)-1);
  end

end