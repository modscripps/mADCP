function Xsmo = boxcar_smoothNONAN(X,n)

% Xsmo = boxcar_smooth(X,n)
%
% This function performs a 1-D boxcar running average of the variable X,
% smoothing the columns with the length of the boxcar filter n.
% it uses the conv2 function and boxcar,
% the ends of X, which are influenced by the zero-padding performed by
% conv2, are substituted with un-smoothed data from X, so that Xsmo
% is the same length of X and doesn't include any crappy data.

% interpolating over NaN's in X *(between columns),
% then replacing after smoothing.
Inn  = find(isnan(X));
Xi   = NaN.*X;
dumm = 1:size(X,2);

for j=1:size(X,1);
  Ibb = find(~isnan(X(j,:)));
  if length(Ibb)>2;
    Xin=interp1(dumm(Ibb),X(j,(Ibb)),dumm,'linear');
    Xi(j,:)=Xin;
    else
    Xi(j,:)=NaN;
  end
end

% Xi=X;

bb=boxcar(n)./sum(boxcar(n));

Xa=conv2(1,bb',Xi,'same');
Xsmo=Xi.*NaN;
Xsmo(:,ceil(n./2):end-ceil(n./2))=Xa(:,ceil(n./2):end-ceil(n./2));

Xsmo(Inn)=NaN;
