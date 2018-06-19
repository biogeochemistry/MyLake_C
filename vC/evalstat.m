%Funktiota, joka laskee virhesysteemit.
%Sen lis‰ksi, ett‰ biasin merkki on v‰‰rinp‰in, n_obs on v‰‰r‰ silloin kun
%havainnoissa on NaN-rivej‰, ja kaikki kaavat, joissa on n_obs, antavat
%v‰‰rin!
% -Jos t‰h‰n saisi korjattua...

%Lis‰tty bias-output 16.7.2015

% A = observation, B = model result
function [R2, rmse, NSE, Pbias, Nbias, rmsenu, bias] = evalstat(A,B)

%A =[4,5,6,7,8;4,5,NaN,7,8]';B =[3 6 5 4 6;3 6 5 4 6]';
n_obs = sum(~isnan(A));
n = size(A,2);
Aav = nansum(A)./n_obs;
%Myˆs mallin keskiarvoa laskettaessa tulee poistaa mallitulokset, joita
%vastaavaa havaintoa ei ole!
Bav = NaN*zeros(1,n);
%keyboard
for i=1:n
    inxx = find(~isnan(A(:,i)));
    Bav(1,i) = nansum(B(inxx,i))./n_obs(i);
end

A_av = NaN*zeros(size(A));
B_av = NaN*zeros(size(A));
for i=1:n
    inxx = find(~isnan(A(:,i)));
    A_av(inxx,i) = A(inxx,i)-Aav(i);
    B_av(inxx,i) = B(inxx,i)-Bav(i);
end

A_av2 = A_av.*A_av;
B_av2 = B_av.*B_av;

A_av2sum = nansum(A_av2);
B_av2sum = nansum(B_av2);

B_A = B-A;
A_B2 = B_A.*B_A;
B_Asum = nansum(B_A);
A_B2sum = nansum(A_B2);

%Correlation coefficient
corrcoeff = nansum(A_av.*B_av)./sqrt(A_av2sum.*B_av2sum);

%Coefficient of determination
R2 = corrcoeff.*corrcoeff;

%Standard deviations of matrices A and B
A_std = sqrt(A_av2sum./n_obs);
B_std = sqrt(B_av2sum./n_obs);

%Root mean square error
rmse = sqrt(A_B2sum./n_obs);

%Nash-Sutcliffe efficiency coefficient
NSE = 1 - A_B2sum./A_av2sum;

%Model bias
bias = B_Asum./n_obs;

%Normalized bias
Nbias = bias./A_std;

%Percent bias
Pbias = 100*B_Asum./nansum(A); 

%Unbiased root mean square difference
rmseu = sqrt(nansum((B_av-A_av).^2)./n_obs);

%Normalised unbiased root mean square difference
rmsenu = sign(B_std-A_std).*rmseu./A_std;

end
