
function Fi = aumFeat(A,D)



S = size(A);
N=S(2);
F=[];
for k = 1:N
[~,idx]=sort(abs(A(find(A(:,k)),k)));
L=length(idx);
x=(A(find(A(:,k)),k).*D(:,find(A(:,k)))')';

for j = 1:L
i = idx(j);

X = full(x);
Sf = X(:,i);

th_ssc=0.001;
th_zc=0.11;
th_wamp = 0.06;

Ff(1)=ZC(th_zc,Sf);
Ff(2)=WAMP(th_wamp,Sf);
Ff(3)=RMS(Sf);
Ff(4)=VAR(Sf);
Ff(5)=DASDV(Sf);
Ff(6)=IAV(Sf);
Ff(7)=MFL(Sf);
Ff(8)=MSR(Sf);
Ff(9)=LS(Sf,2);

F = [F;Ff'];
end
Fi(:,k) = F;
F = [];
end

end


 function F=pooling(A)
F = [];
 K = [1 2 4 5 10];
 for j = 1:5
 k = K(j);
% A = full(As);
% k=4;
S = size(A);
N = S(1);
%P=permn([1 0],k);
%P = P(1:end-1,:);
P = eye(k);
n = N/k;

for i = 1:k
    ind = [1+n*(i-1):n*i];
    f1(i,:) = sum(abs(A(ind,:)));
    f2(i,:) = max(abs(A(ind,:)));
end

 for q = 1:length(P)
     F1(q,:) = sum(f1(find(P(q,:)),:),1);
     F2(q,:) = max(f2(find(P(q,:)),:),[],1);
 end

F = [[F1' F2'] F];

%F = [f1' f2'];
end
end