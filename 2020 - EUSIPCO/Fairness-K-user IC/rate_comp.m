function [r] ...
    = rate_comp(H, p_new)

global  N_t N_r sigma sigma_d  K 


for k=1:K
    for i=1:K
         eval(sprintf('H_%i_%i(:,:)= H((1+N_r*(k-1)):k*N_r,(1+N_t*(i-1)):i*N_t);',k,i)) 
         
    end
end



for k=1:K
    eval(sprintf('p_new_%i(:,:)=p_new(:,:,k);',k))
end


for k=1:K
clear temp
temp=zeros(N_r,N_r);
for i=1:K
    temp = temp ...
        + eval(sprintf('H_%i_%i',k,i)) ...
        *diag(diag(eval(sprintf('p_new_%i',i))))* ...
        eval(sprintf('H_%i_%i',k,i))';
end
eval(sprintf('C_z%i = sigma*eye(N_r,N_r)+sigma_d*temp;',k))

temp=zeros(N_r,N_r);
for i=1:K
    if i~=k
        temp = temp ...
            + eval(sprintf('H_%i_%i',k,i)) ...
            *eval(sprintf('p_new_%i',i))* ...
            eval(sprintf('H_%i_%i',k,i))';
    end
end

eval(sprintf('inr_%i = C_z%i + temp;',k,k))

a(k) =  real(log(det( ...
    eval(sprintf('inr_%i',k)))));

R(k) = real(log(det( ...
    eval(sprintf('inr_%i',k)) ...
    + eval(sprintf('H_%i_%i',k,k))* ...
    eval(sprintf('p_new_%i',k))* ...
    eval(sprintf('H_%i_%i',k,k))' ...
    )) - a(k))/log(2);


cc =  sigma_d*diag(diag( ...
    eval(sprintf('H_%i_%i',k,k))'* ...
    ((eval(sprintf('inr_%i',k)))^(-1))* ...
    eval(sprintf('H_%i_%i',k,k)) ...
    ));

eval(sprintf('c_%i=cc;',k))



for i=1:K
    if i~=k
        clear temp1 temp2
        temp1 = eval(sprintf('H_%i_%i',k,i))'* ...
            ((eval(sprintf('inr_%i',k)))^(-1))* ...
            eval(sprintf('H_%i_%i',k,i));

        temp2 = diag(diag( ...
            eval(sprintf('H_%i_%i',k,i))'* ...
            ((eval(sprintf('inr_%i',k)))^(-1))* ...
            eval(sprintf('H_%i_%i',k,i)) ...
            ));
        %--------------------------------------------
        eval(sprintf('b_%i_%i(:,:) = temp1 + sigma_d*temp2;',k,i));
    end
 end
end


r=min(R);