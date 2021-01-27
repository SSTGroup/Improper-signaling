function [r_uni, r, R, t_real, p_new] = r_MIMO_k(H)



global P  N_t N_r sigma sigma_d coun coun2 K 


a  = zeros(K,1);
R  = zeros(K,1);
p_new = zeros(N_t,N_t,K);

% -----------------------------------------------
% --- Extracting the corresponding channels & ---
% --------- Making equivalent channels ----------
% ----------------------------------------------- 

for k=1:K
    for i=1:K
         eval(sprintf('H_%i_%i(:,:)= H((1+N_r*(k-1)):k*N_r,(1+N_t*(i-1)):i*N_t);',k,i)) 
         
    end
end

tic

for k=1:K
    eval(sprintf('p_new_%i=(P/(N_t))*eye(N_t,N_t);',k,k))
end

% =========================================
% =========================================
% --------------------------------------
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
% -----------------------------------------
r_new=min(R);
r_uni=r_new;
cnstrt2=1;iter_n=0;
while cnstrt2==1 && iter_n<30
    iter_n=iter_n+1;
    for k=1:K
        eval(sprintf('p_pre_%i=p_new_%i;',k,k))
    end
    r_pre = r_new;
    % ------------------------------------------
    clear tp
    % --------------------------
    cvx_begin sdp quiet
       for k=1:K
          eval(sprintf('variable p_%i(N_t,N_t) hermitian semidefinite ',k));
       end
       variable u nonnegative
       
       for k=1:K
           clear temp
            temp=zeros(N_r,N_r);
            for i=1:K
                temp = temp ...
                    + eval(sprintf('H_%i_%i',k,i)) ...
                    *diag(diag(eval(sprintf('p_%i',i))))* ...
                    eval(sprintf('H_%i_%i',k,i))';
            end
            eval(sprintf('C_z_%i = sigma*eye(N_r,N_r)+sigma_d*temp;',k))

            temp=zeros(N_r,N_r);
            for i=1:K
                if i~=k
                    temp = temp ...
                        + eval(sprintf('H_%i_%i',k,i)) ...
                        *eval(sprintf('p_%i',i))* ...
                        eval(sprintf('H_%i_%i',k,i))';
                end
            end

            eval(sprintf('inr_cvx_%i = C_z_%i + temp;',k,k))
       end
       
       tp(1) = real(trace(c_1'*(p_1-p_pre_1))) ...
           + real(trace(b_1_2'*(p_2- p_pre_2)));
       
       if K>2
           for i=3:K
               tp(1) = real(trace(eval(sprintf('b_1_%i',i))'* ...
                (eval(sprintf('p_%i',i))- ... 
                eval(sprintf('p_pre_%i',i))))) ...
                   +tp(1);

           end
       end
       
       for k=2:K
           tp(k)=0;
           for i=1:K
               if i~=k
                   tp(k)=real(trace(eval(sprintf('b_%i_%i',k,i))'* ...
                       (eval(sprintf('p_%i',i))- ... 
                       eval(sprintf('p_pre_%i',i)))))... 
                       + tp(k);
               else
                  tp(k)=real(trace(eval(sprintf('c_%i',k))'* ...
                     (eval(sprintf('p_%i',k))- ... 
                      eval(sprintf('p_pre_%i',k)))))... 
                     + tp(k);
               end
           end
       end


       maximize u
       subject to 
        for k=1:K
            % -----------------------------
            % ---- Power Constraint --------
            % -----------------------------
            real(trace(eval(sprintf('p_%i',k))))<=P;
            
            
            % -----------------------------
            % ---- Rate Constraint --------
            % -----------------------------
            real(log_det( ...
                eval(sprintf('inr_cvx_%i',k)) ...
                + eval(sprintf('H_%i_%i',k,k))* ...
                eval(sprintf('p_%i',k))* ...
                eval(sprintf('H_%i_%i',k,k))' ...
                )) - a(k) - tp(k) >=u;%

        end


    cvx_end
    
    % ------------------------------------------
    for k=1:K
        eval(sprintf('p_new_%i=p_%i;',k,k))
        p_new(:,:,k)=eval(sprintf('p_%i',k));
    end
    % --------------------------------------
    
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
    % -----------------------------------------
    % -----------------------------------------
    r_new=min(R);
    
    if abs(r_new-r_pre)/r_new<10^(-4)
        cnstrt2=0;
    elseif r_new-r_pre<0
        abs(r_new-r_pre)/r_new
        if abs(r_new-r_pre)/r_new>0.01
            coun2=coun2+1;
        end
        cnstrt2=0;%keyboard;
        coun=coun+1;
    end

end



r=min(R);

t_real=toc;

