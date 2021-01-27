function [r_uni, r, R_1, R_2, p_new_1, p_new_2, t_real] = r_MIMO_1(H, p_new_1, p_new_2)



global P1 P2 alpha N_t N_r sigma sigma_d sigma_r coun coun2 %K %sigma_eta sigma_eta_c 

% -----------------------------------------------
% --- Extracting the corresponding channels & ---
% --------- Making equivalent channels ----------
% ----------------------------------------------- 
H_11 = H(1:N_r,1:N_t);

H_12 = H(1:N_r,(N_t+1):2*N_t);

H_21 = H((N_r+1):2*N_r,1:N_t);

H_22 = H((N_r+1):2*N_r,(N_t+1):2*N_t);



tic

if alpha==0 
    R_1=0;
    % =========================================
    C_z2 = sigma*eye(N_r,N_r)+sigma_d*(H_22*diag(diag(p_new_2))*H_22');
    % --------------------------------------
    a_2 =  real(log(det(C_z2)));
    % --------------------------------------
    c_2 =sigma_d*diag(diag(H_22'*((C_z2)^(-1))*H_22));
    % ===========================================
    % ===========================================
    R_2 = real(log(det(C_z2 + H_22*p_new_2*H_22')) - a_2);
    % -----------------------------------------
    r_new=R_2;
    r_uni=r_new;
    cnstrt2=1;iter_n=0;
    while cnstrt2==1 && iter_n<30
        iter_n=iter_n+1;
        r_pre=r_new;
        p_pre_2=p_new_2;
        % ------------------------------------------
        % ------------------------------------------
        cvx_begin sdp quiet
           variable p_2(N_t,N_t) hermitian semidefinite %nonnegative
           variable u nonnegative

           maximize u
           subject to  
            real(trace(p_2))<=P2;
            
            real(log_det(sigma*eye(N_r,N_r)+sigma_d*(H_22*diag(diag(p_2))*H_22') ...
                + H_22*p_2*H_22')) - a_2  ...
                 - real(trace(c_2'*(p_2-p_pre_2)))>=u;%

        cvx_end
        % ------------------------------------------
        p_new_2 = p_2; %p_new_2 = zeros(N_t, N_t);
        % --------------------------------------
        % ===========================================
        % ===========================================
          
         % --------------------------------------
        C_z2 = sigma*eye(N_r,N_r)+sigma_d*(H_22*diag(diag(p_new_2))*H_22');
        % --------------------------------------
        a_2 =  real(log(det(C_z2)));
        % --------------------------------------
        c_2 =sigma_d*diag(diag(H_22'*((C_z2)^(-1))*H_22));
        % ===========================================
        % ===========================================
        R_2 = real(log(det(C_z2 + H_22*p_new_2*H_22')) - a_2);
        % -----------------------------------------
        r_new=R_2;
        
        if abs(r_new-r_pre)/r_new<10^(-4)
            cnstrt2=0;
        elseif r_new-r_pre<0
            abs(r_new-r_pre)/r_new
            if abs(r_new-r_pre)/r_new>0.01
                coun2=coun2+1;
            end
            cnstrt2=0;keyboard;
            coun=coun+1;
        end
    end
    R_2=R_2/log(2);

elseif alpha==1
    R_2=0;
    % =========================================
    p_new_2 = zeros(N_t, N_t);
    % --------------------------------------
    C_z1 = sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_new_1))*H_11');
    % --------------------------------------
    a_1 =  real(log(det(C_z1)));
    % --------------------------------------
    c_1 =sigma_d*diag(diag(H_11'*((C_z1)^(-1))*H_11));
    % ===========================================
    % ===========================================
    R_1 = real(log(det(C_z1 + H_11*p_new_1*H_11')) - a_1)/log(2);
    % -----------------------------------------
    r_new=R_1;
    r_uni=r_new;
    cnstrt2=1;iter_n=0;
    while cnstrt2==1 && iter_n<30
        iter_n=iter_n+1;
        p_pre_1 = p_new_1; 
        r_1_pre = R_1;
        r_pre=r_new;
        % ------------------------------------------
        a_1_pre = a_1;
        % ------------------------------------------
        % ------------------------------------------
        cvx_begin sdp quiet
           variable p_1(N_t,N_t) hermitian semidefinite %nonnegative
           variable u nonnegative

           maximize u
           subject to  
            real(trace(p_1))<=P1;
            
            real(log_det(sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_1))*H_11') ...
                + H_11*p_1*H_11')) - a_1  ...
                 - real(trace(c_1'*(p_1-p_pre_1)))>=u;%

        cvx_end
        % ------------------------------------------
        p_new_1 = p_1; p_new_2 = zeros(N_t, N_t);
        % --------------------------------------
        % ===========================================
        % ===========================================
          
         % --------------------------------------
        C_z1 = sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_new_1))*H_11');
        % --------------------------------------
        a_1 =  real(log(det(C_z1)));
        % --------------------------------------
        c_1 =sigma_d*diag(diag(H_11'*((C_z1)^(-1))*H_11));
        % ===========================================
        % ===========================================
        R_1 = real(log(det(C_z1 + H_11*p_new_1*H_11')) - a_1)/log(2);
        % -----------------------------------------
        r_new = R_1;
        
        if abs(r_new-r_pre)/r_new<10^(-4)
            cnstrt2=0;
        elseif r_new-r_pre<0
            abs(r_new-r_pre)/r_new
            if abs(r_new-r_pre)/r_new>0.01
                coun2=coun2+1;
            end
            cnstrt2=0;keyboard;
            coun=coun+1;
        end
    end
    
    

else
    % =========================================
    % =========================================
    C_z1 = sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_new_1))*H_11' ...
        + H_12*diag(diag(p_new_2))*H_12')+sigma_r*diag(diag(H_11*p_new_1*H_11' ...
        + H_12*p_new_2*H_12'));
    C_z2 = sigma*eye(N_r,N_r)+sigma_d*(H_21*diag(diag(p_new_1))*H_21' ...
        + H_22*diag(diag(p_new_2))*H_22')+sigma_r*diag(diag(H_22*p_new_2*H_22' ...
        + H_21*p_new_1*H_21'));
    
    % --------------------------------------
    a_1 =  real(log(det(C_z1 + H_12*p_new_2*H_12')));
    a_2 =  real(log(det(C_z2 + H_21*p_new_1*H_21')));
    % --------------------------------------
    temp=zeros(N_t,N_t);
    if sigma_r>0
        for i=1:N_t
            for j=1:N_t
                temp2=zeros(N_t,N_t);
                temp2(i,j)=1;
                temp(i,j)=trace(((C_z1 + H_12*p_new_2*H_12')^(-1)) ... 
                    *diag(diag(H_12*temp2*H_12')));
            end
        end
    end
    
    b_1 =  H_12'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_12 ...
        +(sigma_d)*diag(diag(H_12'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_12)) ...
        +sigma_r*temp;
    
    temp=zeros(N_t,N_t);
    if sigma_r>0
        for i=1:N_t
            for j=1:N_t
                temp2=zeros(N_t,N_t);
                temp2(i,j)=1;
                temp(i,j)=trace(((C_z2 + H_21*p_new_1*H_21')^(-1)) ... 
                    *diag(diag(H_21*temp2*H_21')));
            end
        end
    end
    b_2 =  H_21'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_21 ...
        +(sigma_d)*diag(diag(H_21'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_21)) ...
        +sigma_r*temp;
    % --------------------------------------
    
    temp=zeros(N_t,N_t);
    if sigma_r>0
        for i=1:N_t
            for j=1:N_t
                temp2=zeros(N_t,N_t);
                temp2(i,j)=1;
                temp(i,j)=trace(((C_z1 + H_12*p_new_2*H_12')^(-1)) ... 
                    *diag(diag(H_11*temp2*H_11')));
            end
        end
    end
    c_1 =  (sigma_d)*diag(diag(H_11'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_11)) ...
        +sigma_r*temp;
    
    temp=zeros(N_t,N_t);
    if sigma_r>0
        for i=1:N_t
            for j=1:N_t
                temp2=zeros(N_t,N_t);
                temp2(i,j)=1;
                temp(i,j)=trace(((C_z2 + H_21*p_new_1*H_21')^(-1)) ... 
                    *diag(diag(H_22*temp2*H_22')));
            end
        end
    end
    c_2 =  (sigma_d)*diag(diag(H_22'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_22)) ... 
        +sigma_r*temp;
    % ===========================================
    % ===========================================
    R_1 = real(log(det(C_z1 + H_11*p_new_1*H_11' ...
        + H_12*p_new_2*H_12')) - a_1);
    R_2 = real(log(det(C_z2 + H_22*p_new_2*H_22' ...
        + H_21*p_new_1*H_21')) - a_2);
    % -----------------------------------------
    r_new=min(R_1/alpha,R_2/(1-alpha));
    r_uni=r_new;
    cnstrt2=1;iter_n=0;
    while cnstrt2==1 && iter_n<30
        iter_n=iter_n+1;
        p_pre_1 = p_new_1; 
        p_pre_2 = p_new_2; 
        r_1_pre = R_1;
        r_2_pre = R_2;
        r_pre=r_new;
        % ------------------------------------------
        a_1_pre = a_1;
        a_2_pre = a_2;
        cvx_begin sdp quiet
           variable p_1(N_t,N_t) hermitian semidefinite %nonnegative
           variable p_2(N_t,N_t) hermitian semidefinite %nonnegative
           variable u nonnegative

           maximize u
           subject to  
            real(trace(p_1))<=P1;
            real(trace(p_2))<=P2;
            
            real(log_det(sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_1))*H_11' ...
                + H_12*diag(diag(p_2))*H_12') + H_11*p_1*H_11' ...
                + H_12*p_2*H_12' + sigma_r*diag(diag(H_11*p_1*H_11' ...
                + H_12*p_2*H_12')))) - a_1 - real(trace(b_1'*(p_2-p_pre_2)))  ...
                 - real(trace(c_1'*(p_1-p_pre_1)))>=alpha*u;%
             
             real(log_det(sigma*eye(N_r,N_r)+sigma_d*(H_21*diag(diag(p_1))*H_21' ...
                + H_22*diag(diag(p_2))*H_22') + H_22*p_2*H_22' ...
                + H_21*p_1*H_21' + sigma_r*diag(diag( H_22*p_2*H_22' ...
                + H_21*p_1*H_21')))) - a_2 - real(trace(b_2'*(p_1-p_pre_1)))  ...
                 - real(trace(c_2'*(p_2-p_pre_2)))>=(1-alpha)*u;%>=(1-alpha)*u;

        cvx_end
        % ------------------------------------------
        p_new_1 = p_1; p_new_2 = p_2;
        % --------------------------------------
        % ===========================================
        % ===========================================
          
        % --------------------------------------
        C_z1 = sigma*eye(N_r,N_r)+sigma_d*(H_11*diag(diag(p_new_1))*H_11' ...
            + H_12*diag(diag(p_new_2))*H_12')+sigma_r*diag(diag(H_11*p_new_1*H_11' ...
            + H_12*p_new_2*H_12'));
        C_z2 = sigma*eye(N_r,N_r)+sigma_d*(H_21*diag(diag(p_new_1))*H_21' ...
            + H_22*diag(diag(p_new_2))*H_22')+sigma_r*diag(diag(H_22*p_new_2*H_22' ...
            + H_21*p_new_1*H_21'));

        % --------------------------------------
        a_1 =  real(log(det(C_z1 + H_12*p_new_2*H_12')));
        a_2 =  real(log(det(C_z2 + H_21*p_new_1*H_21')));
        % --------------------------------------
        temp=zeros(N_t,N_t);
        if sigma_r>0
            for i=1:N_t
                for j=1:N_t
                    temp2=zeros(N_t,N_t);
                    temp2(i,j)=1;
                    temp(i,j)=trace(((C_z1 + H_12*p_new_2*H_12')^(-1)) ... 
                        *diag(diag(H_12*temp2*H_12')));
                end
            end
        end

        b_1 =  H_12'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_12 ...
            +(sigma_d)*diag(diag(H_12'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_12)) ...
            +sigma_r*temp;

        temp=zeros(N_t,N_t);
        if sigma_r>0
            for i=1:N_t
                for j=1:N_t
                    temp2=zeros(N_t,N_t);
                    temp2(i,j)=1;
                    temp(i,j)=trace(((C_z2 + H_21*p_new_1*H_21')^(-1)) ... 
                        *diag(diag(H_21*temp2*H_21')));
                end
            end
        end
        b_2 =  H_21'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_21 ...
            +(sigma_d)*diag(diag(H_21'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_21)) ...
            +sigma_r*temp;
        % --------------------------------------

        temp=zeros(N_t,N_t);
        if sigma_r>0
            for i=1:N_t
                for j=1:N_t
                    temp2=zeros(N_t,N_t);
                    temp2(i,j)=1;
                    temp(i,j)=trace(((C_z1 + H_12*p_new_2*H_12')^(-1)) ... 
                        *diag(diag(H_11*temp2*H_11')));
                end
            end
        end
        c_1 =  (sigma_d)*diag(diag(H_11'*((C_z1 + H_12*p_new_2*H_12')^(-1))*H_11)) ...
            +sigma_r*temp;

        temp=zeros(N_t,N_t);
        if sigma_r>0
            for i=1:N_t
                for j=1:N_t
                    temp2=zeros(N_t,N_t);
                    temp2(i,j)=1;
                    temp(i,j)=trace(((C_z2 + H_21*p_new_1*H_21')^(-1)) ... 
                        *diag(diag(H_22*temp2*H_22')));
                end
            end
        end
        c_2 =  (sigma_d)*diag(diag(H_22'*((C_z2 + H_21*p_new_1*H_21')^(-1))*H_22)) ... 
            +sigma_r*temp;
        % ===========================================
        % ===========================================
        R_1 = real(log(det(C_z1 + H_11*p_new_1*H_11' ...
            + H_12*p_new_2*H_12')) - a_1);
        R_2 = real(log(det(C_z2 + H_22*p_new_2*H_22' ...
            + H_21*p_new_1*H_21')) - a_2);
        % -----------------------------------------
        % -----------------------------------------
        r_new=min(R_1/alpha,R_2/(1-alpha));
        
        if abs(r_1_pre-R_1)/R_1<10^(-4)&&abs(r_2_pre-R_2)/R_2<10^(-4) ...
                &&abs(r_new-r_pre)/r_new<10^(-4)
            cnstrt2=0;
        elseif r_new-r_pre<0
            abs(r_new-r_pre)/r_new
            if abs(r_new-r_pre)/r_new>0.01
                coun2=coun2+1;
            end
            cnstrt2=0;keyboard;
            coun=coun+1;
        end
%         if r_new-r_pre<0
%             keyboard;
%         end
    end
    
    R_1=R_1/log(2);
    R_2=R_2/log(2);
    
end

%r=min(R_1/alpha,R_2/alpha);

r=min(R_1,R_2);

t_real=toc;

