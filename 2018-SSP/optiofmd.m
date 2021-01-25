   clc,clear;
tic
% =========== Initialization ===============
n=2;
N=10; % Number of subcarriers
P=1; % the primary user's tx power
Qmax = P*N; % The maximum total power
sigma2=1; % Noise Power.
Ave_rate=0;
 RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
for i=1:100
    f=(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));  % Secondary Channel Gain
    ff(:,i)=f(:);
%     f(:)=ff(:,i);
    for cnt=1:N
        F(:,:,cnt)=[f(cnt),0;0,conj(f(cnt))]; % Secondary Channel Matrix.
    end
    
    g=(1/sqrt(2))*sqrt(1)*(randn(1,N)+1i*randn(1,N));  % Secondary Channel Gain
    gg(:,i)=g(:);
%     g(:)=gg(:,i);
    for cnt=1:N
        G(:,:,cnt)=[g(cnt),0;0,conj(g(cnt))]; % Secondary Channel Matrix.
    end



    h=(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));  % Secondary Channel Gain
    hh(:,i)=h(:);
%     h(:)=hh(:,i);
    for cnt=1:N
        H(:,:,cnt)=[h(cnt),0;0,conj(h(cnt))]; % Secondary Channel Matrix.
    end
    
    
    h2=(1/sqrt(2))*(randn(1,N)+1i*randn(1,N));  % Secondary Channel Gain
    hh2(:,i)=h2(:);
%     h2(:)=hh2(:,i);
    for cnt=1:N
        H2(:,:,cnt)=[h2(cnt),0;0,conj(h2(cnt))]; % Secondary Channel Matrix.
    end
    
    
    Gain_S(:,i)=(1+abs(h2).^2)./(abs(f).^2);
    
    
    Rpi=0;
    for cnt=1:N
        Rpi=.5*log2(det(eye(n)+(sigma2*eye(n))^(-1)*H(:,:,cnt)*P*H(:,:,cnt)'))+Rpi;
    end

    Rth=0.7*Rpi; % the minimum acceptable rate for primary user

    %Rpi=1;
    Q0=.00*ones(n,n,N);Rpi=0;
    for cnt=1:N
        Rpi=.5*log2(det(eye(n)+(sigma2*eye(n)+G(:,:,cnt)*Q0(:,:,cnt)*G(:,:,cnt)')^(-1)*H(:,:,cnt)*P*H(:,:,cnt)'))+Rpi;
    end

    % =========================================
    % =========================================
    cnstrt = 1;count=0;cvx_clear;
    while cnstrt == 1 && count<40
        cvx_begin quiet
            variable Q(n,n,N) hermitian semidefinite %toeplitz
            Obj=0;
            for cnt=1:N
                Obj = Obj+sum(0.5*log_det( eye(n) + (sigma2*eye(2)+H2(:,:,cnt)*P*H2(:,:,cnt)')^(-1)*F(:,:,cnt)*Q(:,:,cnt)*F(:,:,cnt)'));
            end
            A=0;B=0;

            for cnt=1:N
                DummyVar = -.5*1/log(2)*G(:,:,cnt)'*(sigma2*eye(n)+G(:,:,cnt)*Q0(:,:,cnt)*G(:,:,cnt)')^(-1)*H(:,:,cnt)*P*H(:,:,cnt)'* ... 
                    (eye(n)+(sigma2*eye(n)+G(:,:,cnt)*Q0(:,:,cnt)*G(:,:,cnt)')^(-1)*H(:,:,cnt)*P*H(:,:,cnt)')^(-1)... 
                    *(sigma2*eye(n)+G(:,:,cnt)*Q0(:,:,cnt)*G(:,:,cnt)')^(-1)*G(:,:,cnt);

                A = trace(  DummyVar*(Q(:,:,cnt)-Q0(:,:,cnt))  )+A;
            end
            for cnt=1:N
                B = B+sum(diag(Q(:,:,cnt)));
            end
            maximize( Obj )
            subject to
                Q(1,1,:) == Q(2,2,:);
                Q(1,2,:) == 0; % To make it proper
                B <= 2*Qmax;
                real(A - Rth+Rpi)>= 0;
                %Z == hermitian_semidefinite( n ) : Q;
        cvx_end
        temp=0;
        for cnt=1:N
            temp = temp+trace(abs(Q0(:,:,cnt)-Q(:,:,cnt)));
        end
        if temp < 0.0001
            cnstrt=0;
        else
            Q0=Q;Rpi=0;
            for cnt=1:N
                Rpi=.5*log2(det(eye(n)+(sigma2*eye(n)+G(:,:,cnt)*Q0(:,:,cnt)*G(:,:,cnt)')^(-1)*H(:,:,cnt)*P*H(:,:,cnt)'))+Rpi;
            end
        end
        count=count+1;

    end
    Power_Secondary(:,:,:,i)=Q;
    Power_Secondary_2(:,i)=Q(1,1,:);
    Improp(:,i) = real(Power_Secondary(2,1,:,i)./Power_Secondary(1,1,:,i));
    Rate_Secondary(i)=Obj;
%     Rate_Secondary_imp(i)=Obj;
    Ave_Rate=Ave_rate+Obj/100;
    i
end
toc