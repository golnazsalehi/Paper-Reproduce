SEP_10 = Both_Methods(10,10,10000);
SEP_4 = Both_Methods(4,4,10000);


semilogy(snr,SEP_4(1,:),'r--x','linewidth',1.5)
hold on
semilogy(snr,SEP_4(2,:),'r-d','linewidth',1.5)
semilogy(snr,SEP_10(1,:),'b-+','linewidth',1.5)
semilogy(snr,SEP_10(2,:),'b-*','linewidth',1.5)
xlabel('rho (dB)');ylabel('SER')
legend('4x4 Channal Inversion','4x4 Regularized Channal Inversion',...
    '10x10 Channal Inversion','10x10 Regularized Channal Inversion')

function  SEP = Both_Methods(K,M,N)
    interval = -10:5:30;
    
    SNR = 10.^((interval)/10);
    N0_vector = 1./SNR;
    SEP = zeros(2,length(SNR));
    for i = 1:length(SNR)
    
       
        for j = 1:N
            x = randsrc(K,1,0:3);
            u = qammod(x,4);
            H = (randn(K,M)+1i*randn(K,M))*sqrt(0.5); % channel
    
    
            alpha = K*N0_vector(i);
            w = sqrt(N0_vector(i)/2)*(randn(1,K)+1i*randn(1,K)); % noise
            w = w.';
    
            s_ChInv = H\u;
            s_ChInv = s_ChInv/sqrt(s_ChInv'*s_ChInv);
    
            s_ReInv =H'/(H*H'+alpha*eye(K))*u;
            s_ReInv = s_ReInv/sqrt(s_ReInv'*s_ReInv);
    
            y_ChInv = H*s_ChInv + w;
            y_ReInv = H*s_ReInv + w;
            
    
            u_ChInv = qamdemod(y_ChInv,4);
            
            u_ReInv = qamdemod(y_ReInv,4);
    
            SEP(1,i) = SEP(1,i) + sum(x ~= u_ChInv);
            SEP(2,i) = SEP(2,i) + sum(x ~= u_ReInv);
        end
    end
    SEP = SEP./(N*M);
end