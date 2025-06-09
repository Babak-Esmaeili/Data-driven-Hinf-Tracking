function [noisy_signal,scaled_noise,noise] = add_noise_func(signal,alpha,snr_dB)
    
    N = length(signal);
    noise = randn(N,1);
    
    filtered_noise = filter([1 -alpha],[1 -1],noise);
    
    scale = (norm(signal)/norm(filtered_noise)) / 10.0^(0.05*snr_dB);
    
    scaled_noise = scale.*filtered_noise;
    
    noisy_signal = signal + scaled_noise;
    
    % Validate scalization
    snr_dB_new = 20*log10(norm(signal)/norm(scaled_noise));
    assert(abs(snr_dB_new - snr_dB) < 1e10*eps);

end