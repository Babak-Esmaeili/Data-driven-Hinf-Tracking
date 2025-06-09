function noise_scale = snr_scale(signal,noise,snr_dB)

    noise_scale = (norm(signal)/norm(noise))*(1/(10.0^(0.05*snr_dB)));

end