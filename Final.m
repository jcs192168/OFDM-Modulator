%ECE 161C Final Project
%Jean-Claude Sleiman
%PID: A12934981
clear all;
close all;

%Part A

%40 MHz sample rate, 1/40 usec/sample, 160 samples/symbol, 4 usec duration
% 1/40 usec/sample, 128 samples/symbol, 3.2 usec transform segment
% spectral channel spacing 1/3.2 usec = 312.5 kHz
x0=[]; % empty array for OFDM with Guard interval
y0=[]; % empty array for OFDM with cyclic prefices
for k=1:100
    fxx=zeros(1,128);
    dat=(floor(4*rand(1,54))-1.5)/1.5+j*(floor(4*rand(1,54))-1.5)/1.5;
    fxx(65+(-27:27))=[dat(1:27) 0 dat(28:54)];
    xx=8*ifft(fftshift(fxx));
    x0=[x0 zeros(1,32) xx];
    y0=[y0 xx(97:128) xx];
end
figure(1)
subplot(3,1,1)
plot(0:640, real(x0(1:641)),'b','linewidth',2)
hold on
for k=1:160:640
    plot(k:k+31,real(x0(k:k+31)),'r','linewidth',2)
end
hold off
grid on
set(gca,'XTick',[0:32:640])
axis([0 640 -1.2 1.2])
title('Real Part 128 Point OFDM Signal, with 32 Sample Guard Intervals')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,2)
plot(0:640, imag(x0(1:641)),'b','linewidth',2)
hold on
for k=1:160:640
    plot(k:k+31,imag(x0(k:k+31)),'r','linewidth',2)
end
hold off
grid on
set(gca,'XTick',[0:32:640])
axis([0 640 -1.2 1.2])
title('Imaginary Part 128 Point OFDM Signal, with 32 Sample Guard Intervals')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,3)
ww=kaiser(2048,8.3)';
ww=50*ww/sum(ww);
fxxx=zeros(1,2048);
mm=0;
for n=1:512:16000-2048
    fx=fft(x0(n:n+2047).*ww);
    fxxx=fxxx+fx.*conj(fx);
   mm=mm+1;
end
fxxx=fxxx/mm;
plot(-0.5:1/2048:0.5-1/2048,fftshift(20*log10(abs(fx))),'b','linewidth',2)
hold on
plot(-0.5:1/2048:0.5-1/2048,fftshift(10*log10(fxxx)),'r','linewidth',2)
hold off
grid on
axis([-0.5 0.5 -70 5])
title('Single Transform Power Spectrum (Blue), and Average Power Spectrum (Red)')
xlabel('Normalized Frequency')
ylabel('Log Magnitude (dB)')

%Part B
figure(2)
subplot(2,1,1)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x0(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(-64:63,real(fx1),'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis([-32 32 -1.5 1.5])
set(gca,'XTick',[-64:8:64])
title('Real Part of FFT (No Distortion)')

subplot(2,1,2)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x0(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(fx1,'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
title('Constellation (No Distortion)')



%Part C
chan=[1 0 0.2 0 0 0 0 0 0 j*0.1];

x1_dat=filter(chan,1,x0);


figure(3)
subplot(2,1,1)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x1_dat(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(-64:63,real(fx1),'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis([-32 32 -1.5 1.5])
set(gca,'XTick',[-64:8:64])
title('Real Part of FFT [Empty Guard Interval](With Channel Distortion)')

subplot(2,1,2)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x1_dat(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(fx1,'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
title('Constellation [Empty Guard Interval](With Channel Distortion)')

%Part D
chan=[1 0 0.2 0 0 0 0 0 0 j*0.1];
x1_dat=filter(chan,1,y0);


figure(4)
subplot(2,1,1)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x1_dat(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(-64:63,real(fx1),'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis([-32 32 -1.5 1.5])
set(gca,'XTick',[-64:8:64])
title('Real Part of FFT [Cyclic Prefices](With Channel Distortion)')

subplot(2,1,2)
plot(0,0)
hold on
for n=1:160:16000-160
    x1=x1_dat(n:n+159);
    fx1=fftshift(fft(x1(33:160)))/8;
    plot(fx1,'ro','linewidth',2,'markersize',6)
end
hold off
grid on
axis('equal')
axis([-1.5 1.5 -1.5 1.5])
title('Constellation [Cyclic Prefices](With Channel Distortion)')


%Part E


% Zadoff-Chu Sequence for Long Preamble
    
    NN=64;     
    rr=((0:NN-1)-1/2).*((0:NN-1)+0.5)/2;
    prb0=exp(j*2*pi*rr/NN);
    fprb0=fft(fftshift(prb0));
    fprb1=zeros(1,128);
    fprb1(65+(-32:31))=fprb0;
  
    prb1=ifft(fftshift(fprb1));
    prb1a=[prb1(97:128) prb1];
    packet = zeros(1,320);
    packet = [prb1a y0];
    
     % Long Preamble Channel probe
    xprb=filter(chan,1,prb1a);
   	xprb=xprb+0.00*(randn(1,160)+j*randn(1,160))/sqrt(2);
    
    f_input =fftshift(fft(prb1a(33:160)/8));
    f_output=fftshift(fft(xprb(33:160)/8));
    f_eq=zeros(1,128);
    f_eq(65+(-32:31))=f_input(65+(-32:31))./f_output(65+(-32:31));
    f_eq(65+32)=f_eq(65-32);
    f_ch=zeros(1,128);
    f_ch(65+(-32:31))=f_output(65+(-32:31))./f_input(65+(-32:31));
    f_ch(65+23)=f_ch(65-32);
     f_chan=fftshift(abs(fft(chan,128)));
    
    figure(5)
    subplot(2,1,1)
    plot((-0.5:1/128:0.5-1/128)*2,fftshift(abs(fft(xprb(33:160)/8))),'linewidth',2) %estimated channel spectrum
    hold on
    plot((-0.5:1/128:0.5-1/128)*2,f_chan,'r--','linewidth',2) %channel spectrum
    hold off
    grid on
    axis([-1 1 0 1.4])
    title('Spectrum Magnitude; Channel Spectrum, Estimated Channel Spectrum','fontsize',14)
    ylabel('Magnitude','fontsize',14)
    xlabel('Frequency','fontsize',14)
   
    subplot(2,1,2)
    plot((-0.5:1/128:0.5-1/128)*2,abs(f_eq),'linewidth',2)
    hold on
    plot((-0.5:1/128:0.5-1/128)*2,f_chan,'r','linewidth',2)
     plot((-0.5:1/128:0.5-1/128)*2,abs(f_chan.*f_eq),'k','linewidth',2)
     hold off
    grid on
    axis([-1 1 0 1.4])
    title('Spectrum Magnitude; Channel Spectrum, Estimated Channel Equalizer, Equalized Spectrum','fontsize',14)
    ylabel('Magnitude','fontsize',14)
    xlabel('Frequency','fontsize',14)


%Part F

 x5=filter(chan,1,y0)+0.00*(randn(1,16000)+j*randn(1,16000))/sqrt(2);
 fwt=zeros(1,128);
fwt(65+(-27:26))=ones(1,54);
fwt=fftshift(fwt);
  figure(6) 
    subplot(2,1,1)
    plot(0,0)
    hold on
    m=1;
    for k=1:160:length(x5)-160
        x6=x5(k:k+159);
        x6a=x6(33:160);
        fx6a=fftshift(fft(x6a));
        
        fx6a=2*fx6a.*fftshift(fwt);
        x6a=ifft(fftshift(fx6a));
        x6b=x6a(1:2:128);
        
        m=m+1;
        plot(-32:31,real(x6b),'ro','linewidth',2,'markersize',4)
    end
    hold off
    grid on
    axis([-64 63 -1.4 1.4])
    title('Real Part Channel Distorted and Demodulated SC-OFDM Time Signal Without Equalizer','fontsize',14)
    xlabel('Time Index','fontsize',14)
    ylabel('Amplitude','fontsize',14)
    
    subplot(2,1,2)
    plot(0,0)
    hold on
    m=1;
    for k=1:160:length(x5)-160
        x6=x5(k:k+159);
        x6a=x6(33:160);
        fx6a=fftshift(fft(x6a));
        
        fx6a=2*fx6a.*f_eq.*fftshift(fwt);
        x6a=ifft(fftshift(fx6a));
        x6b=x6a(1:2:128);
        
        m=m+1;
        plot(-32:31,real(x6b),'ro','linewidth',2,'markersize',4)
    end
    hold off
    grid on
    axis([-64 63 -1.4 1.4])
    title('Real Part Channel Distorted and Demodulated SC-OFDM Time Signal With Equalizer','fontsize',14)
    xlabel('Time Index','fontsize',14)
    ylabel('Amplitude','fontsize',14)


%Part G


cp00=(floor(2*rand(1,16))-0.5)/0.5+j*(floor(2*rand(1,16))-0.5)/0.5;
cp0=reshape([cp00 cp00 cp00 cp00],1,64); 
fcp0=fftshift(fft(cp0));
    
fcp1=zeros(1,128);
fcp1(65+(-32:31))=fcp0;
fcp1(65-32)=sqrt(0.5)*fcp1(65-32);
fcp1(65+32)=fcp1(65-32);
cp1=2*ifft(fftshift(fcp1));
cp1a=[cp1(97:128) cp1];

 cp2=[zeros(1,50) cp1a zeros(1,10)];
    cp2=filter(chan,1,cp2);
    p2a=cp2.*exp(j*2*pi*(0:219)*1/360);
    p2a=p2a+0.01*(randn(1,220)+j*randn(1,220))/sqrt(2);
    
    reg1=zeros(1,33);
    buf1=zeros(1,32);
    buf2=zeros(1,32);
    cor1=zeros(1,220);
    cor2=zeros(1,220);
    cor3=zeros(1,220);
    sqlsh=zeros(1,220);
    for n=1:220
        reg1=[p2a(n) reg1(1:32)];
        prd1=reg1(1)*conj(reg1(33));
        prd2=reg1(33)*conj(reg1(33));
        
        buf1=[prd1 buf1(1:31)];
        buf2=[prd2 buf2(1:31)];
        
        cor1(n)=sum(buf1)/32;
        cor2(n)=sum(buf2)/32;
        cor3(n)=cor1(n)/(0.01+cor2(n));
        sqlsh(n)=(sign(abs(cor3(n))-0.5)+1)/2;
    end
    
    figure(7)
    subplot(3,1,1)
    plot(0:219,abs(cor1),'linewidth',2)
    hold on
     plot(0:219,abs(cor2),'r','linewidth',2)
     plot([82 82],[-0.05 2.0],'k','linewidth',2)
    hold off
    grid on
    axis([0 220 -0.1 2.5])
    title('Magnitude: Input Short Preamble','fontsize',14)
    ylabel('Magnitude','fontsize',14)
    
    subplot(3,1,2)
    plot(0:219,abs(cor3),'.-','linewidth',2)
    hold on
    plot([40 60]+32,[0.5 0.5],'r','linewidth',2)
      plot(0:219,sqlsh,'--r','linewidth',2)
      hold off
    grid on
    axis([0 220 0 1.2])
    title('Magnitude: Normalized Cross Correlation (Cross/Auto)','fontsize',14)
    ylabel('Magnitude','fontsize',14)
    
    subplot(3,1,3)
    plot(0:219,sqlsh.*angle(cor3)*360/(2*pi*32),'linewidth',2)
    hold on
    plot([82 220],[1 1]*1,'r','linewidth',2)    
    hold off
    axis([0 220 0 1.2])
    grid on
    title('Estimate of Frequency Offset: Cross Correlation Angle Divided by 32','fontsize',14)
    ylabel('Magnitude','fontsize',14)
    xlabel('Time Index','fontsize',14)

