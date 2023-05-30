function [out,f,t, fig]=detector(signal,fs,win,thr)

% function [out,f,t]=detector(signal,fs,win,thr)
% 
% signal - measured sensor signal
% fs - measurement frequency [3000]
% win - window size (in seconds) [0.5]
% thr - detection threshold [2] 

signal = (signal - 2048)/2048;

if nargin<4
    thr=2;
end
if nargin<3
    win=0.5;
end
if nargin<2
    fs=3000;
end


% disp(length(signal)/fs)
 fig = figure;
 sgtitle('With Matlab fft implementation')
 
 subplot(3,1,1)
 
 plot([1:length(signal)]/fs,signal)
 diap=max(signal)-min(signal);
 ylabel('voltage [V]')
 axis([win length(signal)/fs min(signal)-0.1*diap max(signal)+0.1*diap]);

 f0=20; 
 [f,t,spctr,acc_en,locs]=spgram(signal,fs,win,f0,thr);
 
 spctr=10*log10(spctr);
 
 subplot(3,1,2)
 imagesc(t,f/17.3148,spctr');
 hold on
 set(gca,'YDir','normal');
 
 for m=1:size(locs,1)
   % plot(t([locs(m,1) locs(m,1)]),[min(acc_en) max(acc_en)],'g')
   % plot(t([locs(m,2) locs(m,2)]),[min(acc_en) max(acc_en)],'g')
   if ~locs(m,4) % first item of the group
       % do nothing
   else
   if ~locs(m,3) % second (we assume last item of the group here
       locs(m,3)=locs(m-1,3);
       locs(m,5)=sum(locs(m-1:m,5));
   end    
       
    if locs(m,5)>locs(m,6)  % car from the left
%     plot(t([locs(m,3) locs(m,3)]),[fl maxfr]/spf,'r')
%     plot(t([locs(m,3)+3 locs(m,3)+3]),[fl maxfr]/spf,'r')
     plot(t([locs(m,3) locs(m,3)+3]),f([locs(m,7) locs(m,7)])/17.3148,'w','LineWidth',4);
    else
%     plot(t([locs(m,4) locs(m,4)]),[fl maxfr]/spf,'r')
%     plot(t([locs(m,4)-3 locs(m,4)-3]),[fl maxfr]/spf,'r')
     plot(t([locs(m,4)-3 locs(m,4)]),f([max(1,locs(m,7)) max(1,locs(m,7))])/17.3148,'w','LineWidth',4);
    end   
   end     
    
end
 
 
 axis([t(1) t(end) f(1)/17.3148 f(end)/17.3148]);
 ylabel('speed [km/h]')
 xlabel('time [s]')
 %colormap('hot')
 
 subplot(3,1,3)
 plot(t,acc_en);
 hold on
 axis([min(t) max(t) 0 1.1*max(acc_en)]);
 ylabel('total spectral power [W/Hz]');
 xlabel('time [s]');
 for l=1:size(locs,1),
    plot(t([locs(l,1) locs(l,1)]),[0 1*max(acc_en)],'r')
    plot(t([locs(l,2) locs(l,2)]),[0 1*max(acc_en)],'r')    
    % speed/direction markers
    
 end
plot([t(1) max(t)],[thr thr],'y');
%out=10*log10(spctr);
out=locs;

end


function [f,t,spctr,acc_en,locs]=spgram(y,fs,win,f0,thr)

% given a signal y of length l, sampling frequency fs, window size (in seconds)
% and frequency cutoff f0 (Hz)
% this function calculates the spectrogram (matrix) spctr, 
% with frequency and time axes in f and t, respectively
    locs=[];
    win=win*fs;
    overlap=win-fs/10;
    % window moved by 0.2*win

    l=length(y);

    dw=2;
    w=5;
    wt=0;
    
    spctr=[];
    alg=1;
    lopp=win;
    i=1;
    k=0;
    
    disp('detection results')
    disp('****************************************************')
    disp(' car no. |    t_a |    t_b | direction |     speed ')
    disp('****************************************************')
    
    dstr=[];
    
    while lopp<=l    
        t(i)=lopp/fs;
        [f,F]=spekter(y,fs,alg,lopp);
        if nargin<4
        % do nothing
        else
            F=F(f>f0);
            f=f(f>f0);
        end
        spctr=[spctr; F'];  % actual spectrum

        % total spectral energy
        acc_en(i)=sum(F);
                
        % actual detection
        if i>1
         if acc_en(i-1)< thr & acc_en(i)>=thr    
           disp(dstr)  
           k=k+1;  % a new row           
           dstr=[];
           dstr= [padd(sprintf('%d',k),8,'left') ' |' padd(sprintf('%.1fs',t(i-1)),7,'left') ' |'];
           locs(k,1)=i-1;  % previous step             
         if k>1 & ~locs(k-1,4)  % right side energy of the previous car not checked yet
             % meaning that it is very close to the current one
         else
           locs(k,3)=i-1-w;
           locs(k,5)=sum(acc_en(max(i-1-w,1):max(i-1-dw,1)));
         end
         end
        
        
         if acc_en(i-1)>=thr & acc_en(i)<thr  

           if k==0  % this means that the stream starts with a cut signal
              k=k+1;  % a new row           
              locs(k,1)=1;           
              dstr= [padd(sprintf('%d',k),8,'left') ' |' padd(sprintf('%.1fs',t(1)),7,'left') ' |'];
              locs(k,3)=1;
              locs(k,5)=0;
           end    
            
           dstr=[dstr padd(sprintf('%.1fs',t(i)),7,'left') ' |'];    
           
           locs(k,2)=i;
           wt=i; % keeper of the track of ending
         end 
         
         % if signal is cut from the end
         
         if l-lopp<=fs*0.1 && ~isempty(locs) && locs(k,2)==0 
             locs(k,2)=i;
             wt=i;
             dstr=[dstr padd(sprintf('%.1fs',t(i)),7,'left') ' |'];    
         end
         
         if wt>0 & i==wt+w  % abiakna lõpus
           locs(k,4)=i;  % locs(l,2)+5
           locs(k,6)=sum(acc_en(i-w+dw:i));  % spectral energy to the right
           wt=0;
           if k>1 & locs(k,2) & ~locs(k-1,4)
               %group
               if locs(k-1,5)>locs(k,6)
                  spd=calcsp(spctr(locs(k-1,3):locs(k-1,3)+w-dw,:));
                  dstr=[dstr padd(sprintf('%s','left'),10,'left') ' |' padd(sprintf('%.1fkm/h',f(spd)/17.3148),10,'left') ' '];      
               else
                  spd=calcsp(spctr(locs(k,4)-w+dw:locs(k,4),:));
                  dstr=[dstr padd(sprintf('%s','right'),10,'left') ' |' padd(sprintf('%.1fkm/h',f(spd)/17.3148),10,'left') ' ']; 
               end
               locs(k,7)=spd;
           end
         
         if locs(k,3) & locs(k,4)
           if locs(k,5)>locs(k,6)
                  spd=calcsp(spctr(locs(k,3):locs(k,3)+w-dw,:));
                  dstr=[dstr padd(sprintf('%s','left'),10,'left') ' |' padd(sprintf('%.1fkm/h',f(spd)/17.3148),10,'left') ' '];      
           else 
                  spd=calcsp(spctr(locs(k,4)-w+dw:locs(k,4),:));
                  dstr=[dstr padd(sprintf('%s','right'),10,'left') ' |' padd(sprintf('%.1fkm/h',f(spd)/17.3148),10,'left') ' ']; 
           end
           locs(k,7)=spd;
         end
        end
        end

        alg=alg+(win-overlap);
        lopp=alg+win-1;
        i=i+1;
        
    end
    disp(dstr)  % write the last result to the table
end

function  spd=calcsp(spctr)

 %finds the peak frequencies in the interval [lbnd,hbnd]
 %spctr=10*log10(spctr);
 [~,mm]=max(spctr');

 % gets a maximum of those
 spd=max(mm);

end

function [f,F]=spekter(signal,fs,alg,lopp)

% function [f,F]=spktr(signal,fs,alg,lopp)
%
% signal - raw signal
% alg - initial position [samples]
% lopp - end position [samples]
% fs - sampling frequency 
%  
% f - frequency vector [Hz] 
% F - spectral energy vector 

    if mod((lopp-alg),2)
        lopp=lopp-1;
    end

    signal=signal(alg:lopp);
    
    % TODO check if windowing helps
    %wd = hanning(length(signal));
    %signal = signal .* wd;
    
    y=fft(signal);
    N=size(y,1);

    L=(N-1)/fs; % size of the slice in seconds
    F=abs(y);
    F=fftshift(F)/N;
    f=[-L*fs/2:L*fs/2]/L;

    F=F((N-1)/2+2:N);
    f=f((N-1)/2+2:N);
    
end

function s=padd(s,n,suund)
   if strcmp(suund,'left')
    s=[blanks(n-length(s)) s];
   else
    s=[s blanks(n-length(s))];
   end
end    