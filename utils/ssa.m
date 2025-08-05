function r = ssa(s, fps, K, hr_filt)

    if nargin < 3
        K = 2:5;
    end

    if nargin < 4
        hr_filt = false;
    end
    % step 1: embedding
    L = fps * 1; % length of embeding signal (1 sec)
    e = zeros(size(s,2) - L + 1, L);
    for i = 1 : size(s,2) - L + 1
        e(i,:) = s(i:i+L-1);
        e(i,:) = e(i,:)-mean(e(i,:));
    end

    % step 2: decomposition (SVD)
    [u,~,~] = svd(e);

    if hr_filt
        hr_filt_k = [];
        u1 = u(:,1);
        PR1 = prpsd(u1, fps, 42, 150, false);
        u2 = u(:,2);
        PR2 = prpsd(u2, fps, 42, 150, false);
        if abs(PR1-PR2) < 5
            hr_filt_k = [hr_filt_k, 1, 2]; % Add 1-2 PCs

            [f3, Y3] = simple_fft(u(:,3), fps);
            [~, ind3]= max(Y3);
            PR3 = f3(ind3)*60;
            if abs(PR3-(PR1+PR2)/2*2)<10 % Add 3-4 if around first harmonics
                hr_filt_k = [hr_filt_k, 3];
            end

            [f4, Y4] = simple_fft(u(:,4), fps);
            [~, ind4]= max(Y4);
            PR4 = f4(ind4)*60;
            if abs(PR4-(PR1+PR2)/2*2)<10
                hr_filt_k = [hr_filt_k, 4];
            end

            [f5, Y5] = simple_fft(u(:,5), fps); % in case 3rd PC is noise, 4-5 would be used
            [~, ind5]= max(Y5);
            PR5 = f5(ind5)*60;
            if abs(PR5-(PR1+PR2)/2*2)<10
                hr_filt_k = [hr_filt_k, 5];
            end

            % [f5, Y5] = simple_fft(u(:,5), fps);
            % [~, ind5]= max(Y5);
            % PR5 = f5(ind5)*60;
            % 
            % [f6, Y6] = simple_fft(u(:,6), fps);
            % [~, ind6]= max(Y6);
            % PR6 = f6(ind6)*60;
            % if abs(PR5-(PR1+PR2)/2*3)<10 % Add 5-6 if around second harmonics
            %     hr_filt_k = [hr_filt_k, 5];
            % end
            % if abs(PR6-(PR1+PR2)/2*3)<10
            %     hr_filt_k = [hr_filt_k, 6];
            % end
        else
            hr_filt_k = [1,2]; % Otherwise, use 1-4 PCs (can be improved) 
        end
        bp = u(:,hr_filt_k)*(u(:,hr_filt_k)'*e);
    else
        bp = u(:,K)*(u(:,K)'*e);        
    end

    % step 3: reconstruction (select components)

    % step 4: diagonal averaging
    r = zeros(1, size(s,2));
    d = zeros(1, size(s,2));

    for i = 1 : size(bp,1)
        r(1, i:i+L-1) = r(1, i:i+L-1) + bp(i,:);
        d(1, i:i+L-1) = d(1, i:i+L-1) + 1;
    end
    r = r./d;

end