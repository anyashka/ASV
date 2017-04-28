function [D_out] = vl_asvcovdet_modified(im, opt, frames_ori,des,isInter,threshold)

nr = opt.nr;
rc_min = opt.rc_min;
rc_max = opt.rc_max;
D_out = [];
frames = frames_ori;
nf = size(frames_ori, 2);
k = -0.2; %constant for threshold calculation
for rc = linspace(rc_min, rc_max, nr)
    
    frames([1 2], :) = frames_ori([1 2], :);
    for rf = 1:nf
        frames(3:6, rf) = reshape([cos(rc) -sin(rc);sin(rc) cos(rc)] * reshape(frames_ori(3:6,rf),2,2),4,1);
    end
    
    ns = opt.ns;
    sc_min = opt.sc_min;
    sc_max = opt.sc_max;
    % Sample scales around detection
    
    f = zeros(6, nf, ns);
    cnt = 0;
    for sc = linspace(sc_min, sc_max, ns)
        cnt = cnt + 1;
        f([1 2], :, cnt) = frames([1 2], :);
        f(3:6, :, cnt) = sc * frames(3:6,:);
    end
    
    
    if strcmp(des,'sift') == 1
        dim = 128;
    elseif strcmp(des,'liop') == 1
        dim = 144;
    elseif strcmp(des,'patch') == 1
        dim = 1681;
    end
    D = zeros(dim,nf,ns);
    for i = 1:ns
        [~,d] = extract(im,des,f(:,:,i));
        D(:,:,i) = d;
    end
    
    D = permute(D,[1,3,2]);
    
    % ------------------------------interpolation -----------------------------------
    if isInter == 1
        for in = 1:size(D,2)-1
            tempIntep = (D(:,in,:) + D(:,in+1,:))/2;
            D = cat(2,D,tempIntep);
        end
    end
    
    % ------------------------------interpolation -----------------------------------
    
    d_out = [];
    for i_f = 1:nf
        accVec = zeros(dim,1);
        
        for c = 1:size(D,2)-1
            
            temp = D(:,c,i_f);
            
            M = bsxfun(@minus,D(:,(c+1):end,i_f),temp);
            M = abs(M);
            %% 1st stage median thresholding
            m = median(M,1);
            deviation = std(M,1);
            
            for t = 1:size(m,2)
                
              switch(threshold)
                  case ThresholdType.Median
                    accVec = accVec + (M(:,t) <=  m(t)); 
                    
                  case ThresholdType.Niblack
                    niblack = m(t) + (k * deviation(t));
                    accVec = accVec + (M(:,t) <=  niblack);
                    
                  case ThresholdType.Sauvola
                    sauvola = m(t) * (1 + k * (deviation(t) / 128 - 1));
                    accVec = accVec + (M(:,t) <=  sauvola);  
              end
            end
            
        end
        
        
        
        accVec = double(accVec);
        d_out = [d_out,accVec];
    end
    
    D_out = [D_out;d_out];
end
end

