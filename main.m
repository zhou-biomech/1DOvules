% Codes for the paper:
% Shi-Xia Yu, Lv-Wen Zhou and et al. Asynchrony of ovule primordia initiation 
% in Arabidopsis. Development 2020 : dev.196618 doi: 10.1242/dev.196618. 
% https://dev.biologists.org/content/early/2020/11/23/dev.196618
%
% Zhou Lvwen: zhoulvwen@nbu.edu.cn

clear; clc
figure('position',[50,50,1600,100])
ovule = Ovule(1:1:50);
ovule.plot(0)
axis([-2,404,-2,2])

caxis([0,2])
dt = 1e-3;
colorbar
for t = 0:dt:190
    ovule.grow;
    ovule.divde;
    ovule.auxin;
    if mod(t,1)==0
        ovule.plot(t);
        drawnow
    end
end


