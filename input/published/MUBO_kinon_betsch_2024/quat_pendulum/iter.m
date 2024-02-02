
iter_fail = [1070.3902, 365.8401, 456.961, 1622.5253, 2262.7679, 925.7382, 847.9203 , 325.8736 , 962.7068, 479.8777 , 300.4326 ,9756.4798 ,9005320.8143, 9006072.6917 , 3130521.0207 ,2562777.5902 ,1985208.5449 ,10102458.429, 28461410.9962 ,12412877442.1243];
iter_suc = [102.9578, 16.0597, 0.031365, 4.3846e-06, 3.3363e-12];

figure()
semilogy(iter_fail);
hold on
semilogy(iter_suc);
ylabel('norm_R')
xlabel('iter')
grid on

matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'num_iter.tikz', 'showInfo', false, 'floatformat', '%.4g');