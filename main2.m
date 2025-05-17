clear all
close all
% mex cec17_func.cpp -DWINDOWS
Xmin=-100;
Xmax=100;
pop_size=100;
iter_max=500;
runs=2;
fhd=str2func('cec17_func');
alg={'KATSA','TSA'};
alg2={'KATSA','TSA'};
for D=[100]
    disp(['D=',num2str(D)])
    for i=1:1%function number
        disp(['fun=',num2str(i)]);
        func_num=i;
        for j=1:runs
            [KATSA(j).gbest,KATSA(j).gbestval,KATSA(j).con,KATSA(j).time]= KATSA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
            [TSA(j).gbest,TSA(j).gbestval,TSA(j).con,TSA(j).time]= TSA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
        end
        for j=1:runs
            for an=1:size(alg,2)
                eval([alg{an} '_con(j,:)=' alg{an} '(j).con;']);
            end
        end
        for an=1:size(alg,2)
            eval([alg{an} '_r(i).conmean = mean(' alg{an} '_con);']);%convergence curve
            eval([alg{an} '_r(i).mean = mean([' alg{an} '.gbestval]);'])%mean result
            eval([alg{an} '_r(i).std = std([' alg{an} '.gbestval]);'])
            eval([alg{an} '_r(i).worst = max([' alg{an} '.gbestval]);'])
            eval([alg{an} '_r(i).best = min([' alg{an} '.gbestval]);'])
            eval([alg{an} '_r(i).time = mean([' alg{an} '.time]);'])
            eval(['meanR_T(i,an)= ' alg{an} '_r(i).mean ;']) %#ok<EVLEQ> %result table
        end
        figure(i)
        hold on
        for an=1:size(alg,2)
            if isequal(alg{an},'KATSA')
                semilogy(KATSA_r(i).conmean,'Marker','^','MarkerIndices',1:50:iter_max,'markersize',3,'LineWidth',1,'Color','r');
            else
                eval(['semilogy(' alg{an} '_r(i).conmean,"LineWidth",1)'])
            end
        end
        legend(alg2,'Location','bestoutside')
        title(['Convergence curve of F',num2str(i)])
        xlabel('Iteration');
        ylabel('Best score obtained so far');
        drawnow
    end
    % save(['data_17_','D',num2str(D)])% save your result if you want
end
disp('program finished');

% for i=1:29
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);i,f(i)
% end