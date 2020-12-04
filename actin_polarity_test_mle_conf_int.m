%F = [1 1 1 1 1 0];

%F = [1 1 1 1 1];

%F = [1 1 0];

%F = [1 0 1 1 1 1];
F = [1 ];

%F = [1 1 1 1 1 1 1 0 0];

p_hat = mean(F);

std_F = sqrt(    p_hat*(1-p_hat)./length(F)    );

CI_l = p_hat - 1.96*std_F;

CI_r = p_hat + 1.96*std_F;

CI = [CI_l CI_r];

CI













% a = 3;
% b = 50;
% 
% c = 0.5;
% d = 1.0;
% 
% res = zeros(3,1000000);
% 
% for k=1:1000000
%     
%     l = (b-a).*rand(1) + a;
%     
%     p = (c-d).*rand(1) + d;
%     
%     filament_length = round(l);
%     
%     num_of_1 = ceil(p .* l);
%     
%     filament = zeros(1,filament_length);
%     filament(1,1:num_of_1) = 1;
%     
%     F = filament;
%     p_hat = mean(F);
%     
%     std_F = sqrt(p_hat*(1-p_hat)./length(F));
%     CI_l = p_hat - 1.96*std_F;
%     CI_r = p_hat + 1.96*std_F;
%     CI = [CI_l CI_r];
%     
%     res(1,k) = filament_length;
%     res(2,k) = p_hat;
%     res(3,k) = CI_l;
%     
%     %disp(filament)
%     
% end
% 
% 
% return;
    

    

% %F = [1 1 1 1 1 0];
% 
% %F = [1 1 1 1 1];
% 
% %F = [1 1 0];
% 
% %F = [1 1 1 0 0];
% 
% %F = [1 1 1 1 1 1 1 0 0];
% 
% 
% F = zeros(1,1067);
% F(1:556) = 1;
% 
% p_hat = mean(F);
% 
% 
% 
% std_F = sqrt(    p_hat*(1-p_hat)./length(F)    );
% 
% CI_l = p_hat - 1.96*std_F;
% 
% CI_r = p_hat + 1.96*std_F;
% 
% CI = [CI_l CI_r];
% 
% CI
% 













% 
% 
% 
% 
% i = 1;
% 
% filament_struct_ref(i).spin_filament_conf = p_hat;
% 
% filament_struct_ref(i).num_of_seg = length(F);
% 
% % Calculate maximum likelihood estimate of the filaments polarity
% filament_struct_ref(i).filament_mle_polarity = filament_struct_ref(i).spin_filament_conf;
%         
% % Calculate standard error of the maximum likelihood estimate
% filament_struct_ref(i).filament_std_error_polarity = sqrt((filament_struct_ref(i).filament_mle_polarity.*(1-filament_struct_ref(i).filament_mle_polarity))./filament_struct_ref(i).num_of_seg);
% 
% % Calculate 95% interval for the maximum likelihood intervall
% filament_struct_ref(i).filament_mle_conf_int_a = filament_struct_ref(i).filament_mle_polarity - 1.96.*filament_struct_ref(i).filament_std_error_polarity;
% filament_struct_ref(i).filament_mle_conf_int_b = filament_struct_ref(i).filament_mle_polarity + 1.96.*filament_struct_ref(i).filament_std_error_polarity;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
