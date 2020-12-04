
% Construct base path
actin_polarity_construct_base_path;

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Process bundles
for k=[1,2,3,7]

    % Load filament structure
    load(FilamentStruct{k},'filament_struct_ref');
    
%     % Plot selected filaments
%     figure;hold on;xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);box on;
%     fig_sel_filaments = gcf;
%     figure(fig_sel_filaments);title(k);
%     cmap = colormap(jet(size(filament_struct_ref,2)));
%     for i=1:size(filament_struct_ref,2)
%          for j=1:size(filament_struct_ref(i).cor_filament,1)-1
%               plot3([filament_struct_ref(i).cor_filament(j,1) filament_struct_ref(i).cor_filament(j+1,1)],...
%                     [filament_struct_ref(i).cor_filament(j,2) filament_struct_ref(i).cor_filament(j+1,2)],...
%                     [filament_struct_ref(i).cor_filament(j,3) filament_struct_ref(i).cor_filament(j+1,3)],'-','Color',cmap(i,:));
%          end
%          plot3([filament_struct_ref(i).cor_filament(1,1) filament_struct_ref(i).cor_filament(size(filament_struct_ref(i).cor_filament,1),1)],...
%                [filament_struct_ref(i).cor_filament(1,2) filament_struct_ref(i).cor_filament(size(filament_struct_ref(i).cor_filament,1),2)],...
%                [filament_struct_ref(i).cor_filament(1,3) filament_struct_ref(i).cor_filament(size(filament_struct_ref(i).cor_filament,1),3)],'k--');
%              
%     end
    
%     % Plot spin plot
%     figure;hold on;box on;
%     fig_indx_spin_plot = gcf;
%     figure(fig_indx_spin_plot);title(k);
%     
%     spin_filament_mean_all = [filament_struct_ref.spin_filament_mean];
%     spin_filament_mean_all = round(spin_filament_mean_all.*100)./100;
%     [~,indx_spin_filament_mean_all_sort] = sort(spin_filament_mean_all);
%     
%     zaehler_f = 1;
%     
%     for i=indx_spin_filament_mean_all_sort
%         
%         %plot([zaehler_f zaehler_f],[1 size(filament_struct_ref(i).cor_filament,1)],'-','Color',[0.7 0.7 0.7],'LineWidth',2);
%         
%         for j=1:size(filament_struct_ref(i).cor_filament,1)
%               
%               if filament_struct_ref(i).spin_filament(j,1) == 0
%                   plot(zaehler_f,j,'bo','MarkerFaceColor','b','MarkerSize',2.5);
%               elseif filament_struct_ref(i).spin_filament(j,1) == 1
%                   plot(zaehler_f,j,'ro','MarkerFaceColor','r','MarkerSize',2.5);
%               end
%         end
%         
%         %text(zaehler_f,j+1,[num2str(     round(100.*filament_struct_ref(i).filament_conf_com_score)./100     ) ' | ' num2str(     round(4.*0.33.*filament_struct_ref(i).length_filament)     )],'Rotation',90,'FontSize',6,'HorizontalAlignment','left','VerticalAlignment','middle','FontWeight','bold');
%         
%         zaehler_f = zaehler_f + 1;
%          
%     end
%     
%     xlim([0 size(spin_filament_mean_all,2)+1]);
    
    % Plot decided filaments
    figure;hold on;axis square;xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);box on;
    fig_dec_filaments = gcf;
    figure(fig_dec_filaments);title(k);
    
    for i=1:size(filament_struct_ref,2)
        
        if filament_struct_ref(i).spin_filament_red_red == 0
         
            for j=1:size(filament_struct_ref(i).cor_filament,1)-1
              plot3([filament_struct_ref(i).cor_filament(j,1) filament_struct_ref(i).cor_filament(j+1,1)],...
                    [filament_struct_ref(i).cor_filament(j,2) filament_struct_ref(i).cor_filament(j+1,2)],...
                    [filament_struct_ref(i).cor_filament(j,3) filament_struct_ref(i).cor_filament(j+1,3)],'b-','LineWidth',1.2);
            end
            
        end
        
        if filament_struct_ref(i).spin_filament_red_red == 1
         
            for j=1:size(filament_struct_ref(i).cor_filament,1)-1
              plot3([filament_struct_ref(i).cor_filament(j,1) filament_struct_ref(i).cor_filament(j+1,1)],...
                    [filament_struct_ref(i).cor_filament(j,2) filament_struct_ref(i).cor_filament(j+1,2)],...
                    [filament_struct_ref(i).cor_filament(j,3) filament_struct_ref(i).cor_filament(j+1,3)],'r-','LineWidth',1.2);
            end
            
        end
        
    end
    view(0,-90);xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);
    
%     % Plot decided filaments with neigbourhoods
%     figure;hold on;xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);box on;
%     fig_eps_filaments = gcf;
%     figure(fig_eps_filaments);title(k);
%     
%     for i=1:size(filament_struct_ref,2)
%         
%         if filament_struct_ref(i).spin_filament_red == 0
%          
%             for j=1:size(filament_struct_ref(i).cor_filament,1)-1
%                 plot3([filament_struct_ref(i).cor_filament(j,1) filament_struct_ref(i).cor_filament(j+1,1)],...
%                       [filament_struct_ref(i).cor_filament(j,2) filament_struct_ref(i).cor_filament(j+1,2)],...
%                       [filament_struct_ref(i).cor_filament(j,3) filament_struct_ref(i).cor_filament(j+1,3)],'b-');
%             end
%             
%             for j=1:size(filament_struct_ref(i).cor_filament,1)
%                 plot3(filament_struct_ref(i).cor_filament(j,1),filament_struct_ref(i).cor_filament(j,2),filament_struct_ref(i).cor_filament(j,3),'b.');
%             end
%             
%         end
%         
%         if filament_struct_ref(i).spin_filament_red == 1
%          
%             for j=1:size(filament_struct_ref(i).cor_filament,1)-1
%                 plot3([filament_struct_ref(i).cor_filament(j,1) filament_struct_ref(i).cor_filament(j+1,1)],...
%                       [filament_struct_ref(i).cor_filament(j,2) filament_struct_ref(i).cor_filament(j+1,2)],...
%                       [filament_struct_ref(i).cor_filament(j,3) filament_struct_ref(i).cor_filament(j+1,3)],'r-');
%             end
%             
%             for j=1:size(filament_struct_ref(i).cor_filament,1)
%                 plot3(filament_struct_ref(i).cor_filament(j,1),filament_struct_ref(i).cor_filament(j,2),filament_struct_ref(i).cor_filament(j,3),'r.');
%             end
%             
%         end
%         
%         for m=1:size(filament_struct_ref(i).cor_filament,1)
%             
%             if size(filament_struct_ref(i).segment_neighbourhood_epsilon(m).spin_neighbourhood,2) == 0
%                 continue;
%             end
%             
%             for n=1:size(filament_struct_ref(i).segment_neighbourhood_epsilon(m).spin_neighbourhood,2)
%                 
%                 plot3([filament_struct_ref(i).cor_filament(m,1) filament_struct_ref(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_filament(1,n)   ).cor_filament(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_segment(1,n)   ,1)],...
%                       [filament_struct_ref(i).cor_filament(m,2) filament_struct_ref(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_filament(1,n)   ).cor_filament(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_segment(1,n)   ,2)],...          
%                       [filament_struct_ref(i).cor_filament(m,3) filament_struct_ref(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_filament(1,n)   ).cor_filament(   filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_segment(1,n)   ,3)],'k-');
%                   
%             end
%             
%         end
%         
%     end
%     view(0,-90);xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);
    
    disp(k);
    
end

