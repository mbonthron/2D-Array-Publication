%%
color_matrix = [    0.1190    0.5472    0.6160
                    0.4984    0.1386    0.4733
                    0.9597    0.1493    0.3517
                    0.3404    0.2575    0.8308
                    0.5853    0.8407    0.5853
                    0.2238    0.2543    0.5497
                    0.7513    0.8143    0.9172
                    0.2551    0.2435    0.2858
                    0.5060    0.9293    0.7572
                    0.6991    0.3500    0.7537
                    0.8909    0.1966    0.3804
                    0.9593    0.2511    0.5678];

counter = 1;

b_array = ["3.75" "5.625" "7.5" "9.375" "11.25" "15"];
L_array = ["L = 75mm" "L = 100mm"];

figure(2); clf; hold on

for b = 1:length(b_array)
    for L = 1:length(L_array)
        file_names = dir("b = "+b_array(b) + "mm " + L_array(L) + "*.xlsx");
        
        if ~isempty(file_names)
            for i = 1:length(file_names)
                    A = readtable(file_names(i).name,ReadVariableNames=false);
                    position = A.Var2; 
                    force = A.Var3;
                
                    startidx = find(force>0.05);
                
                    position = position(startidx:end);
                    force = force(startidx:end);
                
                    if i == 1
                        plot(position-position(1),force,'color',color_matrix(counter,:),"DisplayName","b = "+b_array(b)+"mm "+L_array(L),'LineWidth',2)
                    else
                        plot(position-position(1),force,'color',color_matrix(counter,:),"HandleVisibility","off",'LineWidth',2)
                    end
    
            end
        end
        counter = counter + 1;
    end
end

legend()