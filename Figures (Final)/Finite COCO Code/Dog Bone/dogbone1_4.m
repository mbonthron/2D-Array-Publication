%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 4)
clc
BP = coco_bd_labs('dogbone1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','dogbone1',BP(4));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'dogbone1-4',[],1,parameter_names,computational_domain)

%%
BP = coco_bd_labs('dogbone1-4', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','dogbone1-4',BP(1));
prob = coco_set(prob,'cont','branch','switch'); prob = coco_set(prob,'cont','ItMX', iterations_max); prob = coco_set(prob,'cont','NPR',0); prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min); prob = coco_add_event(prob,'UZ','b',UZpoint); 

coco(prob,'dogbone1-4-1',[],1,parameter_names,computational_domain)

% 
BP = coco_bd_labs('dogbone1-4', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','dogbone1-4',BP(2));
prob = coco_set(prob,'cont','branch','switch'); prob = coco_set(prob,'cont','ItMX', iterations_max); prob = coco_set(prob,'cont','NPR',0); prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min); prob = coco_add_event(prob,'UZ','b',UZpoint); 

coco(prob,'dogbone1-4-2',[],1,parameter_names,computational_domain)

%
BP = coco_bd_labs('dogbone1-4', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','dogbone1-4',BP(4));
prob = coco_set(prob,'cont','branch','switch'); prob = coco_set(prob,'cont','ItMX', iterations_max); prob = coco_set(prob,'cont','NPR',0); prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min); prob = coco_add_event(prob,'UZ','b',UZpoint); 

coco(prob,'dogbone1-4-3',[],1,parameter_names,computational_domain)

prob = coco_prob();
prob = ode_ep2ep(prob,'','dogbone1-4',BP(5));
prob = coco_set(prob,'cont','branch','switch'); prob = coco_set(prob,'cont','ItMX', iterations_max); prob = coco_set(prob,'cont','NPR',0); prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min); prob = coco_add_event(prob,'UZ','b',UZpoint); 

coco(prob,'dogbone1-4-4',[],1,parameter_names,computational_domain)