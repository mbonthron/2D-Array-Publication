%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 1)
clc
BP = coco_bd_labs('simplechain1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','simplechain1',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'simplechain1-1',[],1,parameter_names,computational_domain)