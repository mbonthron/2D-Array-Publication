%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 1)
BP = coco_bd_labs('vtriangle1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','vtriangle1',BP(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'vtriangle1-2',[],1,parameter_names,computational_domain)