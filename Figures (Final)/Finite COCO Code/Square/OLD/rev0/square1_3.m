%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 3)
clc
BP = coco_bd_labs('square1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1',BP(3));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-3',[],1,parameter_names,computational_domain)