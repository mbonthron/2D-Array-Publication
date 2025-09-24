%% ========================================================================
%  BRANCH FROM FIRST RUN - "HB" Point (No. 2)
clc
UZ = coco_bd_labs('square1', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1',UZ(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-6',[],1,parameter_names,computational_domain)