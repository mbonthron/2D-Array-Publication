%% ========================================================================
%  BRANCH FROM FIRST RUN - "HB" Point (No. 1)
clc
UZ = coco_bd_labs('square1', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1',UZ(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5',[],1,parameter_names,computational_domain)

%%
% Run from BP points (1)
BP = coco_bd_labs('square1-5', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-1',[],1,parameter_names,computational_domain)

% Run from BP points (2)
BP = coco_bd_labs('square1-5', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5',BP(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-2',[],1,parameter_names,computational_domain)

%% These two branch points don't contribute anything
% BP = coco_bd_labs('square1-5-2', 'BP'); % labels for SN points in run1
% prob = coco_prob();
% prob = ode_ep2ep(prob,'','square1-5',BP(1));
% prob = coco_set(prob,'cont','branch','switch');
% prob = coco_set(prob,'cont','ItMX', iterations_max);
% prob = coco_set(prob,'cont','NPR',0);
% prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
% prob = coco_add_event(prob,'UZ','b',UZpoint);
% 
% coco(prob,'square1-5-2-1',[],1,parameter_names,computational_domain)
% 
% BP = coco_bd_labs('square1-5-2', 'BP'); % labels for SN points in run1
% prob = coco_prob();
% prob = ode_ep2ep(prob,'','square1-5',BP(3));
% prob = coco_set(prob,'cont','branch','switch');
% prob = coco_set(prob,'cont','ItMX', iterations_max);
% prob = coco_set(prob,'cont','NPR',0);
% prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
% prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-2-2',[],1,parameter_names,computational_domain)

%%
% Run from BP points (3)
BP = coco_bd_labs('square1-5', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5',BP(3));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-3',[],1,parameter_names,computational_domain)


BP = coco_bd_labs('square1-5-3', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5-3',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-3-1',[],1,parameter_names,computational_domain)

BP = coco_bd_labs('square1-5-3', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5-3',BP(7));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-3-2',[],1,parameter_names,computational_domain)

%%
% Run from BP points (3)
BP = coco_bd_labs('square1-5', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5',BP(4));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-4',[],1,parameter_names,computational_domain)


BP = coco_bd_labs('square1-5', 'BP'); % labels for SN points in run1
prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-5',BP(5));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-5-5',[],1,parameter_names,computational_domain)