%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 2)
clc
BP = coco_bd_labs('square1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1',BP(2));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-2',[],1,parameter_names,computational_domain)

% Since we detected HB points - rewrite those as UZ to continue from
bd = coco_bd_read('square1-2');
HBlbls = coco_bd_labs('square1-2', 'HB');

bcrits = zeros(1,length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end

prob = coco_add_event(prob,'UZ','b',bcrits(1));

coco(prob,'square1-2',[],1,parameter_names,computational_domain)

%%
UZ = coco_bd_labs('square1-2', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-2',UZ(5));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-2-1',[],1,parameter_names,computational_domain)

%%
UZ = coco_bd_labs('square1-2', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','square1-2',UZ(6));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'square1-2-2',[],1,parameter_names,computational_domain)