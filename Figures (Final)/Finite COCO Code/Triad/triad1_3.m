%% ========================================================================
%  BRANCH FROM FIRST RUN - Branch Point (No. 1)
clc
BP = coco_bd_labs('triad1', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triad1',BP(3));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triad1-3',[],1,parameter_names,computational_domain)

%% Since we detected HB points - rewrite those as UZ to continue from
bd = coco_bd_read('triad1-3');
HBlbls = coco_bd_labs('triad1-3', 'HB');

bcrits = zeros(1,length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end

prob = coco_add_event(prob,'UZ','b',bcrits(1));

coco(prob,'triad1-3',[],1,parameter_names,computational_domain)

%% Continue from HB points
clc
UZ = coco_bd_labs('triad1-3', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triad1-3',UZ(5));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triad1-3-1',[],1,parameter_names,computational_domain)

UZ = coco_bd_labs('triad1-3', 'UZ'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triad1-3',UZ(6));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triad1-3-2',[],1,parameter_names,computational_domain)

BP = coco_bd_labs('triad1-3-2', 'BP'); % labels for SN points in run1

prob = coco_prob();
prob = ode_ep2ep(prob,'','triad1-3-2',BP(1));
prob = coco_set(prob,'cont','branch','switch');
prob = coco_set(prob,'cont','ItMX', iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',h_max,'h_min',h_min);
prob = coco_add_event(prob,'UZ','b',UZpoint);

coco(prob,'triad1-3-3',[],1,parameter_names,computational_domain)