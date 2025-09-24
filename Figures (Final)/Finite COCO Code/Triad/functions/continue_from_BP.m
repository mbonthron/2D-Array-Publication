function [] = continue_from_BP(root_name,BPnum,extension_name,data)
    BP = coco_bd_labs(root_name, 'BP'); % labels for SN points in run1

    prob = coco_prob();
    prob = ode_ep2ep(prob,'',root_name,BP(BPnum));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', data.iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);
    prob = coco_add_event(prob,'UZ','b',data.UZpoint);
    
    coco(prob,extension_name,[],1,data.parameter_names,data.computational_domain)

end

