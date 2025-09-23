function [] = continue_from_HB(root_name,UZnum,extension_name,data)
    UZ = coco_bd_labs(root_name, 'HB'); % labels for SN points in run1

    prob = coco_prob();
    prob = ode_ep2ep(prob,'',root_name,UZ(UZnum));
    prob = coco_set(prob,'cont','branch','switch');
    prob = coco_set(prob,'cont','ItMX', data.iterations_max);
    prob = coco_set(prob,'cont','NPR',0);
    prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);
    prob = coco_add_event(prob,'UZ','b',data.UZpoint);
    
    coco(prob,extension_name,[],1,data.parameter_names,data.computational_domain)

end

