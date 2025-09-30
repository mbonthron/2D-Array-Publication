function [prob] = start_coco_from_another_point(Avalue,bvalue,run_name,data)
f = @(x,p) COCO_arbitrary_grid_ODE(x,p,data);

prob = coco_prob();
prob = ode_isol2ep(prob,'',f,Avalue,data.parameter_names,[bvalue ; data.initial_parameter_values(2)]);
prob = coco_set(prob,'corr','ItMX',25); 

prob = coco_set(prob,'cont','ItMX', data.iterations_max);
prob = coco_set(prob,'cont','NPR',0);
prob = coco_set(prob,'cont','h_max',data.h_max,'h_min',data.h_min);
prob = coco_add_event(prob,'UZ','b',data.UZpoint);

coco(prob,run_name,[],1,data.parameter_names,data.computational_domain)

end