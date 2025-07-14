%material = 1;
for material= 1:2
if ~isfile("Analytical Prediction "+num2str(material) +".mat")
    snap_through_analytical(material);
end

load("Analytical Prediction "+num2str(material) +".mat");

end