function run_VMPS_mult(p, eta)

for  alpha =
    for period = p
		try
			VMPS_FullSBM(1,alpha,0.1,0);
       	catch err
            getReport(err,'extended')
			fprintf('Problem with the combination: period=%.10g',period)
        end
    end
end

exit

end
