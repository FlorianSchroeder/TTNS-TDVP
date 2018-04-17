classdef TDVPLibHmult
	
	properties
		
	end
		
	methods
		
	end
		
	methods (Static = true)
		
		function w = HAAmultV(V, op, dv, M, parity)
			%% HAAmultV(V)
			% Needs previous calculation of op.HrightA, op.HleftA, op.OprightA, op.OpleftA
			% works with multi-chain model!

			w = HmultVmat(V, op, dv(1),dv(2), M,parity);
		end
		
	end
	
	
end
		
		
	