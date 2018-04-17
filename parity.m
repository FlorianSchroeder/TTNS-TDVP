function pa = parity(mps0,Vmat0,para)
pa_op=parityop(para);
pa =expectationvalue(pa_op,mps0,Vmat0,mps0,Vmat0);
end