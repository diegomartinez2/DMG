using DelimitedFiles

function File_Sort_Function(filename_sp, smear_cm, nsm)
   for ism in 1:nsm
       name = "1.00"
       filename_data = string(filename_sp, '_', name, '.dat')
       data = readdlm(filename_data)
       X = data[:,1]
       Y = data[:,2]
       data_out = data[sortperm(Y, by = x -> (x[1], x[0])), :] # Here is the sorting part
       writedlm(filename_data, data_out)
   end
end

function main(args)
   smear_cm = 1.00
   File_Sort_Function("test", smear_cm, 1)
   return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
   main(ARGS)
end
