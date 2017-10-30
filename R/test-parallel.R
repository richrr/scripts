#cd /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/merged/tmp_del_test
#SGE_Batch -c "Rscript test-parallel-v4.R 1000000" -m 100G -F 100G -r log4core1000Kv4 -q biomed -M rodrrich@oregonstate.edu -P 16
#SGE_Batch -c "Rscript test-parallel-v4.R 100000" -m 100G -F 100G -r log4core100Kv4 -q biomed -M rodrrich@oregonstate.edu -P 16
#SGE_Batch -c "Rscript test-parallel-v4.R 10000" -m 100G -F 100G -r log4core10Kv4 -q biomed -M rodrrich@oregonstate.edu -P 16



#### final summary ###
# using function does NOT slow down the parallelization (code)
# passing rownames to the apply function during inline or separate function call slows code down
# using 4 cores gives most efficiency
# using 8 cores is the fastest
# beyond 8 cores there is not much speed up


args = commandArgs(trailingOnly=TRUE)
rownumb = args[1]

# https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r

sink(paste("outlog", rownumb, "v4",  "txt", sep="."))                     # Begin writing output to file

          
rownumb = as.numeric(rownumb)

d1 = c(1,2,-3)
d2 = c(1,-2,-2)
d3 = c(1,-2,-4)
d = data.frame(d1,d2,d3)
d

#  d1 d2 d3
# 1  1  1  1
# 2  2 -2 -2
# 3 -3 -2 -4


s_df = d[rep(seq_len(nrow(d)), each=rownumb),]

dim(s_df)


correlThreshold = 0
total_numb_input_files = 3
rows_passing_consistency = c()

if(FALSE){
print("For loop all conditions")
system.time(
		for(idx in 1:nrow(s_df)){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            rows_passing_consistency = append(rows_passing_consistency, res)
		}
)



print("Parallel For loop all conditions")
library(doParallel)
cl<-makeCluster(4, type="FORK")
registerDoParallel(cl)
system.time(
        foreach(idx = 1:nrow(s_df), .combine = c) %dopar% {
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
            res
		}
)
stopCluster(cl)

}

#http://r.789695.n4.nabble.com/how-to-pass-more-than-one-argument-to-the-function-called-by-lapply-td895470.html
print("Lapply all conditions")
system.time(
  lapply(1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)



print("ParLapply all conditions 4 cores")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)


print("ParLapply all conditions 8 cores")
library(doParallel)
cl<-makeCluster(8, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)

######## does not give much more speed up than using 8 cores.
if(FALSE){
print("ParLapply all conditions 16 cores")
library(doParallel)
cl<-makeCluster(16, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)
}



############### passing rownames makes it slower ###############
calc_consistency = function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
			res
}


print("Lapply all conditions call function")
system.time(lapply(1:nrow(s_df), calc_consistency, s_df))

print("ParLapply all conditions call function 4 cores")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(parLapply(cl, 1:nrow(s_df), calc_consistency, s_df))
stopCluster(cl)


print("ParLapply all conditions call function 8 cores")
library(doParallel)
cl<-makeCluster(8, type="FORK")
#registerDoParallel(cl)
system.time(parLapply(cl, 1:nrow(s_df), calc_consistency, s_df))
stopCluster(cl)




########## inline function but passing rownames ######
#http://r.789695.n4.nabble.com/how-to-pass-more-than-one-argument-to-the-function-called-by-lapply-td895470.html
print("using rownames with inline function")
print("Lapply all conditions")
system.time(
  lapply(rownames(s_df) , 
  		function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			#rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)



print("ParLapply all conditions 4 cores")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, rownames(s_df) , 
  		function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			#rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)


print("ParLapply all conditions 8 cores")
library(doParallel)
cl<-makeCluster(8, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, rownames(s_df) , 
  		function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			#rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)




############### trying compliling didn't help ################
if(FALSE){
# https://www.r-bloggers.com/faster-higher-stonger-a-guide-to-speeding-up-r-code-for-busy-people/
library(compiler)
calc_consistency_uncompiled = function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
			res
}

calc_consistency <- cmpfun(calc_consistency_uncompiled)
}

if(FALSE){
calc_consistency = function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
			res
}


print("Lapply all conditions call function")
system.time(lapply(rownames(s_df), calc_consistency, s_df))

print("ParLapply all conditions call function 4 cores")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(parLapply(cl, rownames(s_df), calc_consistency, s_df))
stopCluster(cl)


print("ParLapply all conditions call function 8 cores")
library(doParallel)
cl<-makeCluster(8, type="FORK")
#registerDoParallel(cl)
system.time(parLapply(cl, rownames(s_df), calc_consistency, s_df))
stopCluster(cl)
}


############### for loops are slower ###############
if(FALSE){
print("For loop all conditions call function")
system.time(
for(rname in rownames(s_df)){
		rows_passing_consistency = append(rows_passing_consistency, calc_consistency(rname, s_df))
}
)

print("Parallel For loop all conditions call function")
library(doParallel)
cl<-makeCluster(4, type="FORK")
registerDoParallel(cl)
system.time(
foreach(rname=rownames(s_df), .combine = c) %dopar% {
		calc_consistency(rname, s_df)
}
)
stopCluster(cl)

}

sink() # Resume writing output to console





############# older version of this code (v3) #########

if(FALSE){
#cd /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/merged/tmp_del_test
# SGE_Batch -c "Rscript test-parallel-v3.R 1000000" -m 100G -F 100G -r log4core1000Kv3 -q biomed -M rodrrich@oregonstate.edu -P 16
# SGE_Batch -c "Rscript test-parallel-v3.R 100000" -m 100G -F 100G -r log4core100Kv3 -q biomed -M rodrrich@oregonstate.edu -P 16
# SGE_Batch -c "Rscript test-parallel-v3.R 10000" -m 100G -F 100G -r log4core10Kv3 -q biomed -M rodrrich@oregonstate.edu -P 16



args = commandArgs(trailingOnly=TRUE)
rownumb = args[1]

# https://stackoverflow.com/questions/2470248/write-lines-of-text-to-a-file-in-r

sink(paste("outlog", rownumb, "v3",  "txt", sep="."))                     # Begin writing output to file

          
rownumb = as.numeric(rownumb)

d1 = c(1,2,-3)
d2 = c(1,-2,-2)
d3 = c(1,-2,-4)
d = data.frame(d1,d2,d3)
d

#  d1 d2 d3
# 1  1  1  1
# 2  2 -2 -2
# 3 -3 -2 -4


s_df = d[rep(seq_len(nrow(d)), each=rownumb),]

dim(s_df)


correlThreshold = 0
total_numb_input_files = 3
rows_passing_consistency = c()

if(FALSE){
print("For loop all conditions")
system.time(
		for(idx in 1:nrow(s_df)){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            rows_passing_consistency = append(rows_passing_consistency, res)
		}
)



print("Parallel For loop all conditions")
library(doParallel)
cl<-makeCluster(4, type="FORK")
registerDoParallel(cl)
system.time(
        foreach(idx = 1:nrow(s_df), .combine = c) %dopar% {
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
            res
		}
)
stopCluster(cl)

}

#http://r.789695.n4.nabble.com/how-to-pass-more-than-one-argument-to-the-function-called-by-lapply-td895470.html
print("Lapply all conditions")
system.time(
  lapply(1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)



print("ParLapply all conditions")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)


print("ParLapply all conditions")
library(doParallel)
cl<-makeCluster(8, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)


print("ParLapply all conditions")
library(doParallel)
cl<-makeCluster(16, type="FORK")
#registerDoParallel(cl)
system.time(
  parLapply(cl, 1:nrow(s_df) , 
  		function(idx, s_df){
			dfx = as.vector(unlist(s_df[idx, ,drop=T]))
			rname = rownames(s_df)[idx]
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)

			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }

			if((!is.na(neg)) && neg == (total_numb_input_files)){
				res = rname
            }
            res 		
  		}
	, s_df
  )
)
stopCluster(cl)



if(FALSE){

calc_consistency = function(rname, s_df){
			dfx = as.vector(unlist(s_df[rname, ,drop=T]))
			res = NULL
			#print(dfx)
			pos = sum(dfx > correlThreshold)
			neg = sum(dfx < correlThreshold)
			
			if((!is.na(pos)) && pos == (total_numb_input_files)){
				res = rname
            }
			if((!is.na(neg)) && neg == (total_numb_input_files)){
				 res = rname
            }
			return(res)
}

print("Lapply all conditions call function")
system.time(lapply(rownames(s_df), calc_consistency, s_df))

print("ParLapply all conditions call function")
library(doParallel)
cl<-makeCluster(4, type="FORK")
#registerDoParallel(cl)
system.time(parLapply(cl, rownames(s_df), calc_consistency, s_df))
stopCluster(cl)

print("For loop all conditions call function")
system.time(
for(rname in rownames(s_df)){
		rows_passing_consistency = append(rows_passing_consistency, calc_consistency(rname, s_df))
}
)

print("Parallel For loop all conditions call function")
library(doParallel)
cl<-makeCluster(4, type="FORK")
registerDoParallel(cl)
system.time(
foreach(rname=rownames(s_df), .combine = c) %dopar% {
		calc_consistency(rname, s_df)
}
)
stopCluster(cl)

}

sink() # Resume writing output to console



}
##############








