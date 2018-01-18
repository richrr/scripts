#Reads file and created appropriate vector of packages
file1 = "package_input.txt"
rf1 = read.table(file1)
print(rf1)
lop = as.vector(rf1[,1])
print(lop)

#Displays user installed packages 
ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
print(ip, row.names=FALSE)
ep = as.vector(ip[,1])
print(ep)

#Displays packages that are not already installed
pdne = lop[!(lop %in% ep)]
print(pdne)

#Installs new packages
#https://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
if(dir.exists(file.path("~/R/x86_64-redhat-linux-gnu-library", "3.3")) == FALSE){
       dir.create("~/R/x86_64-redhat-linux-gnu-library/3.3", recursive = TRUE)}
#if(length(pdne)) 
install.packages(pdne, destdir = "~/R/x86_64-redhat-linux-gnu-library/3.3")

#Checks if new packages were installed 
ap <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ap) <- NULL
ap <- ap[is.na(ip$Priority),1:2,drop=FALSE]
print(ap, row.names=FALSE)
aep = as.vector(ip[,1])
print(aep)
print(pdne %in% aep)
