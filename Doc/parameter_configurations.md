# Parameter Settings 
The parameter settings can be done in the file "config/config.txt". Each parameter is preceded by a line starting with "#" to introduce it.

## MPI settings
HPCrack utilizes 3D MPI topology configurations. The computational domain also employs a 3D decomposition. To ensure proper decomposition of the problem, the user must specify the number of MPI tasks for each direction. The total number of MPI tasks is the product of these three values.

*e.g.*, 
```
np={2,2,2};
```
means 8 MPI tasks are required. 

To disable MPI (and run the code only in OpenMP or serial version), the user should specify:
```
np={1,1,1};
```


