
**FLOSIC_CUBE2020**
------------
Cube utility generates cube files as a post-processing tool to the 
FLOSIC code.

**Brief instruction**
------------

To use cube utility:
- Compile the code:  
 Make -f Makefile.hpcc  
 Make -f Makefile.hpcc cube  

- Mark the get_cube.py file as executable:  
 chmod +x get_cube.py  

- Load necessary modules:
ml -* iccifort/2019.5.281 impi/2018.5.288 SciPy-bundle/2019.10-Python-3.7.4  
module load matplotlib  
export OMP_NUM_THREADS=8  

**Support** 
------------
For any question regarding the cube utility,  
please contact peral1j@cmich.edu (Dr. Peralta) or nguye8t@cmich.edu
