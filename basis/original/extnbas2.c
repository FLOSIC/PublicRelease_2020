#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#define MAXSIZE 100

double alpha[MAXSIZE];
typedef struct{
  double coef;
  int type;
  int gaussian_index;
} coeff;
coeff coef_list[MAXSIZE];
int number_type_functions[MAXSIZE][2];
int function_info[4];
char ElementsLetterU[][3]={"X", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "UUT", "UUQ", "UUP", "UUH" };
char ElementsLetterL[][3]={"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh" };
char test_string[4];

void init_all(void);
void process_token(char *token,FILE *ptr);
void print_alphas(int size,FILE *ptr);
void print_coefs(int size,FILE *ptr);
void print_func_type(int size,FILE *ptr);
void sort_all(int alphas,int coeffsi,int functions);
void swap_push_coefs(int i,int j);
void swap_push_functions(int i,int j);
void print_isymgen_style(FILE *ptr);

void main(int argc,char *argv[])
{
   char line[80];
   char out_file_name[80];
   char current;
   char s[2] = " ";
   char *token;
   FILE *fptr,*outptr;
   int i,j;

   printf("Processing file:%s\n",argv[1]);
   fptr=fopen(argv[1],"r");
/* Initialize arrays */
   init_all();
/* Open out file */
   strcpy(out_file_name,argv[1]);
   i=0;
   current=out_file_name[i];
   while(current!='.'){
//     printf("Comparing %c i=%d\n",current,i);
     i++;
     current=out_file_name[i];
   }
   i++;
   out_file_name[i]='b';
   i++;
   out_file_name[i]='a';
   i++;
   out_file_name[i]='s';
   i++;
   out_file_name[i]='i';
   i++;
   out_file_name[i]='s';
   i++;
   out_file_name[i]='\0';
   printf("Output file:%s\n",out_file_name);
   outptr=fopen(out_file_name,"w");
/* Process definition file line by line */
   while (fgets(line,80, fptr)!=NULL){

   /* get the first token */
   token = strtok(line, s);
   
   /* walk through other tokens */
   while( token != NULL ) 
   {
      process_token(token,outptr);
      token = strtok(NULL, s);
   }
  }
  fclose(fptr);
  fclose(outptr);
}
int get_atom_number(char *token){
  int i,j,k;
//  printf( "In atom number %s ", token );
  for(i=0;i<144;i++){
     for(k=0;k<3;k++)
         test_string[k]=ElementsLetterU[i][k];
     j=strcmp(test_string,token);
     if(j==0)
       return i;
  }
  for(i=0;i<144;i++){
     for(k=0;k<3;k++)
         test_string[k]=ElementsLetterL[i][k];
     j=strcmp(test_string,token);
     if(j==0)
       return i;
  }
  return -1;
}
void process_token(char *token,FILE *ptr){
  static int atom_number,s_function,p_function,d_function;
  static int global_s,global_p_global_d;
  static int mode,number_functions,current_function,function_index;
  static int gaussians,coefficients,function_type;
  static bool process_sp,process_atom;
  double value;
  int k;

  if(strcmp(token,"\n")==0){
//     printf("dumping empty token\n");
     return;
  }
  if(token[strlen(token)-2] =='\n') token[strlen(token)-2] ='\0';
//  printf("process token %s in mode %d\n",token,mode);
  if(strcmp(token,"****\n")==0){
    mode=1;
  }
  switch (mode){
    case 1:
        if(strcmp(token,"****\n")==0){
//	   printf("found header, going to mode 2\n");
           mode=2;
           if(process_atom){
              process_atom=false;
//              fprintf(ptr,"Gaussians %d\n",gaussians);
	      if(atom_number!=-1){
                 fprintf(ptr,"%d \n",atom_number);
                 for(k=0;k<3;k++)
                   test_string[k]=ElementsLetterL[atom_number][k];
                 fprintf(ptr,"%s \n",test_string);
                 fprintf(ptr,"%d \n",gaussians);
                 function_info[0]=gaussians;
                 fprintf(ptr,"%d %d %d \n",s_function,p_function,d_function);
                 function_info[1]=s_function;
                 function_info[2]=p_function;
                 function_info[3]=d_function;
                 sort_all(gaussians,coefficients,function_index); 
                 print_alphas(gaussians,ptr);
//              print_func_type(function_index,ptr);
//              print_coefs(coefficients,ptr);
                 print_isymgen_style(ptr);
	      }
              s_function=0;
              p_function=0;
              d_function=0;
           }
        }
        else if(strcmp(token,"S")==0){
//           printf("Found S function\n");
           s_function++;
           function_type=1;
//           printf("function type %d\n",function_type);
           process_sp=false;
           mode=4;
        }
        else if(strcmp(token,"P")==0){
//           printf("Found P function\n");
           p_function++;
           function_type=2;
           process_sp=false;
//           printf("function type %d\n",function_type);
           mode=4;
        }
        else if(strcmp(token,"D")==0){
//           printf("Found D function\n");
           d_function++;
           function_type=3;
           process_sp=false;
//           printf("function type %d\n",function_type);
           mode=4;
        }
        else if(strcmp(token,"SP")==0){
//           printf("Found SP function\n");
           process_sp=true;
           function_type=4;
//           printf("function type %d\n",function_type);
           mode=8;
        }
        else if(strcmp(token,"F")==0){
	   printf("Found F function ignoring atom %d\n",atom_number);
	   atom_number=-1;
	}
        else if(strcmp(token,"G")==0){
	   printf("Found G function ignoring atom %d\n",atom_number);
	   atom_number=-1;
	}
        else{
          printf("something wrong here!!\n");
          printf("Token is:%s",token);
        }
        break;
    case 2:
        mode=2;
        atom_number=get_atom_number(token);
//        printf("Found element %d\n",atom_number);
//        fprintf(ptr,"Atom %d\n",atom_number);
//        fprintf(ptr,"%d \n",atom_number);
        gaussians=0;
        process_atom=true;
        coefficients=0;
        function_index=0;
        mode=3;
        break;
    case 3:
//        printf("passing through after found atom\n");
        mode=1;
        break;
    case 4:
        number_functions=atoi(token);
//        printf("number of functions %d\n",number_functions);
        mode=5;
        current_function=0;
        number_type_functions[function_index][1]=number_functions;
        number_type_functions[function_index][2]=function_type;
        function_index++; 
        break;
    case 5:
//        printf("passing trough 1.00\n");
        mode=6;
        break;
    case 6:
//        printf("Getting Alpha %s\n",token);
        alpha[gaussians]=atof(token);
        mode =7;
        break;
    case 7:
//        printf("Getting coefficient %s\n",token);
        coef_list[coefficients].coef=atof(token);
        coef_list[coefficients].gaussian_index=gaussians;
        switch(function_type){
           case 1:
                coef_list[coefficients].type=1;
                break;
           case 2:
                coef_list[coefficients].type=2;
                break;
           case 3:
                coef_list[coefficients].type=3;
                break;
        }
        gaussians++;
        current_function++;
        coefficients++;
        if(current_function==number_functions)
           mode=1;
        else
           mode=6;
        break;
    case 8:
        number_functions=atoi(token);
//        printf("number of SP functions %d\n",number_functions);
//        fprintf(ptr,"number of SP functions %d\n",number_functions);
        mode=9;
        s_function++;
        p_function++;
        current_function=0;
        number_type_functions[function_index][1]=number_functions;
        number_type_functions[function_index][2]=1;
        function_index++;
        number_type_functions[function_index][1]=number_functions;
        number_type_functions[function_index][2]=2;
        function_index++;
        break;
    case 9:
//        printf("passing trough 1.00\n");
        mode=10;
        break;
    case 10:
//        printf("getting Alpha for SP function %s\n",token);
        alpha[gaussians]=atof(token);
        mode=11;
        break;
    case 11:
//        printf("getting frst coefficient for SP function %s\n",token);
        coef_list[coefficients].coef=atof(token);
        coef_list[coefficients].type=1;
        coef_list[coefficients].gaussian_index=gaussians;
        coefficients++;
        mode=12;
        break;
    case 12:
//        printf("getting second coefficient for SP function %s\n",token);
        coef_list[coefficients].coef=atof(token);
        coef_list[coefficients].type=2;
        coef_list[coefficients].gaussian_index=gaussians;
        coefficients++;
        gaussians++;
        current_function++;
        if(current_function==number_functions)
           mode=1;
        else
           mode=10;
        break;
  }
}
// 
// Initialize arrays
//
void init_all(void){
  int i;
  for(i=0;i<MAXSIZE;i++){
      alpha[i]=0.0;
      coef_list[i].coef=0.0;
      coef_list[i].type=0;
      coef_list[i].gaussian_index=0;
      number_type_functions[i][1]=0;
      number_type_functions[i][2]=0;
 }
 for(i=0;i<4;i++)
    function_info[i]=0;
}
//
// Print Alphas
//
void print_alphas(int size,FILE *ptr){
 int i;
 for(i=0;i<size;i++)
     fprintf(ptr,"%1.8f ",alpha[i]);
 fprintf(ptr,"\n");
}
//
// Printf Coefficients
//
void print_coefs(int size,FILE *ptr){
 int i;
 for(i=0;i<size;i++)
    fprintf(ptr,"%1.10E ",coef_list[i].coef);
 fprintf(ptr,"\n");
 for(i=0;i<size;i++)
    fprintf(ptr,"%d ",coef_list[i].type);
 fprintf(ptr,"\n");
 for(i=0;i<size;i++)
    fprintf(ptr,"%d ",coef_list[i].gaussian_index);
 fprintf(ptr,"\n");
}
//
// Print function types
//
void print_func_type(int size,FILE *ptr){
 int i;
 for(i=0;i<size;i++)
//    fprintf(ptr,"%d %d ",number_type_functions[i][1],number_type_functions[i][2]);
    fprintf(ptr,"%d ",number_type_functions[i][1]);
 fprintf(ptr,"\n");
}
//
// Sort all arrays 
//
void sort_all(int alphas,int coeffs,int functions){
 int i,j,k;
 double value;
// printf("Gaussians=%d Coefficients=%d\n",alphas,coeffs);
// Sort alpha values
 for(i=0;i<alphas-1;i++){
    for(j=i+1;j<alphas;j++){
        if(alpha[i]<alpha[j]){
           value=alpha[j];
           alpha[j]=alpha[i];
           alpha[i]=value;
// Adjust coefficients accordingly
           for(k=0;k<coeffs;k++){
              if(coef_list[k].gaussian_index==i)
                coef_list[k].gaussian_index=j;
              else if(coef_list[k].gaussian_index==j)
                coef_list[k].gaussian_index=i;
           }
        }
    }
 }
//Sort coefficients according to function type
 for(i=0;i<coeffs-1;i++)
   for(j=i+1;j<coeffs;j++)
        if(coef_list[i].type>coef_list[j].type)
           swap_push_coefs(i,j);
// Sort function type numbers
 for(i=0;i<functions-1;i++)
   for(j=i+1;j<functions;j++)
      if(number_type_functions[i][2]>number_type_functions[j][2])
        swap_push_functions(i,j);
//
}
//
// Swap two coefficients by pushing out
//
void swap_push_coefs(int i,int j){
 int k;
 coeff temp;

 temp.coef=coef_list[j].coef;
 temp.type=coef_list[j].type;
 temp.gaussian_index=coef_list[j].gaussian_index;
 for(k=j;k>i+1;k--){
    coef_list[k].coef=coef_list[k-1].coef;
    coef_list[k].type=coef_list[k-1].type;
    coef_list[k].gaussian_index=coef_list[k-1].gaussian_index;
 }
 coef_list[i+1].coef=coef_list[i].coef;
 coef_list[i+1].type=coef_list[i].type;
 coef_list[i+1].gaussian_index=coef_list[i].gaussian_index;
 coef_list[i].coef=temp.coef;
 coef_list[i].type=temp.type;
 coef_list[i].gaussian_index=temp.gaussian_index;
}
//
// Swap two functions by pushing out
//
void swap_push_functions(int i,int j){
 int k;
 int temp[2];
 temp[1]=number_type_functions[j][1];
 temp[2]=number_type_functions[j][2];
 for(k=j;k>i+1;k--){
    number_type_functions[k][1]=number_type_functions[k-1][1];
    number_type_functions[k][2]=number_type_functions[k-1][2];
 }
 number_type_functions[i+1][1]=number_type_functions[i][1];
 number_type_functions[i+1][2]=number_type_functions[i][2];
 number_type_functions[i][1]=temp[1];
 number_type_functions[i][2]=temp[2];
}
//
// Print to file following ISYMGEN format
//
void print_isymgen_style(FILE *ptr){
 double entry_line[MAXSIZE];
 double expo,entry;
 int total_functions,i,j,m,block;

 block=0;
 m=0;
 total_functions=function_info[1]+function_info[2]+function_info[3];
// printf("Total functions:%d\n",total_functions);
// printf("Total Alphas:%d\n",function_info[0]);
 for(i=0;i<total_functions;i++){
    for(j=0;j<MAXSIZE;j++)
       entry_line[j]=0.0;
//    printf("doing block %d\n",block);
//    printf("this block has %d functions\n",number_type_functions[i][1]);
    for(j=0;j<number_type_functions[i][1];j++){
//        printf("j=%d\n",j);
//        printf("index here is %d\n",coef_list[m].gaussian_index);
        expo=0.75+((double)coef_list[m].type-1.0)/2.0;
        entry=coef_list[m].coef*pow(alpha[coef_list[m].gaussian_index],expo);
        entry_line[coef_list[m].gaussian_index]=entry;
//        if(coef_list[m].coef<0){
//           printf("type=%d coef=%E alpha=%E entry=%E\n",coef_list[m].type,coef_list[m].coef,alpha[coef_list[m].gaussian_index],entry);
//           printf("i=%d j=%d m=%d gaussian index=%d\n",i,j,m,coef_list[m].gaussian_index);
//           printf("%E\n",entry_line[coef_list[m].gaussian_index]);
//        }
        m++;
    }
    for(j=0;j<function_info[0];j++){
       fprintf(ptr,"%1.10E ",entry_line[j]);
//       printf("%E\n",entry_line[j]);
    }
    fprintf(ptr,"\n");
    block++; 
 }
}
