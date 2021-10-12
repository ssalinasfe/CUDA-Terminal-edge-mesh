#include <vector> 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "consts.h"

#define filespath "input/"
#define filespathoutput "output/"


void read_from_triangle(std::string name, int &pnumber, int &tnumber, double *&points, int *&triangles, int *&neigh, int *&trivertex){
    std::string line;
    std::ifstream nodefile(name + ".node");
    double a1, a2, a3, a4;
    int i = 0;
    
    //std::cout<<"Node file"<<std::endl;
    if (nodefile.is_open())
    {
        nodefile >> pnumber ;
        //std::cout<<pnumber<<std::endl;

        std::getline(nodefile, line); 
        points = (double *)malloc(2*pnumber*sizeof(double));
        while (nodefile >> a1 >> a2 >> a3 >> a4)
        {
            points[2*i + 0] = a2;
            points[2*i + 1] = a3;
            //std::cout<<points[2*i + 0]<<" "<<points[2*i + 1]<<std::endl;
            i++;
            //std::cout<<a2<<" "<<a3<<std::endl;
            
        }
        
    }
    else 
        std::cout << "Unable to open node file"; 

    nodefile.close();


    //std::cout<<"Ele file"<<std::endl;
    std::ifstream elefile(name + ".ele");
    int t1, t2, t3, t4;
    i = 0;
    if(elefile.is_open()){
        elefile >> tnumber ;
        triangles = (int *)malloc(3*tnumber*sizeof(int));
        std::getline(elefile, line); 
        while (elefile >> t1 >> t2 >> t3 >> t4 )
        {
            //std::cout<<t2<<" "<<t3<<" "<<t4<<std::endl;
            triangles[3*i + 0] = t2;
            triangles[3*i + 1] = t3;
            triangles[3*i + 2] = t4;
            //std::cout<<triangles[3*i + 0]<<" "<<triangles[3*i + 1]<<" "<<triangles[3*i + 2]<<std::endl;
            i++;
        }
    }else std::cout << "Unable to open ele file";

    elefile.close();

    //std::cout<<"Neigh file"<<std::endl;
    std::ifstream neighfile(name + ".neigh");
    i = 0;
    if(neighfile.is_open()){
        std::getline(neighfile, line); 
        neigh =(int *)malloc(3*tnumber*sizeof(int));
        while (neighfile >> t1 >> t2 >> t3 >> t4 )
        {
            neigh[3*i + 0] = t2;
            neigh[3*i + 1] = t3;
            neigh[3*i + 2] = t4;
            //std::cout<<t2<<" "<<t3<<" "<<t4<<std::endl;
            i++;
        }
    }else std::cout << "Unable to open neigh file";
    neighfile.close();

    //std::cout<<"Neigh file"<<std::endl;
    std::ifstream trivertexfile(name + ".trivertex");
    i = 0;
    if(trivertexfile.is_open()){
        std::getline(trivertexfile, line); 
        trivertex =(int *)malloc(pnumber*sizeof(int));
        while (trivertexfile >> t1 >> t2)
        {
            trivertex[t1] = t2;
            i++;
        }
    }else std::cout << "Unable to open neigh file";
    trivertexfile.close();
}

/*geomview output*/
void write_geomview(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){

    int i,j;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    strcat(cmd,".off");
    
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr, "{ appearance  {+edge +face linewidth 2} LIST\n");
    fprintf(fptr, "OFF\n");
    fprintf(fptr,"%d %d 0\n", pnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.3f %.3f 0\n", r[2*i + 0], r[2*i + 1]);

  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]);
            i++;
        }
        fprintf(fptr, "\n");
    }

    if(print_triangles){
        
        fprintf(fptr, "{ appearance  {+edge -face linewidth 2} LIST\n");
        int p0, p1,p2;
        for(i = 0; i < tnumber; i++){
            p0 = 3*i + 0;
            p1 = 3*i + 1;
            p2 = 3*i + 2;
            fprintf(fptr,"# %d %d\n", p0, p1);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p1]+0], r[2*triangles[p1]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p1, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p1]+0], r[2*triangles[p1]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p0, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.3f %.3f 0\n%.3f %.3f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");
        }
    }

    fprintf(fptr," }\n");
    if(print_triangles){
        fprintf(fptr," }\n");
        fprintf(fptr," }\n");
    }
    /*
    std::sort( border_point.begin(), border_point.end() );
    border_point.erase( std::unique( border_point.begin(), border_point.end() ), border_point.end() );
    fprintf(fptr,"#Border vertices\n#");
    for ( i=0; i<border_point.size(); i++)
    {
        fprintf(fptr,"%d ", border_point[i]);
    }
    */
    fprintf(fptr,"\n");
    fclose(fptr);
}
