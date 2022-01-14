/*
- No se pude implementar el remove BE secuencial porqué la cantidad de poligonos a guardar en la nueva malla, no hay forma de usar guardar n nuevos poligonos y la vez los triangulos semilla de estos en poligonos. Guardar ambos en registos sobrepasa la capacidad de registros por hilo
*/

#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
//#include "detri2.h"
//#include "polymesh.h"
#include <vector> 
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <algorithm>    // std::min
#include <string>


#include "consts.h"


//cuda
#include "io.cuh"

#include "triang_edge.cuh"
#include "polygon.cuh"
#include "BET_elimination.cuh"


#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

int main(int argc, char* argv[])
{

    //int print_triangles = 0;
	
	int tnumber, pnumber;
	double *r;
	int *triangles;
	int *adj;
    int *seed;
	int *max;
	int *mesh;
	int *ind_poly;
	int *trivertex;
	
	std::string name(argv[1]);
	std::string output(argv[2]);
	//std::cout<<name<<std::endl;
	read_from_triangle(name, pnumber, tnumber, r, triangles, adj, trivertex);
	//std::cout << " " << tnumber << " " << pnumber << "\n";

    //tnumber = Tr->tnumber;
    //pnumber = Tr->pnumber;
    //r = (double *)malloc(2*tnumber*sizeof(double));
    //adj =(int *)malloc(3*tnumber*sizeof(int));
    //triangles = (int *)malloc(3*tnumber*sizeof(int));
	max = (int *)malloc(tnumber*sizeof(int));
	seed = (int *)malloc(tnumber*sizeof(int));
	mesh = (int *)malloc(3*tnumber*sizeof(int));
	ind_poly = (int *)malloc(tnumber*sizeof(int));

	

	//Cuda functions
    // Initialize device pointers.
    double *cu_r;
	int *cu_triangles;
	int *cu_adj;
    int *cu_seed;
	int *cu_max;
	int *cu_mesh;
	int *cu_mesh_aux;
	int *cu_ind_poly;
	int *cu_ind_poly_aux;
	int *cu_trivertex;

	//std::cout << "Solicitando memoria cuda"<<std::endl;

	// Allocate device memory.
	gpuErrchk( cudaMalloc((void**) &cu_max,          tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_seed,         tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_r,            2*pnumber*sizeof(double)) );
	gpuErrchk( cudaMalloc((void**) &cu_triangles,    3*tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_adj,          3*tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_mesh,         3*tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_mesh_aux,     3*tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_ind_poly,     tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_ind_poly_aux, tnumber*sizeof(int)) );
	gpuErrchk( cudaMalloc((void**) &cu_trivertex,    pnumber*sizeof(int)) );

	//std::cout << "Solicitada memoria cuda"<<std::endl;

		
    // Transfer arrays to device.
	
    gpuErrchk( cudaMemcpy(cu_r, r,                   2*pnumber*sizeof(double), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_triangles, triangles,   3*tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_adj, adj,               3*tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_seed, seed,    		 tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_max, max,               tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_mesh, mesh,             3*tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_ind_poly, ind_poly,     tnumber*sizeof(int), cudaMemcpyHostToDevice) );
	gpuErrchk( cudaMemcpy(cu_trivertex, trivertex,   pnumber*sizeof(int), cudaMemcpyHostToDevice) );
	
	//std::cout<<"copiada memoria cuda"<<std::endl;

	//se consigue el indice de la malla i_mesh
	unsigned long long int i_mesh = 0;
	unsigned long long int *cu_i_mesh = 0;
	gpuErrchk( cudaMalloc((void**) &cu_i_mesh, sizeof(unsigned long long int)) );
	gpuErrchk( cudaMemcpy(cu_i_mesh, &i_mesh, sizeof(unsigned long long int), cudaMemcpyHostToDevice) );
	
	//se consigue el indice de ind_poly
	unsigned long long int i_ind_poly = 0;
	unsigned long long int *cu_i_ind_poly = 0;
	gpuErrchk( cudaMalloc((void**) &cu_i_ind_poly, sizeof(unsigned long long int)) );
	gpuErrchk( cudaMemcpy(cu_i_ind_poly, &i_ind_poly, sizeof(unsigned long long int), cudaMemcpyHostToDevice) );
	
	int is_there_bet = 1;
	int *cu_is_there_bet;
	gpuErrchk( cudaMalloc((void**) &cu_is_there_bet, sizeof(int)) );
	//cudaMemcpy(cu_is_there_bet, &is_there_bet, sizeof(int), cudaMemcpyHostToDevice);

	int enumber = 3*tnumber;

	//https://stackoverflow.com/questions/47822784/calculating-grid-and-block-dimensions-of-a-kernel
	//int numThreads = 128;  // max register per block is 65536, 65536/512
	//int numBlocks  = (int)tnumber/numThreads;
	//int numBlocks  = (tnumber + (numThreads-1))/numThreads;
	//int numBlocks_edge  = (enumber + (numThreads-1))/numThreads;
	//std::cout<<"Llamando con "<<numBlocks<<" bloques y "<<numThreads<<std::endl;
	int numBlocks  = 65535;
	int numBlocks_edge = 65535;
	int numThreads = 1024;
	//Inicializar seeds y trivertex
	initialize_memory<<<numBlocks, numThreads>>>(cu_seed, tnumber);
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto t1 = std::chrono::high_resolution_clock::now();

	auto tb_label =std::chrono::high_resolution_clock::now();	
	//Label phase
	//Etiquetar el más largo;
	//std::cout<<"Inicia label longest"<<std::endl;
	auto tb_label_max = std::chrono::high_resolution_clock::now();
	label_longest_edges<<<numBlocks, numThreads>>>(cu_max, cu_r, cu_triangles, tnumber);
	gpuErrchk( cudaDeviceSynchronize() );
	auto te_label_max = std::chrono::high_resolution_clock::now();
	

	//Encontrar un triangulo semilla asociado al arco terminal
	auto tb_label_seed = std::chrono::high_resolution_clock::now();
	//get_seeds<<<numBlocks, numThreads>>>(cu_max, cu_triangles, cu_adj, cu_seed, tnumber);
	//std::cout<<"inicia get seeds"<<std::endl;
	get_seeds<<<numBlocks_edge, numThreads>>>(cu_max, cu_triangles, cu_adj, cu_seed, enumber);
	gpuErrchk( cudaDeviceSynchronize() );
	auto te_label_seed = std::chrono::high_resolution_clock::now();

	auto tb_label_non_frontier = std::chrono::high_resolution_clock::now();
	//Etiquetar label frontier-edges
	//label_frontier_edges<<<numBlocks, numThreads>>>(cu_max, cu_disconnect, cu_triangles, cu_adj, tnumber);
	label_frontier_edges<<<numBlocks_edge, numThreads>>>(cu_max, cu_triangles, cu_adj, enumber);
	gpuErrchk( cudaDeviceSynchronize() );
	auto te_label_non_frontier = std::chrono::high_resolution_clock::now();
	//std::cout<<"terminado label frontier"<<std::endl;

	auto te_label =std::chrono::high_resolution_clock::now();


	//Se ordenan las semillas
	//cudaMemcpy(adj, cu_adj,3*tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	//for (i = 0; i < tnumber; i++)
	//	//std::cout<<adj[3*i+0]<<" "<<adj[3*i+1]<<" "<<adj[3*i+2]<<"\n";

	//cudaMemcpy(seed, cu_seed,tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	//int num_region = 0;
	//for (i = 0; i < tnumber; i++)
	//{	
	//	if(seed[i] == TRUE){
	//		seed[num_region] = i;
	//		num_region++;
	//	}
	//}
	//for (i = 0; i < num_region; i++)
	//	//std::cout<<seed[i]<<" ";
	////std::cout<<"\nregiones = "<<num_region<<std::endl;


	auto tb_travel = std::chrono::high_resolution_clock::now();
	generate_mesh<<<numBlocks, numThreads>>>(cu_triangles, cu_adj, cu_r, cu_seed, cu_mesh, tnumber, cu_i_mesh, cu_ind_poly, cu_i_ind_poly);
	auto te_travel = std::chrono::high_resolution_clock::now();
	gpuErrchk( cudaDeviceSynchronize() );
	
	gpuErrchk( cudaMemcpy(&i_mesh, cu_i_mesh, sizeof(int), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(unsigned long long int), cudaMemcpyDeviceToHost) );

	
	gpuErrchk( cudaMemcpy(mesh, cu_mesh, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(ind_poly, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(seed, cu_seed, tnumber*sizeof(int), cudaMemcpyDeviceToHost));


	int num_poly;
	////std::cout<<"\n num poly: "<<i_ind_poly<<std::endl;
	//std::cout<<"Iniciando mesh reparation num poly: "<<i_ind_poly<<std::endl;
	int counter = 0;
	//cudaMemcpy(cu_ind_poly_aux, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToDevice);
	auto tb_reparation = std::chrono::high_resolution_clock::now();
	while(is_there_bet)
	{
		gpuErrchk( cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(unsigned long long int), cudaMemcpyDeviceToHost) );

		num_poly = i_ind_poly;
		
		i_mesh = 0;
		i_ind_poly = 0;
		is_there_bet = 0;

		gpuErrchk( cudaMemcpy(cu_i_mesh, &i_mesh, sizeof(unsigned long long int), cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(cu_i_ind_poly, &i_ind_poly, sizeof(unsigned long long int), cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(cu_is_there_bet, &is_there_bet, sizeof(int), cudaMemcpyHostToDevice) );

		if(counter%2 == 0){
			polygon_reparation2<<<numBlocks, numThreads>>>(cu_mesh, cu_mesh_aux, num_poly, cu_ind_poly, cu_ind_poly_aux, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			//polygon_reparation<<<numBlocks, numThreads>>>(cu_mesh, cu_mesh_aux, num_poly, cu_ind_poly, cu_ind_poly_aux, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			////std::cout<<"mesh esta en cu_mesh_aux"<<std::endl;
		}else{
			polygon_reparation2<<<numBlocks, numThreads>>>(cu_mesh_aux, cu_mesh, num_poly, cu_ind_poly_aux, cu_ind_poly, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			//polygon_reparation<<<numBlocks, numThreads>>>(cu_mesh_aux, cu_mesh, num_poly, cu_ind_poly_aux, cu_ind_poly, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			////std::cout<<"mesh esta en cu_mesh"<<std::endl;
		}
		
		gpuErrchk( cudaDeviceSynchronize() );

		counter++;
		gpuErrchk( cudaMemcpy(&is_there_bet, cu_is_there_bet, sizeof(int), cudaMemcpyDeviceToHost) );	
		////std::cout<<"has_bet? "<<is_there_bet<<", counter: "<<counter<<std::endl;	
	}
	
	auto te_reparation = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	if(counter%2 != 0){
	//	//std::cout<<"mesh esta en cu_mesh_aux"<<std::endl;
		gpuErrchk( cudaMemcpy(mesh, cu_mesh_aux, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost) );
		gpuErrchk( cudaMemcpy(ind_poly, cu_ind_poly_aux, tnumber*sizeof(int), cudaMemcpyDeviceToHost) );
	}else{
	//	//std::cout<<"mesh esta en cu_mesh"<<std::endl;
		gpuErrchk(cudaMemcpy(mesh, cu_mesh, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost) );
		gpuErrchk(cudaMemcpy(ind_poly, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToHost) );
	}

	gpuErrchk( cudaMemcpy(&i_mesh, cu_i_mesh, sizeof(unsigned long long int), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(unsigned long long int), cudaMemcpyDeviceToHost) );

	write_geomview(output, r, triangles, pnumber, tnumber, i_mesh, mesh, seed, i_ind_poly, 0);
	std::cout<<"File output in "<<output + ".off"<<std::endl;
	std::cout << std::setprecision(3) << std::fixed;
    std::cout <<"pnumber tnumber num_reg talgorithm tlabel tlabel_max tlabel_seed tlabel_non_frontier ttravel ttreparation"<<std::endl;
	std::cout<<pnumber<<" "<<tnumber<<" "<<i_ind_poly;
	uint total =std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count();
	uint t_label = std::chrono::duration_cast<std::chrono::milliseconds>(te_label - tb_label).count();
	uint t_label_max = std::chrono::duration_cast<std::chrono::milliseconds>(te_label_max - tb_label_max).count();
	uint t_label_seed = std::chrono::duration_cast<std::chrono::milliseconds>(te_label_seed - tb_label_seed).count();
	uint t_label_non_frontier = std::chrono::duration_cast<std::chrono::milliseconds>(te_label_non_frontier - tb_label_non_frontier).count();
	uint t_travel = std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel ).count();
	uint t_reparation = std::chrono::duration_cast<std::chrono::milliseconds>(te_reparation - tb_reparation ).count();
	std::cout<<" "<<t_label + t_travel + t_reparation;
	std::cout<<" "<<t_label;
	std::cout<<" "<<t_label_max;
	std::cout<<" "<<t_label_seed;
	std::cout<<" "<<t_label_non_frontier;
	std::cout<<" "<<t_travel;
	std::cout<<" "<<t_reparation;
	std::cout<<std::endl;
/*
  	//imprimir polginos
	std::cout<<"\n num poly: "<<i_ind_poly<<", i_mesh: "<<i_mesh<<std::endl;

	int k;
	for(i = 0; i < i_ind_poly; i++){
		std::cout<<"("<<ind_poly[i]<<") "<<mesh[ind_poly[i]]<<": ";
		for(k = 0; k < mesh[ind_poly[i]]; k++){
			std::cout<< mesh[ind_poly[i] + 1 + k]<<" ";
		}
		std::cout<<std::endl;
	}
*/
	/*
    i = 0;
    while(i < i_mesh){
        length_poly = mesh[i];
        std::cout<<"("<<i<<") "<<length_poly<<": ";
		i++;
        for(j=0; j < length_poly;j++){
            std::cout<< mesh[i]<<" ";
			i++;
        }
        std::cout<<std::endl;
    }
	
	for(i = 0; i < i_ind_poly; i++)	
		std::cout<< ind_poly[i]<<" ";
	std::cout<<std::endl;
	
	for(i = 0; i < i_mesh; i++)	
		std::cout<< mesh[i]<<" ";
	std::cout<<std::endl;
	*/

	free(r);
	free(triangles);
	free(adj);
	free(seed );
	free(mesh);
	free(max);
	free(ind_poly);

	gpuErrchk( cudaFree(cu_r) );
	gpuErrchk( cudaFree(cu_triangles) );
	gpuErrchk( cudaFree(cu_adj) );
	gpuErrchk( cudaFree(cu_seed) );
	gpuErrchk( cudaFree(cu_mesh) );
	gpuErrchk( cudaFree(cu_max) );
	gpuErrchk( cudaFree(cu_i_mesh) );
	gpuErrchk( cudaFree(cu_ind_poly) );
	gpuErrchk( cudaFree(cu_mesh_aux) );
	gpuErrchk( cudaFree(cu_ind_poly_aux) );
	gpuErrchk( cudaFree(cu_trivertex) );
	return EXIT_SUCCESS;
}
    

