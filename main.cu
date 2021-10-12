/*
- No se pude implementar el remove BE secuencial porqué la cantidad de poligonos a guardar en la nueva malla, no hay forma de usar guardar n nuevos poligonos y la vez los triangulos semilla de estos en poligonos. Guardar ambos en registos sobrepasa la capacidad de registros por hilo
*/

#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include "detri2.h"
#include "polymesh.h"
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



int main(int argc, char* argv[])
{

    int nparam = 3;
    //char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("test.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506randompoints.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506equilateral.node")};
    int print_triangles = 0;
    char* ppath;
    //char* ppath = const_cast<char*> ("test");
    //TMesh *Tr = new TMesh(nparam, params);    
	//auto tb_delaunay = std::chrono::high_resolution_clock::now();
	//TMesh *Tr = new TMesh(argc, argv);    	
	//auto te_delaunay = std::chrono::high_resolution_clock::now();
    //Tr->print();
    
	
	int tnumber, pnumber, i;
	double *r;
	int *triangles;
	int *adj;
    int *seed;
	int *max;
	int *mesh;
	int *disconnect;
	int *ind_poly;

	std::string name(argv[1]);
	std::cout<<name<<std::endl;
	read_from_triangle(name, pnumber, tnumber, r, triangles, adj);
	std::cout << " " << tnumber << " " << pnumber << "\n";

    //tnumber = Tr->tnumber;
    //pnumber = Tr->pnumber;
    //r = (double *)malloc(2*tnumber*sizeof(double));
    //adj =(int *)malloc(3*tnumber*sizeof(int));
    //triangles = (int *)malloc(3*tnumber*sizeof(int));
	max = (int *)malloc(tnumber*sizeof(int));
	disconnect = (int *)malloc(3*tnumber*sizeof(int));
	seed = (int *)malloc(tnumber*sizeof(int));
	mesh = (int *)malloc(3*tnumber*sizeof(int));
	ind_poly = (int *)malloc(3*tnumber*sizeof(int));

	//Cuda functions
    // Initialize device pointers.
    double *cu_r;
	int *cu_triangles;
	int *cu_adj;
    int *cu_seed;
	int *cu_max;
	int *cu_disconnect;
	int *cu_mesh;
	int *cu_mesh_aux;
	int *cu_ind_poly;
	int *cu_ind_poly_aux;
	int *cu_trivertex;

	// Allocate device memory.
	cudaMalloc((void**) &cu_max, tnumber*sizeof(int));
	cudaMalloc((void**) &cu_seed, tnumber*sizeof(int));
	cudaMalloc((void**) &cu_disconnect, 3*tnumber*sizeof(int));
	cudaMalloc((void**) &cu_r, 2*tnumber*sizeof(double));
	cudaMalloc((void**) &cu_triangles, 3*tnumber*sizeof(int));
	cudaMalloc((void**) &cu_adj, 3*tnumber*sizeof(int));
	cudaMalloc((void**) &cu_mesh, 3*tnumber*sizeof(int));
	cudaMalloc((void**) &cu_mesh_aux, 3*tnumber*sizeof(int));
	cudaMalloc((void**) &cu_ind_poly, tnumber*sizeof(int));
	cudaMalloc((void**) &cu_ind_poly_aux, tnumber*sizeof(int));
	cudaMalloc((void**) &cu_trivertex, pnumber*sizeof(int));

	/* Llamada a detr2 */
	{
	/*
    int idx =0;
    //copiar arreglo de vertices
    //std::cout<<"pnumber "<<pnumber<<std::endl;
    for (i = 0; i < Tr->trimesh->ct_in_vrts; i++) {
        if (!Tr->trimesh->io_keep_unused) { // no -IJ
            if (Tr->trimesh->in_vrts[i].typ == UNUSEDVERTEX) continue;
        }
        r[2*i + 0]= Tr->trimesh->in_vrts[i].crd[0];
        r[2*i + 1]= Tr->trimesh->in_vrts[i].crd[1];
        //std::cout<<idx<<" ("<<r[2*i + 0]<<", "<<r[2*i + 1]<<") "<<std::endl;
        Tr->trimesh->in_vrts[i].idx = idx;
        idx++;
    }
    idx = 0;
    for (int i = 0; i < Tr->trimesh->tr_tris->used_items; i++) {
        detri2::Triang* tri = (detri2::Triang *) Tr->trimesh->tr_tris->get(i);
        if (tri->is_deleted()) continue;
        if (tri->is_hulltri()) {
            tri->idx = -1;
        } else {
            tri->idx = idx;
            idx++;
        }
    }

    //std::cout<<"tnumber: "<<Tr->trimesh->tr_tris->objects - Tr->trimesh->ct_hullsize<<std::endl;
    idx = 0;
    for (int i = 0; i < Tr->trimesh->tr_tris->used_items; i++)
    {
        
        detri2::Triang* tri = (detri2::Triang *) Tr->trimesh->tr_tris->get(i);
        if (tri->is_deleted() || tri->is_hulltri()) continue;
        triangles[3*idx+0] = tri->vrt[0]->idx;
        triangles[3*idx+1] = tri->vrt[1]->idx;
        triangles[3*idx+2] = tri->vrt[2]->idx;
        adj[3*idx+ 0] = tri->nei[0].tri->idx;
        adj[3*idx+ 1] = tri->nei[1].tri->idx;
        adj[3*idx+ 2] = tri->nei[2].tri->idx;
        //std::cout<<idx<<" | "<<triangles[3*idx+0]<<" "<<triangles[3*idx+1]<<" "<<triangles[3*idx+2]<<" | ";
        //std::cout<<adj[3*idx+ 0]<<" "<<adj[3*idx+ 1]<<" "<<adj[3*idx+ 2]<<" | "<<std::endl;
        idx++;
    }
	delete Tr;
	*/
	}

		
    // Transfer arrays to device.
    cudaMemcpy(cu_r, r,                   2*tnumber*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_triangles, triangles,   3*tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_adj, adj,               3*tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_seed, seed,    		  tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_max, max,               tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_disconnect, disconnect, 3*tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_mesh, mesh,             3*tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_ind_poly, ind_poly,    tnumber*sizeof(int), cudaMemcpyHostToDevice);
	
	//se consigue el indice de la malla i_mesh
	int i_mesh = 0;
	int *cu_i_mesh;
	cudaMalloc((void**) &cu_i_mesh, sizeof(int));
	cudaMemcpy(cu_i_mesh, &i_mesh, 1*sizeof(int), cudaMemcpyHostToDevice);
	
	//se consigue el indice de ind_poly
	int i_ind_poly = 0;
	int *cu_i_ind_poly;
	cudaMalloc((void**) &cu_i_ind_poly, sizeof(int));
	cudaMemcpy(cu_i_ind_poly, &i_ind_poly, 1*sizeof(int), cudaMemcpyHostToDevice);
	
	int is_there_bet = 1;
	int *cu_is_there_bet;
	cudaMalloc((void**) &cu_is_there_bet, sizeof(int));
	//cudaMemcpy(cu_is_there_bet, &is_there_bet, 1*sizeof(int), cudaMemcpyHostToDevice);

	int enumber = 3*tnumber;

	//https://stackoverflow.com/questions/47822784/calculating-grid-and-block-dimensions-of-a-kernel
	int numThreads = 128;  // max register per block is 65536, 65536/512
	//int numBlocks  = (int)tnumber/numThreads;
	int numBlocks  = (tnumber + (numThreads-1))/numThreads;
	int numBlocks_edge  = (enumber + (numThreads-1))/numThreads;


	//Inicializar seeds y trivertex
	initialize_memory<<<numBlocks, numThreads>>>(cu_seed, cu_trivertex, cu_triangles, tnumber);
	cudaDeviceSynchronize();
	
	auto t1 = std::chrono::high_resolution_clock::now();

	auto tb_label =std::chrono::high_resolution_clock::now();	
	//Label phase
	//Etiquetar el más largo;
	std::cout<<"Inicia label longest"<<std::endl;
	auto tb_label_max = std::chrono::high_resolution_clock::now();
	label_longest_edges<<<numBlocks, numThreads>>>(cu_max, cu_r, cu_triangles, tnumber);
	cudaDeviceSynchronize();
	auto te_label_max = std::chrono::high_resolution_clock::now();
	

	//Encontrar un triangulo semilla asociado al arco terminal
	auto tb_label_seed = std::chrono::high_resolution_clock::now();
	//get_seeds<<<numBlocks, numThreads>>>(cu_max, cu_triangles, cu_adj, cu_seed, tnumber);
	std::cout<<"inicia get seeds"<<std::endl;
	get_seeds<<<numBlocks_edge, numThreads>>>(cu_max, cu_triangles, cu_adj, cu_seed, enumber);
	cudaDeviceSynchronize();
	auto te_label_seed = std::chrono::high_resolution_clock::now();

	auto tb_label_non_frontier = std::chrono::high_resolution_clock::now();
	//Etiquetar label frontier-edges
	//label_frontier_edges<<<numBlocks, numThreads>>>(cu_max, cu_disconnect, cu_triangles, cu_adj, tnumber);
	label_frontier_edges<<<numBlocks_edge, numThreads>>>(cu_max, cu_disconnect, cu_triangles, cu_adj, enumber);
	cudaDeviceSynchronize();
	auto te_label_non_frontier = std::chrono::high_resolution_clock::now();
	std::cout<<"terminado label frontier"<<std::endl;
	//Desconectar frontier-edges
	//disconnect_edges<<<numBlocks, numThreads>>>(cu_adj, cu_disconnect, tnumber);
	//disconnect_edges<<<numBlocks_edge, numThreads>>>(cu_adj, cu_disconnect, enumber);
	cudaDeviceSynchronize();
	std::cout<<"terminado disconnect"<<std::endl;

	auto te_label =std::chrono::high_resolution_clock::now();


	//Se ordenan las semillas
	//cudaMemcpy(adj, cu_adj,3*tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	//for (i = 0; i < tnumber; i++)
	//	std::cout<<adj[3*i+0]<<" "<<adj[3*i+1]<<" "<<adj[3*i+2]<<"\n";

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
	//	std::cout<<seed[i]<<" ";
	//std::cout<<"\nregiones = "<<num_region<<std::endl;


	auto tb_travel = std::chrono::high_resolution_clock::now();
	generate_mesh<<<numBlocks, numThreads>>>(cu_triangles, cu_adj, cu_r, cu_seed, cu_mesh, tnumber, cu_i_mesh, cu_ind_poly, cu_i_ind_poly);
	std::cout<<"terminado mesh generation"<<std::endl;
	cudaDeviceSynchronize();
	auto te_travel = std::chrono::high_resolution_clock::now();
	
	
	cudaMemcpy(&i_mesh, cu_i_mesh, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(int), cudaMemcpyDeviceToHost);

	
	cudaMemcpy(mesh, cu_mesh, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(ind_poly, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(seed, cu_seed, tnumber*sizeof(int), cudaMemcpyDeviceToHost);


	int num_poly;
	//std::cout<<"\n num poly: "<<i_ind_poly<<std::endl;
	int counter = 0;
	//cudaMemcpy(cu_ind_poly_aux, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToDevice);
	auto tb_reparation = std::chrono::high_resolution_clock::now();
	while(is_there_bet)
	{
		cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(int), cudaMemcpyDeviceToHost);

		num_poly = i_ind_poly;
		
		i_mesh = 0;
		i_ind_poly = 0;
		is_there_bet = 0;

		cudaMemcpy(cu_i_mesh, &i_mesh, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(cu_i_ind_poly, &i_ind_poly, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(cu_is_there_bet, &is_there_bet, sizeof(int), cudaMemcpyHostToDevice);

		if(counter%2 == 0){
			polygon_reparation2<<<numBlocks, numThreads>>>(cu_mesh, cu_mesh_aux, num_poly, cu_ind_poly, cu_ind_poly_aux, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			//std::cout<<"mesh esta en cu_mesh_aux"<<std::endl;
		}else{
			polygon_reparation2<<<numBlocks, numThreads>>>(cu_mesh_aux, cu_mesh, num_poly, cu_ind_poly_aux, cu_ind_poly, cu_triangles, tnumber, cu_adj, cu_r, cu_i_mesh, cu_i_ind_poly, cu_trivertex, cu_is_there_bet);
			//std::cout<<"mesh esta en cu_mesh"<<std::endl;
		}
		
		cudaDeviceSynchronize();

		counter++;
		cudaMemcpy(&is_there_bet, cu_is_there_bet, 1*sizeof(int), cudaMemcpyDeviceToHost);	
		//std::cout<<"has_bet? "<<is_there_bet<<", counter: "<<counter<<std::endl;	
	}
	
	auto te_reparation = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	if(counter%2 != 0){
	//	std::cout<<"mesh esta en cu_mesh_aux"<<std::endl;
		cudaMemcpy(mesh, cu_mesh_aux, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(ind_poly, cu_ind_poly_aux, tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	}else{
	//	std::cout<<"mesh esta en cu_mesh"<<std::endl;
		cudaMemcpy(mesh, cu_mesh, 3*tnumber*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(ind_poly, cu_ind_poly, tnumber*sizeof(int), cudaMemcpyDeviceToHost);
	}

	cudaMemcpy(&i_mesh, cu_i_mesh, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&i_ind_poly, cu_i_ind_poly, sizeof(int), cudaMemcpyDeviceToHost);

	write_geomview(name, r, triangles, pnumber, tnumber, i_mesh, mesh, seed, i_ind_poly, 0);

	std::cout << std::setprecision(3) << std::fixed;
    std::cout <<"pnumber tnumber num_reg talgorithm tlabel tlabel_max tlabel_seed tlabel_non_frontier ttravel ttreparation"<<std::endl;
	std::cout<<pnumber<<" "<<tnumber<<" "<<i_ind_poly;
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label - tb_label).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label_max - tb_label_max).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label_seed - tb_label_seed).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_label_non_frontier - tb_label_non_frontier).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel ).count();
	std::cout<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(te_reparation - tb_reparation ).count();

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

	cudaFree(cu_r);
	cudaFree(cu_triangles);
	cudaFree(cu_adj);
	cudaFree(cu_seed);
	cudaFree(cu_mesh);
	cudaFree(cu_max);
	cudaFree(cu_i_mesh);
	cudaFree(cu_disconnect);
	cudaFree(cu_ind_poly);
	cudaFree(cu_mesh_aux);
	cudaFree(cu_ind_poly_aux);
	return EXIT_SUCCESS;
}
    

