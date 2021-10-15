#include "consts.h"
#include "triang_edge.cuh"
#include "polygon.cuh"
#include <stdio.h>

//Search a triangle asociated to a barrier-edge that contains vertex v
//This use trivertex to find a triangle asociated to v and travel through  adjacents triangles until find one with frontie-edges
//Input:  index vertex v, array of triangles, array of neigh, number of triangles, array that asociated each triangle to a vertex
//output: index of triangle that have one frontier-edge that contains v
__device__ int search_triangle_by_vertex_with_FrontierEdge_from_trivertex(int v, int *triangles, int *adj, int tnumber, int* trivertex){
	int t = trivertex[v];
	int origen = -1;
	int j, aux;
	while (1)
	{
		for (j = 0; j < 3; j++){
			//If the triangle contains v and has 2 fronter-edges
			if(triangles[3*t +j] == v  && ( adj[3*t + ((j + 1)%3)] == -1 || adj[3*t + ((j + 2)%3)] == -1 ))
				return t;
		}
		//avanza al siguiente triangulo
		aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v, triangles, adj);
        origen = aux;

	}	
	return -1;
}


//Search a triangle asociated to a barrier-edge that contains vertex v
//This doesnt use trivertex to find a triangle asociated to v, instead of, search in the list of the triangles one that containts v as frontier-edge, so the cost is o(n) per search
//Input:  index vertex v, array of triangles, array of neigh, number of triangles, array that asociated each triangle to a vertex
//output: index of triangle that have one frontier-edge that contains v
__device__ int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber){
	int i,j;
	for (i = 0; i < tnumber; i++)
		for (j = 0; j < 3; j++){
			//If the triangle contains v and has 2 fronter-edges
			if(triangles[3*i +j] == v  && ( adj[3*i + ((j + 1)%3)] == -1 || adj[3*i + ((j + 2)%3)] == -1 )){
				//debug_print("v %d |t %d | Triangles %d %d %d | ADJ  %d %d %d\n", v, i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
				return i;
			}
		}
				
	//fprintf(stderr,"%s:%d:%s() No se encontro triangulo que contiene el vertice %d \n",__FILE__,  __LINE__, __func__, v);
    //exit(0);
	return -1;
}


//Given a triangle i, return the triangle adjacent to the triangle origen that containts the vertex v
__device__ int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj){
	int t0, t1,t2;
	int a0, a1, a2;

	t0 = triangles[3*i + 0];
	t1 = triangles[3*i + 1];
	t2 = triangles[3*i + 2];

	a0 = adj[3*i + 0];
	a1 = adj[3*i + 1];
	a2 = adj[3*i + 2];

	//debug_print("origen %d, actual %d  | Triangles %d %d %d | ADJ  %d %d %d\n", origen,i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2] );

	if(t1 == v && origen == a0)
			return t2;
	else  if(t2 == v && origen == a0)
			return t1;
	else  if(t0 == v && origen == a1)
			return t2;
	else  if(t2 == v && origen == a1)
			return t0;
	else  if(t0 == v && origen == a2)
			return t1;
	else  if(t1 == v && origen == a2)
			return t0;	
	
//	fprintf(stderr,"%s:%d:%s()No se pudo encontrar el vertice anterior para partición \n",__FILE__,  __LINE__, __func__);
    //sexit(0);
	return -1;
}


//Given a triangle i, return the triangle no adjacent to the triangle origen that containts the vertex v
__device__ int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj){
	int t0, t1,t2;
	int a0, a1, a2;

	t0 = triangles[3*i + 0];
	t1 = triangles[3*i + 1];
	t2 = triangles[3*i + 2];

	a0 = adj[3*i + 0];
	a1 = adj[3*i + 1];
	a2 = adj[3*i + 2];

	//debug_print("v %d origen %d, actual %d  | Triangles %d %d %d | ADJ  %d %d %d\n", v, origen,i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);

	if(a0 != NO_ADJ && t1 == v && origen != a0)
			return t2;
	else if(a0 != NO_ADJ && t2 == v && origen != a0)
			return t1;		
	else if(a1 != NO_ADJ && t0 == v && origen != a1)
			return t2;
	else if(a1 != NO_ADJ && t2 == v && origen != a1)
			return t0;
	else if(a2 != NO_ADJ && t0 == v && origen != a2)
			return t1;
	else if(a2 != NO_ADJ && t1 == v && origen != a2)
			return t0;

	/*caso particular en poligonos grandes, ya no hay más triangulos para avanzar */
	if(get_adjacent_triangle_share_endpoint(i, origen, v, triangles, adj) == -2){
		//debug_msg("No se encuentran más triangulos para avanzar\n");
		return -2;
	}

//	fprintf(stderr,"%s:%d:%s()No se pudo encontrar el vertice siguiente para partición \n",__FILE__,  __LINE__, __func__);
    //exit(0);
	return -1;
}

// advance i triangles arround vertex endpoint
__device__ int advance_i_adjacents_triangles_share_endpoint(int adv, int t, int &origen, int endpoint, int *p, int *adj){
	int aux;
	while(adv > 0){
	//	printf("%d %d\n", t, origen) ;
		aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, endpoint, p, adj);
        origen = aux;
		adv--;
	}
	//printf("%d %d\n", t, origen) ;
	return t;
}

__device__ int get_edge(int i, int u, int v, int *p){
	int j, ind1,ind2;
	for(j = 0; j < 3; j++){
		ind1 = 3*i + j;
		ind2 = 3*i + (j+1)%3;
		////debug_print("%d %d %d %d %d\n", ind1, ind2, ind3, (p[ind1] == u || p[ind2] == u), (p[ind1] == v || p[ind2] == v));
		if( (p[ind1] == u || p[ind2] == u) && (p[ind1] == v || p[ind2] == v))
			return (j+2)%3;
	}
	//fprintf(stderr, "ERROR get_edge: No se encontro el edge %d - %d del triangulo %d", u,v,i);
	//exit(0);
	return -1;

}

__device__ int optimice_middle_edge_no_memory(int *t_original, int v_be, int *triangles, int *adj){
    
    int aux, origen,adv;

    int t_incident;
    int t = *t_original;
    t_incident = t;
    adv = 0;
    origen = -1; 
    int t_next = -1, t_prev;
    while (1)
    {
     //   debug_print("%d %d %d\n", *t_original, t, origen) ;
        adv++;
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        if (t<0)
            break;
    }
    //debug_print("%d %d %d\n", *t_original, t, origen );
    //print_poly(t_incident, i);
    if(adv == 1){
        *t_original = t_incident;
        //return search_prev_vertex_to_split(t_incident[0], v_be, -1, triangles, adj);
        for(int j = 0; j < 3; j++){
            if(triangles[3*t_incident + j] == v_be)
                return triangles[3*t_incident + (j+1)%3];
        }
    }
    if(adv%2 == 0){ //if the triangles surrounding the BET are even 
        adv = adv/2 - 1;
        origen = -1;
        t_prev = advance_i_adjacents_triangles_share_endpoint(adv,t_incident, origen, v_be, triangles, adj);
        *t_original = t_prev;
        //Choose the edge in common with the two middle triangles   
        t_next = get_adjacent_triangle_share_endpoint(t_prev, origen, v_be, triangles, adj);
   //     debug_print("search_prev adv %d t_next %d v_be %d t_prev %d \n", adv, t_next, v_be, t_prev);    
        return search_prev_vertex_to_split(t_next, v_be, t_prev, triangles, adj);
    }else{   
        //if the triangles surrounding the BET are odd, edges are even
        //Choose any edge of the triangle in the middle; prov is choose due to this always exists
        adv = adv/2;
        t_next = advance_i_adjacents_triangles_share_endpoint(adv,t_incident, t_prev, v_be, triangles, adj);
        *t_original = t_next;
        //t_next = get_adjacent_triangle_share_endpoint(t_prev, -1, v_be, triangles, adj);
 //       debug_print("search_next adv %d t_next %d v_be %d t_prev %d \n", adv, t_next, v_be, t_prev);
        return search_next_vertex_to_split(t_next, v_be, t_prev, triangles, adj);
    }   
}

__device__ int has_bet(int *poly, int length_poly){
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            return 1;
    }
    return 0;
}

__device__ int get_vertex_BarrierEdge(int *poly, int length_poly){
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            return poly[(i+1) %length_poly];
    }
    //fprintf(stderr,"num_BE %d\n", count_BarrierEdges(poly, length_poly));
    //fprintf(stderr,"%s:%d:%s(): No se encontro vertice BarrierEdge\n",__FILE__,  __LINE__, __func__);
    //exit(0);
    return -1;
}



__device__ int generate_polygon_from_BET_removal(int i, int * poly, int * triangles, int * adj, double *r, int ind_poly){
	//    int ind_poly = 0;
		
		int initial_point = 0;
		int end_point = 0;
		
		int t0;
		int t1;	
		int t2;
		int ind0;
		int ind1;
		int ind2;
		int continuous;
		int k, j, aux;
		int origen;
	
		int num_FrontierEdges = count_FrontierEdges(i, adj);
		//debug_print("Generando polinomio con triangulo %d FE %d\n", i, num_FrontierEdges);
		/*si tiene 3 se agregan y se corta el ciclo*/
		if (num_FrontierEdges == 3) {
			//debug_print("T %d Tiene 3 Frontier edge, se guardan así\n", i);
			poly[ind_poly] = triangles[3 * i + 0];
			ind_poly++;
			poly[ind_poly] = triangles[3 * i + 1];
			ind_poly++;
			poly[ind_poly] = triangles[3 * i + 2];
			ind_poly++;
	
			//visited[i] = TRUE;
			return ind_poly;
		} else if(num_FrontierEdges == 2) {
			//debug_print("T %d Tiene 2 Frontier edge, es oreja, se usa como semilla para generar el poly\n", i);
			/*si tiene dos FE se agregan y se empieza el ciclo*/
			for(j = 0; j<3; j++){
				ind0 = 3*i + j;
				ind1 = 3*i + (j+1)%3;
				ind2 = 3*i + (j+2)%3;
				if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
					poly[ind_poly] = triangles[ind1];
					ind_poly++;
					poly[ind_poly] = triangles[ind2];
					ind_poly++;
	
					initial_point = triangles[ind1];
					end_point = triangles[ind0];  
				}
			}
		}else if (num_FrontierEdges == 1){
			//debug_print("T %d Tiene 1 Frontier edge,se usa como FE initial\n", i);
			/*si tiene dos FE se agregan y se empieza el ciclo*/
			for(j = 0; j<3; j++){
				if(adj[3*i + j] == NO_ADJ){
					poly[ind_poly] = triangles[3*i + (j+1)%3];
					ind_poly++;
					initial_point = triangles[3*i + (j+1)%3];
	
					end_point = triangles[3*i + (j+2)%3];  
				}
			}
		}else {
			end_point = triangles[3*i + 0];
			initial_point = triangles[3*i + 0];
		}
		
		
		/*se marca como visitado */
		//visited[i] = TRUE;
		num_FrontierEdges = 0;
		k = i;
		aux = k;
		k = get_adjacent_triangle_share_endpoint(k, k, end_point, triangles, adj); /* cambia el indice */
		continuous = is_continuous(k, end_point, triangles);
		origen = aux;
	//        //debug_print("k %d origen %d, conti %d\n", k, origen, continuous);
		////debug_print("T_inicial %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
		////debug_print("initial_point %d endpoint %d | T_sig %d\n", initial_point, end_point, k);
	
		int triangugulo_initial = i;
		while (initial_point != end_point || triangugulo_initial != k) {
	
			/*se marca el triangulo visto como visitado y se suma al area del poligono */
			//searchandremove(hashtable_seed[hashy(k)], k);
		  //  visited[k] = TRUE;
			t0 = adj[3 * k + 0];
			t1 = adj[3 * k + 1];
			t2 = adj[3 * k + 2];
	
			num_FrontierEdges = count_FrontierEdges(k, adj);
			//debug_print("FE %d | origen %d t %d | Triangles %d %d %d | ADJ  %d %d %d\n", num_FrontierEdges, origen, k, triangles[3*k + 0], triangles[3*k + 1], triangles[3*k + 2], adj[3*k + 0], adj[3*k + 1], adj[3*k + 2]);
			if(origen == -2)
				return -3;
			if (num_FrontierEdges == 2 && continuous != -1) {
				/* ///////////////////si tiene 2 frontier edge se agregan a poly //////////////////////////////////// */
	
				if (t0 == NO_ADJ && t1 == NO_ADJ) {
					/*si endpoint es continua a t0  y t0-t1 son fe*/
					if (continuous == 1) {
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
	
						end_point = triangles[3 * k + 0];
	
					} else if (continuous == 0) {
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
	
						end_point = triangles[3 * k + 1];
					}
				} else if (t2 == NO_ADJ && t0 == NO_ADJ) {
					/*si endpoint es continua a t2  y t2-t0 son fe*/
					if (continuous == 0) {
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
	
						end_point = triangles[3 * k + 2];
	
					} else if (continuous == 2) {
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
	
						end_point = triangles[3 * k + 0];
	
					}
				} else if (t1 == NO_ADJ && t2 == NO_ADJ) {
					/*si endpoint es continua a t1 y t1-t2 son fe*/
					if (continuous == 2) {
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
	
						end_point = triangles[3 * k + 1];
	
					} else if (continuous == 1) {
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
	
						end_point = triangles[3 * k + 2];
	
					}
				}
	
				aux = k;
				k = get_adjacent_triangle_share_endpoint(k, -1, end_point, triangles, adj); /* se le permite volver al triangulo anterior */
				continuous = is_continuous(k, end_point, triangles);
				origen = aux;
	
			} else if (num_FrontierEdges == 1 && continuous != -1) {
				/* ///////////////////si solo se tiene 1 frontier edge //////////////////////////////////// */
				if (t0 == NO_ADJ) {
					/*si endpoint es continua a t0  y t0 es fe*/
					if (continuous == 1) {
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
	
						end_point = triangles[3 * k + 2];
	
					} else if (continuous == 2) {
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
	
						end_point = triangles[3 * k + 1];
	
					}
				} else if (t2 == NO_ADJ) {
					/*si endpoint es continua a t2  y t2 es fe*/
					if (continuous == 0) {
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
	
						end_point = triangles[3 * k + 1];
	
					} else if (continuous == 1) {
						poly[ind_poly] = triangles[3 * k + 1];
						ind_poly++;
	
						end_point = triangles[3 * k + 0];
	
					}
				} else if (t1 == NO_ADJ) {
					/*si endpoint es continua a t1  y t1 es fe*/
					if (continuous == 2) {
						poly[ind_poly] = triangles[3 * k + 2];
						ind_poly++;
	
						end_point = triangles[3 * k + 0];
	
					} else if (continuous == 0) {
						poly[ind_poly] = triangles[3 * k + 0];
						ind_poly++;
	
						end_point = triangles[3 * k + 2];
	
					}
				}
				/*si es continuo y tiene 1 fe no puede volver, ind si se guarda  o no*/
				aux = k;
				k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
				origen = aux;
				continuous = is_continuous(k, end_point, triangles);
			} else {
				/*si no es continuo no puede regresar de donde venía*/
				aux = k;
				k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
				continuous = is_continuous(k, end_point, triangles);
				origen = aux;
			}
	
		}
		
		return ind_poly;
	}



//Split a polygon in two and save it in the same array
__device__ int split_poly(int * poly, int length_poly, int * triangles, int * adj, double *r, int* cu_trivertex, int &pos2_poly, int tnumber){

	int v_be, v_other, t1,t2, ipoly, ipoly_after;


	v_be = get_vertex_BarrierEdge(poly, length_poly);
	//t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
	t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, cu_trivertex);
	v_other = optimice_middle_edge_no_memory(&t1, v_be, triangles, adj);
	t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);

	adj[3*t1 + get_shared_edge(t1, v_be, v_other, triangles)] = NO_ADJ;
	adj[3*t2 + get_shared_edge(t2, v_be, v_other, triangles)] = NO_ADJ;

	ipoly = 0;

	ipoly_after = generate_polygon_from_BET_removal(t1, poly,triangles,adj,r,ipoly+1); 
	poly[ipoly] = ipoly_after - ipoly - 1; // calculate lenght poly and save it before their vertex
	ipoly = ipoly_after;

	pos2_poly = ipoly;
	
	ipoly_after = generate_polygon_from_BET_removal(t2, poly,triangles,adj,r,ipoly+1); 
	poly[ipoly] = ipoly_after - ipoly - 1; // calculate lenght poly and save it before their vertex
	//ipoly = ipoly_after + ipoly ; //lenght of element in poly;

	return ipoly_after;
}

//Split polygons with BET in two polygons and save again in cu_mesh, if a polygon hasn't bet, this is relloc in 
//INPUT
//cu_mesh: polygon mesh with bet
//num_poly: Número de poligonos en ind_poly
//ind_poly: index first element of the poly in cu_mesh
//cu_triangles, cu_adj, cu_r: Elements of Delaunay triangulation
//range_mesh: atomic variable to save elements in mesh
//range_ind_poly: atomic variable to save elements in ind_poly
//is_there_bet: atomic variable to check if there are bet or not
//OUTPUT
//cu_mesh_aux: new polygon mesh
//is_there_bet: atomic variable to check if there are bet or not
__global__ void polygon_reparation(int* cu_mesh, int* cu_mesh_aux, int num_poly, int *cu_ind_poly, int *cu_ind_poly_aux, int *cu_triangles, int tnumber, int *cu_adj, double *cu_r, unsigned long long int *cu_i_mesh, unsigned long long int* cu_i_ind_poly, int* cu_trivertex, int *is_there_bet){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    //int i_mesh, i_mesh2, bet, length_poly, poly[100], k, pos2_poly;
	int bet = 0, length_poly, poly[100], k, pos2_poly, i_mesh, i_ind_poly;	
	
	if(i < num_poly){

		//save polygon in temporal array
		
		i_ind_poly = cu_ind_poly[i];
		length_poly = cu_mesh[i_ind_poly];
		
		for(k = 0; k < length_poly; k++){
			poly[k] = cu_mesh[i_ind_poly + 1 + k];
		
		}
		bet = has_bet(poly, length_poly);
		//bet = 0;
		if(bet){//if has bet
			//printf("Poly to split %d: ", length_poly);
			for(k = 0; k < length_poly; k++){
				poly[k] = cu_mesh[i_ind_poly + 1 + k];
			//	printf("%d ", poly[k]);
			}
			//printf("\n");

			length_poly = split_poly(poly, length_poly, cu_triangles, cu_adj, cu_r, cu_trivertex, pos2_poly, tnumber);
			i_mesh = atomicAdd(cu_i_mesh, length_poly); //imesh indice inicial a guardar
			i_ind_poly = atomicAdd(cu_i_ind_poly, 2);
			for(k = 0; k < length_poly; k++)
				cu_mesh_aux[i_mesh + k] = poly[k];
			
			//printf("resultado %d: ", length_poly);
			//for(k = 0; k < length_poly; k++){
			//	printf("%d ", poly[k]);
			//}
			//printf("\n");
			cu_ind_poly_aux[i_ind_poly] = i_mesh;
			cu_ind_poly_aux[i_ind_poly+1] = i_mesh + poly[0] + 1;
			//printf("i_mesh_1: %d - i_ind_poly1: %d | i_mesh2: %d - i_ind_poly2: %d\n",i_mesh, i_ind_poly - 1, i_mesh + poly[0] + 1, i_ind_poly);
		}else{// if no bet, then just preprare the polygon to save in array
			for(k = length_poly; k > 0 ; k--)
				poly[k] = poly[k-1];
			poly[0] = length_poly;

			i_mesh = atomicAdd(cu_i_mesh, length_poly+1); //imesh indice inicial a guardar
			i_ind_poly = atomicAdd(cu_i_ind_poly, 1);

			for(int k = 0; k < length_poly + 1; k++)
				cu_mesh_aux[i_mesh + k] = poly[k];
			cu_ind_poly_aux[i_ind_poly] = i_mesh;
		}

		atomicAdd(is_there_bet, bet);
	}	
}

//Split polygons with BET in two polygons and save again in cu_mesh, if a polygon hasn't bet, this is relloc in 
//INPUT
//cu_mesh: polygon mesh with bet
//num_poly: Número de poligonos en ind_poly
//ind_poly: index first element of the poly in cu_mesh
//cu_triangles, cu_adj, cu_r: Elements of Delaunay triangulation
//range_mesh: atomic variable to save elements in mesh
//range_ind_poly: atomic variable to save elements in ind_poly
//is_there_bet: atomic variable to check if there are bet or not
//OUTPUT
//cu_mesh_aux: new polygon mesh
//is_there_bet: atomic variable to check if there are bet or not
__global__ void polygon_reparation2(int* cu_mesh, int* cu_mesh_aux, int num_poly, int *cu_ind_poly, int *cu_ind_poly_aux, int *cu_triangles, int tnumber, int *cu_adj, double *cu_r, unsigned long long int *cu_i_mesh, unsigned long long int* cu_i_ind_poly, int* cu_trivertex, int *is_there_bet){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    //int i_mesh, i_mesh2, bet, length_poly, poly[100], k, pos2_poly;
	int bet = 0, length_poly, k, i_mesh, i_ind_poly, i_ind_poly_aux;	
	int v_be, v_other, t1,t2, ipoly, ipoly_after;
	int x,y;
	v_be = -1;
	if(i < num_poly){
		//save polygon in temporal array
		i_ind_poly = cu_ind_poly[i];
		length_poly = cu_mesh[i_ind_poly];
		//printf("poly %d lenght %d\n", i_ind_poly, length_poly);
		for (k = 0; k < length_poly; k++)
		{
			x = k % length_poly;
			y = (k+2) % length_poly;
			if (cu_mesh[i_ind_poly + 1 + x] == cu_mesh[i_ind_poly + 1 + y]){
				v_be = cu_mesh[i_ind_poly+1 + (k+1)%length_poly];
				bet=1;
				break;
			}
		}

		//bet = 0;
		if(bet){//if has bet

			i_mesh = atomicAdd(cu_i_mesh, length_poly + 1 + 2+ 1); // +1 lenght poly, +2 new size of sum two poly, + 1 leng_poly2
			i_ind_poly = atomicAdd(cu_i_ind_poly, 2);
			
			
			//v_be = get_vertex_BarrierEdge(poly, length_poly);
			//t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, cu_triangles, cu_adj, tnumber, cu_trivertex);
			t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, cu_triangles, cu_adj, tnumber);
			v_other = optimice_middle_edge_no_memory(&t1, v_be, cu_triangles, cu_adj);
			t2 = get_adjacent_triangle(t1, v_other, v_be, cu_triangles, cu_adj);

			cu_adj[3*t1 + get_shared_edge(t1, v_be, v_other, cu_triangles)] = NO_ADJ;
			cu_adj[3*t2 + get_shared_edge(t2, v_be, v_other, cu_triangles)] = NO_ADJ;

			ipoly = i_mesh;

			ipoly_after = generate_polygon_from_BET_removal(t1, cu_mesh_aux,cu_triangles, cu_adj,cu_r,ipoly+1); 
			cu_mesh_aux[ipoly] = ipoly_after - ipoly - 1; // calculate lenght poly and save it before their vertex
			ipoly = ipoly_after;
			
			ipoly_after = generate_polygon_from_BET_removal(t2, cu_mesh_aux, cu_triangles,cu_adj,cu_r,ipoly+1); 
			cu_mesh_aux[ipoly] = ipoly_after - ipoly - 1; // calculate lenght poly and save it before their vertex

			cu_ind_poly_aux[i_ind_poly] = i_mesh;
			cu_ind_poly_aux[i_ind_poly+1] = i_mesh + cu_mesh_aux[i_mesh] + 1;
			
		}else{// if no bet, then just preprare the polygon to save in array
			i_mesh = atomicAdd(cu_i_mesh, length_poly+1); //imesh indice inicial a guardar
			i_ind_poly_aux = atomicAdd(cu_i_ind_poly, 1);
			for(k = 0; k < length_poly + 1; k++)
				cu_mesh_aux[i_mesh + k] = cu_mesh[i_ind_poly + k];
			cu_ind_poly_aux[i_ind_poly_aux] = i_mesh;
		}

		atomicAdd(is_there_bet, bet);
	}	
}


/*
//Split polygons with BET in two polygons and save again in cu_mesh, if a polygon hasn't bet, this is relloc in 
//INPUT
//cu_mesh: polygon mesh with bet
//num_poly: Número de poligonos en ind_poly
//ind_poly: index first element of the poly in cu_mesh
//cu_triangles, cu_adj, cu_r: Elements of Delaunay triangulation
//range_mesh: atomic variable to save elements in mesh
//range_ind_poly: atomic variable to save elements in ind_poly
//is_there_bet: atomic variable to check if there are bet or not
//OUTPUT
//cu_mesh_aux: new polygon mesh
//is_there_bet: atomic variable to check if there are bet or not
__global__ void polygon_reparation3(int* cu_mesh, int* cu_mesh_aux, int num_poly, int *cu_ind_poly, int *cu_ind_poly_aux, int* cu_seed, int *cu_triangles, int tnumber, int *cu_adj, double *cu_r, int *cu_i_mesh, int* cu_i_ind_poly, int *is_there_bet){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    //int i_mesh, i_mesh2, bet, length_poly, poly[100], k, pos2_poly;
	int bet = 0, length_poly, poly[100], k, pos2_poly, i_mesh, i_ind_poly;	
	int seeds[20], s = 0;
	int t1, v_other, t2, v_be;
	if(i < num_poly){

		//save polygon in temporal array
		i_ind_poly = cu_ind_poly[i];
		length_poly = cu_mesh[i_ind_poly];
		
		for(k = 0; k < length_poly; k++){
			poly[k] = cu_mesh[i_ind_poly + 1 + k];
		}

		int x, y;
		for (k = 0; k < length_poly; k++)
		{
			x = k % length_poly;
			y = (k+2) % length_poly;
			if (poly[x] == poly[y]){

				v_be = poly[(k+1) %length_poly];
				t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, cu_triangles, cu_adj, tnumber);
				v_other = optimice_middle_edge_no_memory(&t1, v_be, cu_triangles, cu_adj);
				t2 = get_adjacent_triangle(t1, v_other, v_be, cu_triangles, cu_adj);

				cu_adj[3*t1 + get_shared_edge(t1, v_be, v_other, cu_triangles)] = NO_ADJ;
				cu_adj[3*t2 + get_shared_edge(t2, v_be, v_other, cu_triangles)] = NO_ADJ;	

				seeds[]

			}
				
		}

		if(bet){//if has bet
			for(k = 0; k < length_poly; k++){
				poly[k] = cu_mesh[i_ind_poly + 1 + k];
			}
			length_poly = split_poly(poly, length_poly, cu_triangles, cu_adj, cu_r, pos2_poly, tnumber);
			i_mesh = atomicAdd(cu_i_mesh, length_poly); //imesh indice inicial a guardar
			i_ind_poly = atomicAdd(cu_i_ind_poly, 2);
			for(k = 0; k < length_poly; k++)
				cu_mesh_aux[i_mesh + k] = poly[k];
			
			cu_ind_poly_aux[i_ind_poly] = i_mesh;
			cu_ind_poly_aux[i_ind_poly+1] = i_mesh + poly[0] + 1;
			
		}else{// if no bet, then just preprare the polygon to save in array
			for(k = length_poly; k > 0 ; k--)
				poly[k] = poly[k-1];
			poly[0] = length_poly;

			i_mesh = atomicAdd(cu_i_mesh, length_poly+1); //imesh indice inicial a guardar
			i_ind_poly = atomicAdd(cu_i_ind_poly, 1);

			for(int k = 0; k < length_poly + 1; k++)
				cu_mesh_aux[i_mesh + k] = poly[k];
			cu_ind_poly_aux[i_ind_poly] = i_mesh;
		}

		atomicAdd(is_there_bet, bet);
	}	
}
*/