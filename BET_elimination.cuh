__device__ int search_triangle_by_vertex_with_FrontierEdge_from_trivertex(int v, int *triangles, int *adj, int tnumber, int* trivertex);
__device__ int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber);
__device__ int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
__device__ int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
__device__ int advance_i_adjacents_triangles_share_endpoint(int adv, int t, int origen, int endpoint, int *p, int *adj);
__device__ int get_edge(int i, int u, int v, int *p);
__device__ int has_bet(int *poly, int length_poly);
__device__ int split_poly(int * poly, int length_poly, int * triangles, int * adj, double *r,  int &pos2_poly, int tnumber);


__global__ void polygon_reparation(int* cu_mesh, int* cu_mesh_aux, int num_poly, int *cu_ind_poly, int *cu_ind_poly_aux, int *cu_triangles, int tnumber, int *cu_adj, double *cu_r, int *cu_i_mesh, int* cu_i_ind_poly, int *is_there_bet);
__global__ void polygon_reparation2(int* cu_mesh, int* cu_mesh_aux, int num_poly, int *cu_ind_poly, int *cu_ind_poly_aux, int *cu_triangles, int tnumber, int *cu_adj, double *cu_r, int *cu_i_mesh, int* cu_i_ind_poly, int* cu_trivertex, int *is_there_bet);