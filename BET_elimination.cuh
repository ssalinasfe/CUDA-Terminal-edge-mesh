__device__ int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber);
__device__ int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
__device__ int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
__device__ int advance_i_adjacents_triangles_share_endpoint(int adv, int t, int origen, int endpoint, int *p, int *adj);
__device__ int get_edge(int i, int u, int v, int *p);