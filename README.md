# CUDA Terminal-edge region mesh generation

<p align="center">
 <img src="https://github.com/ssalinasfe/CUDA-Terminal-edge-mesh/blob/main/images/cuycuy.png" width="80%">
</p>

New algorithm to generate polygonal meshes of arbitrary shape based on a terminal-edge regions.

The algorithm needs a triangulation as input, this triangulation is represent using 3 files:

- **[.node](https://www.cs.cmu.edu/~quake/triangle.node.html)**: Point file
- **[.ele](https://www.cs.cmu.edu/~quake/triangle.ele.html)**: Triangle file
- **[.neigh](https://www.cs.cmu.edu/~quake/triangle.neigh.html)**: adjacency list

The argv of the algorithm are

```./PolyllaCUDA  <name file without extention> <output filename> ```

Example of a triangulation to use is in example folder

```./PolyllaCUDA  square_3000000.1 off_file ```



## Other versions

- Original sequential: https://github.com/ssalinasfe/Polylla-Mesh
- Sequential with POO: https://github.com/ssalinasfe/Polylla-Mesh-DCEL

