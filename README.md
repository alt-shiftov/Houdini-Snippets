# VEX-Snippets

## Basic 3D Graphics Functions (Intersection / Nearest Distance / Projection) 

### [Intersection Functions](VEX/Intersection)
I was trying to explore how intersections work in 3D graphics, so I made a collection of basic functions designed to find intersections in 3D space.
These functions cover various intersection scenarios, including:

- [`Line-Line`](VEX/Intersection/intersection_functions.h#L5)
- [`Segment-Segment`](VEX/Intersection/intersection_functions.h#L66)
- [`Segment-Curve`](VEX/Intersection/intersection_functions.h#L152)
- [`Line-Plane`](VEX/Intersection/intersection_functions.h#L210)
- [`Plane-Plane`](VEX/Intersection/intersection_functions.h#L291)
- [`Segment-Triangle`](VEX/Intersection/intersection_functions.h#L337)
- [`Triangle-Triangle`](VEX/Intersection/intersection_functions.h#L406)
- [`Segment-Planar Quad`](VEX/Intersection/intersection_functions.h#L485)
- [`Segment-Bilinear Patch (both planar and non-planar quad)`](VEX/Intersection/intersection_functions.h#L577)

Based on this research I developed the <img src="readme_images/icons/booleanfracture.svg" width=15px> ![Shatterizer](https://github.com/alt-shiftov/Shatterizer) which improves precision of boolean fracturing as well as the <img src="readme_images/icons/intersectionanalysis.svg" width=15px> ![Intersection Analysis Custom](Nodes/sop_ALT.intersectionanalysis_custom.1.0.hdalc) node.
<br><br> <img src="readme_images/intersection_vis.gif" width=450px> <br>


### [Nearest Distance (Geometric Tools)](VEX/NearestDistance)
This collection of functions is designed to find the nearest distances between various geometric entities. <br>
Based on functions from [Geometric Tools](https://github.com/davideberly/GeometricTools/):

- [`dist_point_segment`](VEX/NearestDistance/neardistance.h#L11)
- [`dist_segment_segment`](VEX/NearestDistance/neardistance.h#L62)
- [`dist_line_line`](VEX/NearestDistance/neardistance.h#L305)
- [`dist_line_segment`](VEX/NearestDistance/neardistance.h#L344)
- [`dist_line_triangle`](VEX/NearestDistance/neardistance.h#L421)
- [`dist_pt_triangle`](VEX/NearestDistance/neardistance.h#L520)
- [`dist_segment_triangle`](VEX/NearestDistance/neardistance.h#L972)
- [`dist_triangle_triangle`](VEX/NearestDistance/neardistance.h#L1031)


### [Projection](VEX/Projection)
Projection simply involves finding the minimal distance.
This can be achieved using [Near Distance Functions](VEX/NearestDistance) (Geometric Tools Functions).

However, I found another implementations that might also be useful.
These functions project a point onto: [`line`](VEX/Projection/projection_functions.h#L90), [`segment`](VEX/Projection/projection_functions.h#L68), [`plane`](VEX/Projection/projection_functions.h#L53), [`triangle`](VEX/Projection/projection_functions.h#L121), [`planar`](VEX/Projection/projection_functions.h#L187) and [`non-planar quad`](VEX/Projection/projection_functions.h#L358).
<br><br>

### ![Hedges, Winding](VEX/hedge_winding_functions.h), and ![Edges](VEX/edge_functions.h)
- Snippets for ![iterating over half-edges (hedges)](VEX/hedge_winding_functions.h#L22).
- Finding the ![winding number](VEX/hedge_winding_functions.h#L182) (determine if a point is inside or outside 2D/3D geometry).
- ![Comparing winding](VEX/hedge_winding_functions.h#L110) direction between two primitives.
- ![Edge-related functions](VEX/edge_functions.h):
  - `primedges()` – Get edges of a primitive.
  - `pointedges()` – Get edges connected to a point.
  - `isedgeinprim()` – Check if an edge belongs to a specific primitive.
  - `edgeprims()` – Find primitives that share an edge.
  - `nearedgestopoint()` – Find the nearest edge to a point.
 
<img src="readme_images/winding2d.jpg" width=400px> <br>



### ![Lines, Curves, and UVWs](VEX/lines_functions.h)
- ![`resamplebylength()`](VEX/resample_functions.h#L101) – Resamples a curve (similar to the Resample SOP but as a function).
- ![`resamplebylength_pin()`](VEX/resample_functions.h#L113) - Resample a curve while keeping certain points fixed (pinned).
- ![`floodfill() + polyline_neighbours()`](VEX/lines_functions.h#L114) - Finds **polyline connectivity** based on constraint points, assigning different classes to segments at the Detail level.
- Finding points along a curve between UV coordinates.
- Finding nearby vertices based on UV coordinates.
- Finding a primitive line by two points.

# HIPs

## Solvecurve() / Solveik()
- ![Example file](Examples/solvecurve_solveik.hiplc) demonstrating the use of `solvecurve()` and `solveik()` functions in VEX.
<br><br> <img src="readme_images/solvecurve_solveik.gif" width=600px> <br>

<br><br>


# Nodes

## Edge Preservation
Houdini’s basic nodes don't always preserve edge groups well during geometry operations. This collection of nodes is designed to handle geometry while maintaining edge groups.

These nodes were developed during the creation of my **custom Boolean fracture system**, which can be found in the ![Shatterizer](https://github.com/alt-shiftov/Shatterizer) repo.

- <img src="readme_images/icons/grouptransfer.svg" width=20px> ![**Transfer Edge**](Nodes/sop_ALT.transferedge.1.0.hdalc): More accurate edge group transfer compared to Transfer Groups. <br>
- <img src="readme_images/icons/fuse.svg" width=20px> ![**Fuse Preserve Edges**](Nodes/sop_ALT.fuse_preserveedges.1.0.hdalc): Fuses points while keeping edge groups intact. <br>
- <img src="readme_images/icons/remesh.svg" width=20px> ![**Remesh Preserve Edges**](Nodes/sop_ALT.remesh_preserveedges.1.0.hdalc): Remeshes geometry while transferring and preserving edge groups. <br>
- <img src="readme_images/icons/groupcopy.svg" width=20px> ![**Copy Edge Groups**](Nodes/sop_ALT.copyedgegroup.1.0.hdalc): Copies edge groups, with an option to copy by ID. <br>


## Intersection and Shattering
These nodes were created to better understand intersection operations in Houdini.

- <img src="readme_images/icons/intersectionanalysis.svg" width=20px> ![**Intersection Analysis Custom:**](Nodes/sop_ALT.intersectionanalysis_custom.1.0.hdalc) A custom version of the **Intersection Analysis node**, which can be slightly faster than the original. Like the original node, it can store `@sourceprim`, `@sourceinput`, and `@sourceprimuv` attributes. 
Additionally, it can **split lines** based on the primitives they belong to and store information about which **edges** points are associated with.
<br> <img src="readme_images/intersectionanalysis.jpg" width=500px> <br>


- <img src="readme_images/icons/boolean.svg" width=20px> ![**Surface Shatter**](Nodes/sop_ALT.surfaceshatter.1.0.hdalc): A basic version of a **Boolean shatter** node using Surface Mode. This node is a work in progress, and I'm continuing to improve it. It was developed as part of my exploration into creating something similar to a Boolean Shatter operation.



## Handling Pieces

- <img src="readme_images/icons/explodedview.svg" width=20px> ![**Separate Pieces**](Nodes/sop_ALT.separatepieces.1.0.hdalc): This node is essential for **processing multiple pieces**. It separates pieces from each other using their piece attribute, allowing you to iterate through them simultaneously. This is useful for tasks like **generating captures** between identical pieces, **transferring attributes**, or **fusing by pieces**. It is much faster than using a foreach loop or piece groups. This method is based on the "Point Deform" node.
<br><img src="readme_images/separate_pieces.jpg" width=600px><br>

- <img src="readme_images/icons/bound.svg" width=20px> ![**Bound Pieces**](Nodes/sop_ALT.boundpiece.1.0.hdalc) - generate bound for pieces
- <img src="readme_images/icons/matchsize.svg" width=20px> ![**Match Size Piece**](Nodes/sop_ALT.matchsize_piece.1.0.hdalc) - make Match Size for pieces


## Useful for simulations

- <img src="readme_images/icons/smooth.svg" width=20px> ![**Dejitter**](Nodes/sop_ALT.dejiiter.1.0.hdalc): Removes jittering from animated points over time, especially useful for handling jittery simulations like cloth.
<br> <img src="readme_images/dejitter_vis.gif" width=600px> <br>

- <img src="readme_images/icons/fuse.svg" width=20px> **Fuse Generate** / **Fuse Extract**: A pair of nodes for fusing points and extracting values from the fused points.
	- ![**Fuse Generate**](Nodes/sop_ALT.fusegenerate.1.0.hdalc): Fuses points and generates attributes for use with the Fuse Extract node.
	- ![**Fuse Extract**](Nodes/sop_ALT.fuseextract.1.0.hdalc): Extracts attributes from fused geometry and applies them to unfused geometry. These nodes are useful for simulations where geometry needs to be fused, such as simulating cloth or wires. After the simulation, the animated attributes are transferred back to the original unfused geometry.
<br> <img src="readme_images/fuse_generate_fuse_extract.jpg" width=400px>
## Deformation and Skeleton Export

- ![**RBD to FBX**](Nodes/sop_ALT.rbd_to_fbx.1.0.hdalc): Converts rigid body animations into skeletal animations, allowing export in FBX format. This is useful for transferring RBD simulations into skeletal animations for game engines.

- ![**Deform to FBX**](Nodes/sop_ALT.deform_to_fbx.1.0.hdalc): Converts deforming geometry into skeletal animations, exportable in FBX format. Ideal for simulations such as cloth or soft bodies, with support for tearing cloth from Vellum. It can handle multiple pieces and offers four capture methods.

- ![**Wire to FBX**](Nodes/sop_ALT.wire_to_fbx.1.0.hdalc): Converts wire simulations from Vellum into skeletal animations for export in FBX format. Useful for wires, hair, string, and similar simulations. This node supports both intact and broken wires and works with multiple separate wires. Offers three capture methods.

- <img src="readme_images/icons/bonedeform.svg" width=20px> ![**pCapt to boneCapt**](Nodes/sop_ALT.pCapt_to_boneCapture.1.1.hdalc): Converts point capture data to bone capture data (`@pCaptPts` and `@pCaptWeights` to `@boneCapture`). This node allows to generate capture data between geometry and a skeleton using Point Deform, then convert it to **bone capture deformation**. Result will be work with the **Bone Deform** node.
<br><img src="readme_images/pointCapture to Bone Capture.jpg" width=450px><br>

- ![**Spline Deformer**](Nodes/sop_ALT.splinedeformer.1.0.hdalc): Deforms geometry based on a spline, capable of handling multiple pieces. It uses **primitive UVW** and **point deformation** and has an OpenCL implementation for better performance.
<br> <img src="readme_images/spline_deformer.gif" width=600px> <br>


## Additional Nodes

- <img src="readme_images/icons/uvtransfer.svg" width=20px> ![**UV Transfer**](Nodes/sop_ALT.uvtransfer.1.0.hdalc): Transfers UVs from one geometry to another. It's similar to the **Labs UV Transfer** node but much, much faster.
<br><br> <img src="readme_images/uvtransfer.jpg" width=600px> <br>
<img src="readme_images/uvtransfer_vs_labs.jpg" width=600px> <br>

- <img src="readme_images/icons/voronoifracture.svg" width=20px> ![**RBD Guided Voronoi**](Nodes/sop_ALT.rbdguidedvoronoi.1.0.hdalc): This node enhances Voronoi fracturing by using curves or planes as guides. It was developed during the search for a more precise solution to fracture geometry. While Voronoi fractures are accurate and stable, they often lack flexibility. This node makes the Voronoi fracture process more controllable.
<br><br> <img src="readme_images/rbd_guided_voronoi.jpg" width=600px> <br>


- <img src="readme_images/icons/view_vertex_order.svg" width=20px> ![**View Hedges**](Nodes/sop_ALT.view_hedges.1.0.hdalc): A node for visualizing half-edges (hedges). It can generate visualizer markers to display hedge numbers.
<br><br> <img src="readme_images/view_hedges.jpg" width=600px> <br>


