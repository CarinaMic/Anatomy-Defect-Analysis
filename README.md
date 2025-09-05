# Anatomy-Defect-Analysis

The Anatomy Defect Analysis toolbox integrates several specialized geometry-analysis classes – such as curvature, boundary, transformation, volume generating, intersection & distribution, and approximation – to provide an end-to-end solution for anatomical defect assessment.  
It contains a main script that structures the entire workflow: import, preprocess, create anatomical volume models, approximate, compare, and visualize anatomical volume models.

---

## Key Features
1. **Modular Integration:** Combines Curvature, Boundary, Transformation, IntersectionDistribution and Approximation classes under one unified interface  
2. **Main Script Management:** A central `main.m` controls the execution of the full anatomical analysis pipeline  
3. **Comprehensive Workflow:** From mesh import and preprocessing to multi-method analysis and visualization  

---

## How to Use
1. **Structure:** Use `main.m` to run the complete defect-analysis workflow. Individual class methods can also be invoked independently in custom workflows.  
2. **Import Data:** Load anatomical STL models into MATLAB. Ensure supporting functions (e.g., *inpolyhedron*) are on the path.  
3. **Execute Workflow:** Run `main.m` to perform sequential steps:  
   - Curvature analysis  
   - Boundary and transformation for alignment  
   - Volume generating  
   - Intersection and distribution analysis  
   - Approximation of volume models (multi-sphere models)  
4. **Visualize Results:** Output includes voxel-based comparisons, mesh overlays, and colored defect maps for intuitive interpretation.  

---

## Use Cases
- Orthopaedic research on bone deformities  
- Comparative anatomical modeling across patient data sets  
- Workflow prototyping for geometry processing in biomedical engineering  

---

## Benefits
- **All-in-One Framework:** Project unites multiple classes for anatomical defect analysis  
- **Reusability & Flexibility:** Classes remain modular and can be reused individually or collectively  
- **Proven Methodologies:** Use of already validated analytical tools to ensure reliability  

---

## Keywords
`anatomy`, `defect analysis`, `pelvis`, `orthopaedics`, `biomechanics`, `mesh processing`, `STL`, `triangulation`, `curvature`, `bounding box`, `scaling`, `transformation`, `mesh repair`, `hole filling`, `shrink wrap`, `alpha shape`, `intersection`, `sphere approximation`, `clumped spheres`, `volume analysis`, `rigid body alignment`, `Kabsch algorithm`, `distance mapping`, `boundary detection`, `point cloud`, `DBSCAN`

---

## Compatibility
These class is compatible with MATLAB and can be easily integrated into existing workflows to generate multi-sphere models from anatomical models.  
It is optimised for **MATLAB 2024a** but should be compatible with the most recent versions of MATLAB.

---

## Disclaimer
This MATLAB class is provided on the MATLAB File Exchange for educational and research purposes.  
Users should ensure that the class meets their specific analysis requirements and may need to adapt it accordingly.  
The code is provided *"as-is,"* and the author assumes no responsibility for its use or any consequences thereof.

---

## References
- M. P. Do Carmo, *Differential Geometry of Curves and Surfaces*, 2nd ed. A K Peters/CRC Press, Mar 2010. [Online].  
- M. Meyer, M. Desbrun, P. Schröder, and A. H. Barr, *‘‘Discrete Differential-Geometry Operators for Triangulated 2-Manifolds BT - Visualization and Mathematics III,’’* Visualization and Mathematics III, pp. 35–57, 2003. [Online].  
- K. Moreland, *‘‘Why we use bad color maps and what you can do about it,’’* Human Vision and Electronic Imaging 2016, HVEI 2016, pp. 262–267, 2016.  
- K. Subburaj, B. Ravi, and M. G. Agarwal, *‘‘3D shape reasoning for identifying anatomical landmarks,’’* Ph.D. dissertation, 2008.  
- K. Subburaj, B. Ravi, and M. Agarwal, *‘‘Automated identification of anatomical landmarks on 3D bone models reconstructed from CT scan images,’’* Computerized Medical Imaging and Graphics, vol. 33, no. 5, pp. 359–368, 2009.  
- W. Wilke, *‘‘Segmentierung und Approximation großer Punktwolken,’’* Dissertation, Technische Universität Darmstadt, 2002. [Online].  
- Y. Zhang, M. Fjeld, M. Fratarcangeli, A. Said, and S. Zhao, *‘‘Affective colormap design for accurate visual comprehension in industrial tomography,’’* Sensors, vol. 21, no. 14, pp. 1–20, 2021.  
- Johannes Korsawe (2024). [Minimal Bounding Box](https://www.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box), MATLAB Central File Exchange.  
- Svenja Reimer (2024). [calc_OriBoundingBox(data)](https://www.mathworks.com/matlabcentral/fileexchange/64417-calc_oriboundingbox-data), MATLAB Central File Exchange.  
- Anton Semechko (2024). [Exact minimum bounding spheres and circles](https://github.com/AntonSemechko/Bounding-Spheres-And-Circles), GitHub.  
- C. M. Micheler, J. J. Lang, N. J. Wilhelm, I. Lazic, F. Hinterwimmer, C. Fritz, R. V. Eisenhart-Rothe, M. F. Zäh, and R. H. Burgkart, *‘‘Scaling Methods of the Pelvis without Distortion for the Analysis of Bone Defects,’’* Current Directions in Biomedical Engineering, vol. 8, no. 2, pp.797–800, 2022.  
- W. Kabsch, *‘‘A solution for the best rotation to relate two sets of vectors,’’* Acta Crystallographica Section A, vol. 32, no. 5, pp. 922–923, Sep 1976. [Online].  
- W. Kabsch, *‘‘A discussion of the solution for the best rotation to relate two sets of vectors,’’* Acta Crystallographica Section A, vol. 34, no. 5, pp.827–828, Sep 1978. [Online].  
- E. Schreiber, [Kabsch algorithm](https://de.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm), MATLAB Central File Exchange, 2025. [Online].  
- N. Douillet, [Discrete contour mesh patch](https://de.mathworks.com/matlabcentral/fileexchange/78901-discrete-contour-mesh-patch-2d-3d), MATLAB Central File Exchange, 2025. [Online].  
- Sven, [inpolyhedron - Are points inside a triangulated volume?](https://de.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume), MATLAB Central File Exchange, 2024. [Online].  
- S. de Wolski, *‘‘File Exchange Pick of the Week - INPOLYHEDRON,’’* 2013. [Online]. Available: [https://blogs.mathworks.com/pick/2013/09/06/inpolyhedron/](https://blogs.mathworks.com/pick/2013/09/06/inpolyhedron/)  
- S. Haeri, *‘‘Optimisation of blade type spreaders for powder bed preparation in Additive Manufacturing using DEM simulations,’’* Powder Technology, vol. 321, pp. 94–104, 2017. [Online]. Available: [http://dx.doi.org/10.1016/j.powtec.2017.08.011](http://dx.doi.org/10.1016/j.powtec.2017.08.011)  
- S. Haeri, [sihaeri/DEM-ClumpedSphere](https://de.mathworks.com/matlabcentral/fileexchange/67754-sihaeri-dem-clumpedsphere), MATLAB Central File Exchange, 2025. [Online].  

---

## Cite As
Carina Micheler (2025). [Anatomy Defect Analysis](https://www.mathworks.com/matlabcentral/fileexchange/181346-anatomy-defect-analysis), MATLAB Central File Exchange.
