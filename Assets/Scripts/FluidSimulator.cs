using System.Collections;
using System.Collections.Generic;
using UnityEngine;
public class FluidSimulator : MonoBehaviour {
	public int gridSizeX;
	public int gridSizeY;
	public int gridSizeZ;
	public float cellSize;
	public int linearSolverTimes = 20;
	public float diff;
	public float viscosity;

	private float[,,] density;
	private float[,,] densityPrev;
	public Vector3 relativeOrigin;
	private Vector3[,,] forceGrid;
	private Vector3[,,] forceGridPrev;

	void vectorAddSource (Vector3[,,] x, Vector3[,,] s, float deltaTime){
		int i, j, k, iMax = gridSizeX + 2, jMax = gridSizeY + 2, kMax = gridSizeZ + 2;
		for (i = 0; i <	 iMax; i++) {
			for (j = 0; j < jMax; j++) {
				for (k = 0; k < kMax; k++) {
					x [i, j, k] += deltaTime * s [i, j, k];
				}
			}
		}
	}
	void addSource( float[,,] x, float[,,] s, float dt){
		int i, j, k, iMax = gridSizeX + 2, jMax = gridSizeY + 2, kMax = gridSizeZ + 2;
		for (i = 0; i <	 iMax; i++) {
			for (j = 0; j < jMax; j++) {
				for (k = 0; k < kMax; k++) {
					x [i, j, k] += dt * s [i, j, k];
				}
			}
		}
	}

	void vectorDiffuse(Vector3[,,] x, Vector3[,,] x0, float diff, float deltaTime){
		float a = deltaTime * diff * gridSizeX * gridSizeY * gridSizeZ;
		vectorLinearSolve (x, x0, a, (1 + 6 * a), true);
	}
	void diffuse(int b, float[,,] x, float[,,] x0, float diff, float dt){
		float a = dt * diff * gridSizeX * gridSizeY * gridSizeZ;
		linearSolve (b, x, x0, a, (1 + 6 * a));
	}

	void vectorSetBoundary(Vector3[,,] x, bool bounce){
		int i, j, n = bounce ? -1:1;
		for (i = 1; i <= gridSizeY; i++) {
			for (j = 1; j <= gridSizeZ; j++) {
				x [1, i, j].x = -x [1, i, j].x*n;
				x [gridSizeX + 1, i, j].x = -x [gridSizeX, i, j].x*n;

			}
		}
		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeZ; j++) {
				x [i, 1, j].y = -x [i, 1, j].y*n;
				x [i, gridSizeY+1, j].y = -x [i, gridSizeY, j].y*n;
			}
		}

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				x [i, j, 1].z = -x [i, j, 1].z;
				x [i, j, gridSizeZ + 1].z = -x [i, j, gridSizeZ].z*n;
			}
		}

		x [0, 0, 0] = (
			x [1, 0, 0] + 
			x [0, 1, 0] + 
			x [0, 0, 1]
		) / 3;

		x [0, 0, gridSizeZ + 1] = (
			x [1, 0, gridSizeZ + 1] + 
			x [0, 1, gridSizeZ + 1] + 
			x [0, 0, gridSizeZ]
		) / 3;

		x [0, gridSizeY + 1, 0] = (
			x [1, gridSizeY + 1, 0] + 
			x [0, gridSizeY, 0] + 
			x [0, gridSizeY + 1, 1]
		) / 3;

		x [gridSizeX + 1, 0, 0] = (
			x [gridSizeX, 0, 0] + 
			x [gridSizeX + 1, 1, 0] + 
			x [gridSizeX + 1, 0, 1]
		) / 3;

		x [gridSizeX + 1, gridSizeY + 1, 0] = (
			x [gridSizeX + 1, gridSizeY + 1, 0] + 
			x [gridSizeX + 1, gridSizeY + 1, 1] + 
			x [gridSizeX + 1, gridSizeY, 0]
		) / 3;

		x [gridSizeX + 1, 0, gridSizeZ + 1] = (
			x [gridSizeX, 0, gridSizeZ + 1] + 
			x [gridSizeX + 1, 1, gridSizeZ + 1] + 
			x [gridSizeX + 1, 0, gridSizeZ]
		) / 3;

		x [0, gridSizeY + 1, gridSizeZ + 1] = (
			x [1, gridSizeY + 1, gridSizeZ + 1] + 
			x [0, gridSizeY, gridSizeZ + 1] + 
			x [0, gridSizeY + 1, gridSizeZ]
		) / 3;

		x [gridSizeX + 1, gridSizeY + 1, gridSizeZ + 1] = (
			x [gridSizeX, gridSizeY + 1, gridSizeZ + 1] + 
			x [gridSizeX + 1, gridSizeY, gridSizeZ + 1] + 
			x [gridSizeX + 1, gridSizeY + 1, gridSizeZ]) / 3;
	}
	void setBoundary(int b, float[,,] x){
		int i, j;
		for (i = 1; i <= gridSizeY; i++) {
			for (j = 1; j <= gridSizeZ; j++) {
				x [1, i, j] = b == 1 ? x [1, i, j] : -x [1, i, j];
				x [gridSizeX + 1, i, j] = b == 1 ? x [gridSizeX, i, j] : -x [gridSizeX, i, j];

			}
		}

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeZ; j++) {
				x [i, 1, j] = b == 2 ? x [i, 1, j] : -x [i, 1, j];
				x [i, gridSizeY + 1, j] = b == 2 ? x [i, gridSizeY, j] : -x [i, gridSizeY, j];
			}
		}

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				x [i, j, 1] = b == 3 ? x [i, j, 1] : -x [i, j, 1];
				x [i, j, gridSizeZ + 1] = b == 3 ? x [i, j, gridSizeZ] : -x [i, j, gridSizeZ];
			}
		}

		x [0, 0, 0] = (
			x [1, 0, 0] + 
			x [0, 1, 0] + 
			x [0, 0, 1]
		) / 3;

		x [0, 0, gridSizeZ + 1] = (
			x [1, 0, gridSizeZ + 1] + 
			x [0, 1, gridSizeZ + 1] + 
			x [0, 0, gridSizeZ]
		) / 3;
		
		x [0, gridSizeY + 1, 0] = (
			x [1, gridSizeY + 1, 0] + 
			x [0, gridSizeY, 0] + 
			x [0, gridSizeY + 1, 1]
		) / 3;
		
		x [gridSizeX + 1, 0, 0] = (
			x [gridSizeX, 0, 0] + 
			x [gridSizeX + 1, 1, 0] + 
			x [gridSizeX + 1, 0, 1]
		) / 3;

		x [gridSizeX + 1, gridSizeY + 1, 0] = (
			x [gridSizeX + 1, gridSizeY + 1, 0] + 
			x [gridSizeX + 1, gridSizeY + 1, 1] + 
			x [gridSizeX + 1, gridSizeY, 0]
		) / 3;

		x [gridSizeX + 1, 0, gridSizeZ + 1] = (
			x [gridSizeX, 0, gridSizeZ + 1] + 
			x [gridSizeX + 1, 1, gridSizeZ + 1] + 
			x [gridSizeX + 1, 0, gridSizeZ]
		) / 3;

		x [0, gridSizeY + 1, gridSizeZ + 1] = (
			x [1, gridSizeY + 1, gridSizeZ + 1] + 
			x [0, gridSizeY, gridSizeZ + 1] + 
			x [0, gridSizeY + 1, gridSizeZ]
		) / 3;

		x [gridSizeX + 1, gridSizeY + 1, gridSizeZ + 1] = (
			x [gridSizeX, gridSizeY + 1, gridSizeZ + 1] + 
			x [gridSizeX + 1, gridSizeY, gridSizeZ + 1] + 
			x [gridSizeX + 1, gridSizeY + 1, gridSizeZ]) / 3;

	}

	void vectorAdvect(Vector3[,,] forceGrid, Vector3[,,] forceGridPrev, float deltaTime){
		int i = gridSizeX + 1, j = gridSizeY + 1, k = gridSizeZ + 1, i0, j0, k0, i1, j1, k1;
		float x, y, z, r0, s0, t0, r1, s1, t1;

		while (i-->1) {
			while (j-->1) {
				while (k-->1) {

					//newPos = oldPos - time*velocity
					x = Mathf.Clamp(i - deltaTime*forceGrid [i, j, k].x, 0, gridSizeX + 1);
					y = Mathf.Clamp(j - deltaTime*forceGrid [i, j, k].y, 0, gridSizeY + 1);
					z = Mathf.Clamp(k - deltaTime*forceGrid [i, j, k].z, 0, gridSizeZ + 1);

					//each newPos cell value is a result of a proportional addition of each of its forward neighbors.
					i0 = Mathf.RoundToInt (x);
					i1 = i0 - 1;
					j0 = Mathf.RoundToInt (y);
					j1 = j0 - 1;
					k0 = Mathf.RoundToInt (z);
					k1 = k0 - 1;

					s1 = x - i0;
					s0 = 1 - s1;
					t1 = y - j0;
					t0 = 1 - t1;
					r1 = z - k0;
					r0 = 1 - r1;

					forceGrid [i, j, k] = s0 * 
						(
							t0 * (
								r0 * forceGridPrev [i0, j0, k0] + 
								r1 * forceGridPrev [i0, j0, k1]
							) + 
							t1 * (
								r0 * forceGridPrev [i0, j1, k0] + 
								r1 * forceGridPrev [i0, j1, k1]
							)
						) + 
						s1 * (
							t0 * (
								r0 * forceGridPrev [i1, j0, k0] + 
								r1 * forceGridPrev [i1, j0, k1]
							) + 
							t1 * (
								r0 * forceGridPrev [i1, j1, k0] + 
								r1 * forceGridPrev [i1, j1, k1]
							)
						);
				}
			}
		}

		vectorSetBoundary (forceGrid, true);
	}
	void advect( int b, float[,,] d, float[,,] d0, float[,,] u, float[,,] v, float[,,] w, float dt){
		int i = gridSizeX + 1, j = gridSizeY + 1, k = gridSizeZ + 1, i0, j0, k0, i1, j1, k1;
		float x, y, z, r0, s0, t0, r1, s1, t1;

		while (i-->1) {
			while (j-->1) {
				while (k-->1) {

					//newPos = oldPos - time*velocity
					x = Mathf.Clamp(i - dt*u [i, j, k], 0, gridSizeX + 1);
					y = Mathf.Clamp(j - dt*v [i, j, k], 0, gridSizeY + 1);
					z = Mathf.Clamp(k - dt*w [i, j, k], 0, gridSizeZ + 1);

					//each newPos cell value is a result of a proportional addition of each of its forward neighbors.
					i0 = Mathf.RoundToInt (x);
					i1 = i0 - 1;
					j0 = Mathf.RoundToInt (y);
					j1 = j0 - 1;
					k0 = Mathf.RoundToInt (z);
					k1 = k0 - 1;

					s1 = x - i0;
					s0 = 1 - s1;
					t1 = y - j0;
					t0 = 1 - t1;
					r1 = z - k0;
					r0 = 1 - r1;

					d [i, j, k] = s0 * 
						(
							t0 * (
								r0 * d0 [i0, j0, k0] + 
								r1 * d0 [i0, j0, k1]
							) + 
							t1 * (
								r0 * d0 [i0, j1, k0] + 
								r1 * d0 [i0, j1, k1]
							)
						) + 
						s1 * (
							t0 * (
								r0 * d0 [i1, j0, k0] + 
								r1 * d0 [i1, j0, k1]
							) + 
							t1 * (
								r0 * d0 [i1, j1, k0] + 
								r1 * d0 [i1, j1, k1]
							)
						);
				}
			}
		}

		setBoundary (b, d);
	}

	void densityStep(float[,,] x, float[,,] x0, float[,,]u, float[,,] v, float[,,] w, float diff, float dt){
		addSource (x, x0, dt);
		Swap<float[,,]> (x, x0);
		diffuse (0, x, x0, diff, dt);
		Swap<float[,,]> (x, x0);
		advect (0, x, x0, u, v, w, dt);
	}

	void vectorVelocityStep(Vector3[,,] forceGrid, Vector3[,,] forceGridPrev, float viscosity, float deltaTime){
		vectorAddSource (forceGrid, forceGridPrev, deltaTime);
		vectorDiffuse (forceGridPrev, forceGrid, viscosity, deltaTime);
		vectorProject (forceGridPrev);
		vectorAdvect (forceGrid, forceGridPrev, deltaTime);
		vectorProject (forceGrid);
		Debug.Log (string.Format("a: {0}\nb: {1}", test (forceGrid), test (forceGridPrev)));
	}
	void velocityStep(float[,,]u, float[,,]v, float[,,]w, float[,,] u0, float[,,] v0, float[,,] w0, float viscosity, float dt){
		addSource (u, u0, dt);
		addSource (v, v0, dt);
		addSource (w, w0, dt);
		Swap<float[,,]> (u, u0);
		Swap<float[,,]> (v, v0);
		Swap<float[,,]> (w, w0);
		diffuse (1, u, u0, viscosity, dt);
		diffuse (2, v, v0, viscosity, dt);
		diffuse (3, w, w0, viscosity, dt);
		project (u, v, w, u0, v0);
		Swap<float[,,]> (u, u0);
		Swap<float[,,]> (v, v0);
		Swap<float[,,]> (w, w0);
		advect (1, u, u0, u0, v0, w0, dt);
		advect (2, v, v0, u0, v0, w0, dt);
		advect (3, w, w0, u0, v0, w0, dt);
		project (u, v, w, u0, v0);

	}

	void vectorProject(Vector3[,,] forceGrid){
		float[,,] p = getNewGrid<float> (), div = getNewGrid<float>();

		int i, j, k;
		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					div [i, j, k] = -(
						(forceGrid [i + 1, j, k].x - forceGrid [i - 1, j, k].x)/gridSizeX + 
						(forceGrid [i, j + 1, k].y - forceGrid [i, j - 1, k].y)/gridSizeY + 
						(forceGrid [i, j, k + 1].z - forceGrid [i, j, k - 1].z)/gridSizeZ
					)/3;
					p [i, j, k] = 0f;
				}
			}
		}

		setBoundary(0, div);
		setBoundary(0, p);

		linearSolve(0, p, div, 1, 6);

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					forceGrid [i, j, k] = forceGrid[i,j,k] - new Vector3 (
						gridSizeX * (p [i + 1, j, k] - p [i - 1, j, k])/2,
						gridSizeY * (p [i, j + 1, k] - p [i, j - 1, k])/2,
						gridSizeZ * (p [i, j, k + 1] - p [i, j, k - 1])/2
					);
				}
			}
		}

		vectorSetBoundary (forceGrid, true);

	}
	void project(float[,,] u, float[,,] v, float[,,] w, float[,,] p, float[,,] div){
		int i, j, k;
		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					div [i, j, k] = -(
						(u [i + 1, j, k] - u [i - 1, j, k])/gridSizeX + 
						(v [i, j + 1, k] - v [i, j - 1, k])/gridSizeY + 
						(w [i, j, k + 1] - w [i, j, k - 1])/gridSizeZ
					)/3;
					p [i, j, k] = 0;
				}
			}
		}

		setBoundary (0, div);
		setBoundary(0, p);

		linearSolve(0, p, div, 1, 6);

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					u [i, j, k] -= 0.5f * gridSizeX * (p [i + 1, j, k] - p [i - 1, j, k]);
					v [i, j, k] -= 0.5f * gridSizeY * (p [i, j + 1, k] - p [i, j - 1, k]);
					w [i, j, k] -= 0.5f * gridSizeZ * (p [i, j, k + 1] - p [i, j, k - 1]);
				}
			}
		}
		setBoundary (1, u);
		setBoundary (2, v);
		setBoundary (3, w);
	}
	void vectorLinearSolve(Vector3[,,] x, Vector3[,,] x0, float a, float c, bool bounce){
		int n, i, j, k;
		for (n = 0; n < linearSolverTimes; n++) {
			for (i = 1; i <= gridSizeX; i++) {
				for (j = 1; j <= gridSizeY; j++) {
					for (k = 1; k <= gridSizeZ; k++) {
						x[i,j,k] = (x0[i,j,k] + a*(
							x[i-1,j,k]+
							x[i+1,j,k]+
							x[i,j-1,k]+
							x[i,j+1,k]+
							x[i,j,k-1]+
							x[i,j,k+1]
						))/c;
					}
				}
			}
			vectorSetBoundary (x, bounce);
		}
	}
	void linearSolve(int b, float[,,] x, float[,,] x0, float a, float c){
		int n, i, j, k;
		for (n = 0; n < linearSolverTimes; n++) {
			for (i = 1; i <= gridSizeX; i++) {
				for (j = 1; j <= gridSizeY; j++) {
					for (k = 1; k <= gridSizeZ; k++) {
						x[i,j,k] = (x0[i,j,k] + a*(
							x[i-1,j,k]+
							x[i+1,j,k]+
							x[i,j-1,k]+
							x[i,j+1,k]+
							x[i,j,k-1]+
							x[i,j,k+1]
						))/c;
					}
				}
			}
			setBoundary (b, x);
		}
	}
	void Start(){
		Init ();
	}
	void Init(){
		density = getNewGrid<float> ();
		densityPrev = getNewGrid<float> ();
		forceGrid = getNewGrid<Vector3> ();
		forceGridPrev = getNewGrid<Vector3> ();
		//maxSize = Mathf.Max (gridSizeX, Mathf.Max (gridSizeY, gridSizeZ));
		relativeOrigin = transform.position - new Vector3 ((gridSizeX+2) * cellSize / 2, (gridSizeY+2) * cellSize / 2, (gridSizeZ+2) * cellSize / 2);
		InitGrid (density, 1f);
		InitGrid (densityPrev, 1f);
		InitVectorGrid (forceGrid, Vector3.zero);
		InitVectorGrid (forceGridPrev, Vector3.zero);
	}
	void InitGrid (float[,,] grid, float value){
		int i, j, k, iMax = grid.GetLength(0), jMax = grid.GetLength(1), kMax = grid.GetLength(2);
		for (i = 0; i < iMax; i++) {
			for (j = 0; j < jMax; j++) {
				for (k = 0; k < kMax; k++) {
					grid [i, j, k] = value;
				}
			}
		}
	}
	void InitVectorGrid (Vector3[,,] grid, Vector3 value){
		int i, j, k, iMax = grid.GetLength(0), jMax = grid.GetLength(1), kMax = grid.GetLength(2);
		for (i = 0; i < iMax; i++) {
			for (j = 0; j < jMax; j++) {
				for (k = 0; k < kMax; k++) {
					grid [i, j, k] = new Vector3(value.x, value.y, value.z);
				}
			}
		}
	}
	T[,,] getNewGrid<T>(){
		return new T[gridSizeX + 2, gridSizeY + 2, gridSizeZ + 2];
	}
	public void addForce(Vector3 position, Vector3 force){
		int x = getCell (relativeOrigin.x, position.x, gridSizeX + 1),
		y = getCell (relativeOrigin.y, position.y, gridSizeY + 1), 
		z = getCell (relativeOrigin.z, position.z, gridSizeZ + 1);
		forceGridPrev [x, y, z] = force;
	}
	public Vector3 getForce(Vector3 position){
		int x = getCell (relativeOrigin.x, position.x, gridSizeX + 1),
		y = getCell (relativeOrigin.y, position.y, gridSizeY + 1), 
		z = getCell (relativeOrigin.z, position.z, gridSizeZ + 1);
		return forceGrid[x,y,z];
	}
	int getCell(float origin, float position, int max){
		return Mathf.Clamp(Mathf.RoundToInt ((position - origin)/cellSize), 0, max);
	}
	void Swap<T>(T a, T b){
		T tmp = a;
		a = b;
		b = tmp;
	}
	void Update(){
		if (forceGrid == null || forceGridPrev == null) {
			Init ();
		}
		vectorVelocityStep (forceGrid, forceGridPrev, viscosity, Time.deltaTime);
	}
	void OnDrawGizmos(){
		int i, j, k;
		Vector3 source, dest, vec;
		for (i = 0; i<gridSizeX+2;i++){
			for (j = 0; j<gridSizeY+2;j++){
				for (k = 0; k<gridSizeZ+2;k++){
					if (forceGrid != null) {
						source = new Vector3 (i * cellSize, j * cellSize, k * cellSize) + relativeOrigin;
						vec = forceGrid [i, j, k];
						dest = source + vec/200;
						if (vec.magnitude > 0.1f) {
							//Gizmos.color = Color.blue;
							//Gizmos.DrawLine(source, dest);
							Gizmos.color = Color.red;
							Gizmos.DrawLine (source, source + forceGridPrev [i, j, k] / 200);
								
							//Gizmos.DrawSphere(dest, 0.05f);
						}
					}
				}
			}
		}
	}
	float test(Vector3[,,] a){
		float retval = 0f;
		foreach (Vector3 v in a) retval += v.magnitude;
		return retval;
	}
}
