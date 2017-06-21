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

	private float[,,] u;
	private float[,,] v;
	private float[,,] w;
	private float[,,] uPrev;
	private float[,,] vPrev;
	private float[,,] wPrev;
	private float[,,] density;
	private float[,,] densityPrev;
	private int maxSize;
	public Vector3 relativeOrigin;

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

	void diffuse(int b, float[,,] x, float[,,] x0, float diff, float dt){
		float a = dt * diff * maxSize * maxSize * maxSize;
		linearSolve (b, x, x0, a, (1 + 6 * a));
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
				x [i, gridSizeY+1, j] = b == 2 ? x [i, gridSizeY, j] : -x [i, gridSizeY, j];
			}
		}

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				x [i, j, 1] = b == 3 ? x [i, j, 1] : -x [i, j, 1];
				x [i, j, gridSizeZ + 1] = b == 3 ? x [i, j, gridSizeZ] : -x [i, j, gridSizeZ];
			}
		}

		x [0, 0, 0] = (x [1, 0, 0] + x [0, 1, 0] + x [0, 0, 1]) / 3;

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

	void advect( int b, float[,,] d, float[,,] d0, float[,,] u, float[,,] v, float[,,] w, float dt){
		int i, j, k, i0, j0, k0, i1, j1, k1;
		float x, y, z, r0, s0, t0, r1, s1, t1, dt0 = dt*maxSize;

		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					x = Mathf.Clamp(i - dt0*u [i, j, k], 0.5f, gridSizeX + 0.5f);
					y = Mathf.Clamp(j - dt0*v [i, j, k], 0.5f, gridSizeX + 0.5f);
					z = Mathf.Clamp(k - dt0*w [i, j, k], 0.5f, gridSizeX + 0.5f);

					i0 = Mathf.RoundToInt (x);
					i1 = i0 + 1;
					j0 = Mathf.RoundToInt (y);
					j1 = j0 + 1;
					k0 = Mathf.RoundToInt (z);
					k1 = k0 + 1;

					s1 = x - i0;
					s0 = 1 - s1;
					t1 = y - j0;
					t0 = 1 - t1;
					r1 = z - k0;
					r0 = 1 - r1;

					d [i, j, k] = 
						s0 * (
							t0 * (
								r0*d0 [i0, j0, k0] + 
								r1*d0[i0, j0, k1]
							) + 
							t1 * (
								r0*d0 [i0, j1, k0] + 
								r1*d0[i0, j1, k1]
							)
						) + 
						s1 * (
							t0 * (
								r0*d0 [i1, j0, k0] + 
								r1*d0[i1, j0, k1]
							) + 
							t1 * (
								r0*d0 [i1, j1, k0] + 
								r1*d0[i1, j1, k1]
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

	void project(float[,,] u, float[,,] v, float[,,] w, float[,,] p, float[,,] div){
		int i, j, k;
		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					div [i, j, k] = -(
						(u [i + 1, j, k] - u [i - 1, j, k])/gridSizeX + 
						(v [i, j - 1, k] - v [i, j + 1, k])/gridSizeY + 
						(w [i, j, k - 1] + w [i, j, k + 1])/gridSizeZ
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
		u = getNewGrid ();
		v = getNewGrid ();
		w = getNewGrid ();
		uPrev = getNewGrid ();
		vPrev = getNewGrid ();
		wPrev = getNewGrid ();
		density = getNewGrid ();
		densityPrev = getNewGrid ();
		maxSize = Mathf.Max (gridSizeX, Mathf.Max (gridSizeY, gridSizeZ));
		relativeOrigin = transform.position - new Vector3 ((gridSizeX+2) * cellSize / 2, (gridSizeY+2) * cellSize / 2, (gridSizeZ+2) * cellSize / 2);
		InitGrid (u, 0f);
		InitGrid (v, 0f);
		InitGrid (w, 0f);
		InitGrid (uPrev, 0f);
		InitGrid (vPrev, 0f);
		InitGrid (wPrev, 0f);
		InitGrid (density, 1f);
		InitGrid (densityPrev, 1f);
	}
	void InitGrid (float[,,] grid, float value){
		int i, j, k;
		for (i = 1; i <= gridSizeX; i++) {
			for (j = 1; j <= gridSizeY; j++) {
				for (k = 1; k <= gridSizeZ; k++) {
					grid [i, j, k] = value;
				}
			}
		}
	}
	float[,,] getNewGrid(){
		return new float[gridSizeX + 2, gridSizeY + 2, gridSizeZ + 2];
	}
	public void addForce(Vector3 position, Vector3 force){
		Vector3 relativePosition = position - relativeOrigin;
		int x = Mathf.Clamp (Mathf.RoundToInt (relativePosition.x), 0, gridSizeX+1),
		y = Mathf.Clamp (Mathf.RoundToInt (relativePosition.y), 0, gridSizeY+1), 
		z = Mathf.Clamp (Mathf.RoundToInt (relativePosition.z), 0, gridSizeZ+1);
		u [x, y, z] = force.x;
		v [x, y, z] = force.y;
		w [x, y, z] = force.z;
	}
	public Vector3 getForce(Vector3 position){
		Vector3 relativePosition = position - relativeOrigin;
		int x = Mathf.Clamp (Mathf.RoundToInt (relativePosition.x), 0, gridSizeX+1), 
		y = Mathf.Clamp (Mathf.RoundToInt (relativePosition.y), 0, gridSizeY+1), 
		z = Mathf.Clamp (Mathf.RoundToInt (relativePosition.z), 0, gridSizeZ+1);

		return new Vector3 (u [x, y, z], v [x, y, z], w [x, y, z]);
	}
	void Swap<T>(T a, T b){
		T tmp = a;
		a = b;
		b = tmp;
	}
	void Update(){
		velocityStep (u, v, w, uPrev, vPrev, wPrev, viscosity, Time.deltaTime);
		densityStep (density, densityPrev, u, v, w, diff, Time.deltaTime);
	}
	void OnDrawGizmos(){
		int i, j, k;
		Vector3 source;
		for (i = 0; i<gridSizeX+2;i++){
			for (j = 0; j<gridSizeY+2;j++){
				for (k = 0; k<gridSizeZ+2;k++){
					source = new Vector3 (i * cellSize, j * cellSize, k * cellSize) + relativeOrigin;
					if (u != null && v != null && w != null) {
						Gizmos.color = Color.blue;
						Gizmos.DrawLine(source, source + new Vector3(u[i,j,k], v[i,j,k], w[i,j,k]).normalized);
					}
				}
			}
		}
	}
	void Check(float[,,] x){
		int count = 0;
		float sum = 0f;
		foreach (float f in x) {
			if (f > 0) {
				sum += f;
				count++;
			}
		}
		Debug.Log (string.Format("Count: {0}\nSum: {1}", count, sum));
	}
}
