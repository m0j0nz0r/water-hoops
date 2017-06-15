using UnityEngine;
using System.Collections;

public class Tubulence : MonoBehaviour {
	public float turbulence;
	public float torqueFactor;
	private ConstantForce force;
	private Random r;
	// Use this for initialization
	void Start () {
		force = this.GetComponent<ConstantForce> ();
		r = new Random ();
	}
	
	// Update is called once per frame
	void Update () {
		float t = turbulence;
		float tf = t * torqueFactor;
		force.force = new Vector3 ((Random.value - .5f) * t, (Random.value - .5f)  * t, (Random.value - .5f)  * t);
		force.torque = new Vector3 ((Random.value - .5f)  * tf, (Random.value - .5f)  * tf, (Random.value - .5f)  * tf);
	}
}
