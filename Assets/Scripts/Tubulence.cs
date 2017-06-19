using UnityEngine;
using System.Collections;

public class Tubulence : MonoBehaviour {
	public float turbulence;
	public float torqueFactor;
	private ConstantForce force;
	// Use this for initialization
	void Start () {
		force = this.GetComponent<ConstantForce> ();
	}
	
	// Update is called once per frame
	void Update () {
		float t = turbulence;
		float tf = t * torqueFactor;
		force.relativeForce = new Vector3 ((Random.value - .5f) * t, (Random.value - .5f)  * t, (Random.value - .5f)  * t);
		force.relativeTorque = new Vector3 ((Random.value - .5f)  * tf, (Random.value - .5f)  * tf, (Random.value - .5f)  * tf);
	}
}
