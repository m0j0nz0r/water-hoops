using UnityEngine;
using System.Collections;

public class Tubulence : MonoBehaviour {
	public float turbulence;
	public float torqueFactor;
	private ConstantForce force;
	private FluidSimulator simulator;
	// Use this for initialization
	void Start () {
		force = this.GetComponent<ConstantForce> ();
		simulator = FindObjectOfType<FluidSimulator> ();
	}
	
	// Update is called once per frame
	void Update () {
		float t = turbulence;
		float tf = t * torqueFactor;
		force.relativeForce = new Vector3 ((Random.value - .5f) * t, (Random.value - .5f)  * t, (Random.value - .5f)  * t);
		force.relativeTorque = new Vector3 ((Random.value - .5f)  * tf, (Random.value - .5f)  * tf, (Random.value - .5f)  * tf);
		force.force = simulator.getForce (transform.position);
	}
}
