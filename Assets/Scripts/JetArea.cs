using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class JetArea : MonoBehaviour {
	public KeyCode trigger;
	public Vector3 force;
	private FluidSimulator simulator;
	void Start(){
		simulator = FindObjectOfType<FluidSimulator> ();
	}
	void Update(){
		simulator.addForce (transform.position, Vector3.Lerp(force, force.magnitude*(Random.rotation.eulerAngles.normalized - Vector3.forward/2 - Vector3.right/2), 0.5f));
	}
}
