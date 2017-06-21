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
		if (Input.GetKeyDown (trigger)) {
			simulator.addForce (transform.position, force);
		}
		if (Input.GetKeyUp (trigger)) {
			simulator.addForce (transform.position, Vector3.zero);
		}
	}
}
