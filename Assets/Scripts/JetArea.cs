using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class JetArea : MonoBehaviour {
	public float jetSpeed = 5f;
	public float torqueModifier = 0.2f;

	void OnTriggerStay(Collider collider){
		if (Input.GetKey (KeyCode.Space)) {
			collider.attachedRigidbody.AddForce (Vector3.up*jetSpeed, ForceMode.Acceleration);
		};
	}
}
