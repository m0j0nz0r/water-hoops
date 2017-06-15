using UnityEngine;
using System.Collections;

public class Gamecontroller : MonoBehaviour {
	public Transform ring;
	public int ringCount = 15;
	public float xMin;
	public float xMax;
	public float yMin;
	public float yMax;
	// Use this for initialization
	void Start () {
	}
	
	// Update is called once per frame
	void Update () {
		if (ringCount > 0 && Random.value < 0.025f) {
			ringCount--;
			Instantiate (ring, new Vector3 (Random.Range (xMin, xMax), Random.Range (yMin, yMax), 0), Random.rotation);
		}
	
	}
}
