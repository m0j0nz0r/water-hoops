using UnityEngine;
using System.Collections;

public class ColorRandomizer : MonoBehaviour {
	// Use this for initialization
	public float colorDifference = 0.2f;
	void Start () {
		bool colored = false;
		float r = 0f, g = 0f, b = 0f;

		while (!colored) {
			if (Mathf.Abs (r - g) < colorDifference && Mathf.Abs (r - b) < colorDifference && Mathf.Abs (g - b) < colorDifference) {
				r = Random.value;
				g = Random.value;
				b = Random.value;
			} else {
				colored = true;
			}
		}
		Color color = new Color (r, g, b);
		Light light = GetComponent<Light> ();
		Renderer renderer = GetComponent<Renderer> ();

		light.color = color;
		renderer.material.color = color;

	}
	
	// Update is called once per frame
	void Update () {
	
	}
}
