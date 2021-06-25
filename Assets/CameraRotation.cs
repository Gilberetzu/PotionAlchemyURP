using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraRotation : MonoBehaviour
{
    public float speed = 1;
    public float staticDuration = 2;
    public float startTime = 0;
    private Transform _transform;
    // Start is called before the first frame update
    void Start()
    {
        _transform = gameObject.transform;
        startTime = Time.time;
        _transform.rotation = Quaternion.Euler(0, -45, 0);
    }

    // Update is called once per frame
    void Update()
    {
        if(Time.time - startTime - staticDuration > 0 ){
            float offsetTime = Time.time - startTime - staticDuration;
            _transform.rotation = Quaternion.Euler(Mathf.Sin(offsetTime * speed) * 90, Mathf.Sin(offsetTime * speed/2) * 90 - 45, 0);
        }
    }
}
