using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteAlways]
public class TimeProvider : MonoBehaviour
{
    [Range(0, 100)]
    public float animTime;
    private MaterialPropertyBlock _mpb;
    public MaterialPropertyBlock mpb
    {
        get
        {
            if(_mpb is null)
            {
                _mpb = new MaterialPropertyBlock();
            }
            return _mpb;
        }
    }

    private MeshRenderer _renderer;

    private void OnEnable()
    {
        _renderer = GetComponent<MeshRenderer>();
    }

    // Update is called once per frame
    void Update()
    {
        mpb.SetFloat("_shaderTime", Time.time);
        mpb.SetFloat("_animTime", animTime);
        _renderer.SetPropertyBlock(mpb);
    }
}
