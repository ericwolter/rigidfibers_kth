kernel void vadd(
    const global fiberfloat *a,
    const global fiberfloat *b,
    global fiberfloat *c,
    const fiberuint count)
{
    fiberuint i = get_global_id(0);

    if(i < count) 
    {
        c[i] = a[i] + b[i];
    }
}
