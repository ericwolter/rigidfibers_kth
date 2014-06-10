kernel void vadd(
    const global fiberfloat *a,
    const global fiberfloat *b,
    global fiberfloat *c,
    const size_t count)
{
    size_t i = get_global_id(0);

    if(i < count) 
    {
        c[i] = a[i] + b[i];
    }
}
