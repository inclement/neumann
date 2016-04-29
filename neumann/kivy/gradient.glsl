void main(void)
{
   float x = gl_FragCoord.x;
   float y = gl_FragCoord.y;
   gl_FragColor = vec4( x / resolution.x, y / resolution.y, 0.0, 1.0);
}