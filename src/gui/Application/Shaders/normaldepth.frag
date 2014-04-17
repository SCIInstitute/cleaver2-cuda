varying vec3 v_normal, v_position;
    uniform mat3 view_rot;
    uniform mat4 view;
    uniform float near, far;

    vec2 encode_normal(vec3 normal)
    {
        float f = sqrt(8.0*normal.z+8.0);
        return normal.xy / f + 0.5;
    }

    vec2 encode_depth(vec3 position){
        float depth = (length(position)-near)/far;
        depth = depth*255.0*255.0;
        return vec2(
            mod(depth, 255.0)/255.0,
            floor(depth/255.0)/255.0
        );
    }

    void main(void){
        vec3 normal = view_rot * normalize(v_normal);
        vec3 position = (view * vec4(v_position, 1.0)).xyz;
        gl_FragColor = vec4(encode_normal(normal), encode_depth(position));
    }
